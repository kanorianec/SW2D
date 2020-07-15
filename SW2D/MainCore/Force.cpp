#include "Raschet.h"
#include "technical.h"
#include <cmath>
#include "Constants.h"

void Raschet::SetWindSpeed(forceType fT, double WindFrictionCoefficient, double period, double xWindConst, double yWindConst)
{
	windForcing = fT;
	windFrictionCoef = WindFrictionCoefficient;

	switch (windForcing)
	{
	case NO_FORCE:
		break;
	case CONST_FORCE:
		xWind = new double[Nx*Ny]();
		yWind = new double[Nx*Ny]();
		for (int k = 0; k < Nx*Ny; k++)
		{
			xWind[k] = xWindConst;
			yWind[k] = yWindConst;
		}
		break;
	case REAL_FORCE:
		timeWind = T_begin - 1;
		timeWindPeriod = period;

		FWindX.open("windSpeed/xWind.dat", std::ios::binary);
		FWindY.open("windSpeed/yWind.dat", std::ios::binary);

		if (!FWindX.is_open() || !FWindY.is_open())
			std::cout << "Error: windSpeed/x(y)Wind was not found!" << endl;
		else
		{
			xWind = new double[Nx*Ny]();
			yWind = new double[Nx*Ny]();
		}
		break;
	default:
		cout << "Error, wrong force type parameter " << fT << "!" << endl;
		break;
	}
}


void Raschet::Recalc_forces_parallel()
{

	if (TideForcing)
	{
		// Block of tide force initialization =============================

		// convert time to convenient format
		struct tm* tt;
		time_t to_convert = dT + Time_elapsed + RaschetTime;

		to_convert += 0;// 6 * 3600 * 24;
		tt = gmtime(&to_convert);

		int day = tt->tm_mday;
		int month = tt->tm_mon + 1;
		int year = tt->tm_year + 1900;
		int hour = tt->tm_hour;
		int minute = tt->tm_min;
		int second = tt->tm_sec;

		//if (Time_elapsed >= HourMark * 3600)
		//{
		// Recalculate ephemeris for Sun and Moon
		//Ephemeris::setLocationOnEarth(lat0, lonN); // for ephemeris calculation

		SolarSystemObject pSun = Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
		SolarSystemObject pMoon = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);

		// designations in variables: S - Sun, M - Moon

		Sbeta = rad * pSun.equaCoordinates.dec;
		Mbeta = rad * pMoon.equaCoordinates.dec;

		Sa = rad * pSun.equaCoordinates.ra * 360 / 24;
		Ma = rad * pMoon.equaCoordinates.ra * 360 / 24;

		ST = rad * Ephemeris::apparentSideralTime(day, month, year, hour, minute, second) * 360 / 24;

		//cout << Ephemeris::apparentSideralTime(day, month, year, hour, minute, second) * 360 / 24 + 62 - pSun.equaCoordinates.ra * 360 / 24 << endl;
		//pause();

		SDist = pSun.distance;
		MDist = pMoon.distance;

		//	HourMark++;
		//}	

		// ==================================

		//double OmS, OmM;
		//OmS = SD/gc * (SunAxis / SDist) * (SunAxis / SDist) * (SunAxis / SDist) * (3*(sin(latN * rad)*sin(latN * rad) - 1/3)*(sin(Sbeta)*sin(Sbeta) - 1/3) +  sin(2 * latN * rad)*sin(2 * Sbeta)*(cos(ST - Sa + 2 * lonN * rad)) + cos(latN * rad)*cos(latN * rad)*cos(Sbeta)*cos(Sbeta)*(cos(2 * (ST - Sa + 2 * lonN * rad))));
		//OmM += MD/gc * (MoonAxis / MDist) * (MoonAxis / MDist) * (MoonAxis / MDist) * (3 * (sin(latN * rad)*sin(latN * rad) - 1 / 3)*(sin(Mbeta)*sin(Mbeta) - 1 / 3) + sin(2 * latN * rad)*sin(2 * Mbeta)*(cos(ST - Ma + 2 * lonN * rad)) + cos(latN * rad)*cos(latN * rad)*cos(Mbeta)*cos(Mbeta)*(cos(2 * (ST - Ma + 2 * lonN * rad))));

		//	cout << ST / rad << " " << Sa / rad << " " << 2 * lat0 << " " << (ST - Sa) / rad + 2 * lat0 << " " << asctime(gmtime(&to_convert)) << endl;
		//		cout << Sbeta/rad << " " << SD / gc * 3 * (sin(latN * rad)*sin(latN * rad) - 1 / 3)*(sin(Sbeta)*sin(Sbeta) - 1 / 3) << endl;
		//		cout << Mbeta / rad << " " << MD / gc * 3 * (sin(latN * rad)*sin(latN * rad) - 1 / 3)*(sin(Mbeta)*sin(Mbeta) - 1 / 3) << endl;
		//cout << SD / gc << " " << sin(2 * latN * rad) << " " <<sin(2 * Sbeta)/*(cos(ST - Sa + 2 * lonN * rad))*/ << endl;
		//cout << SD / gc /* (SunAxis / SDist) * (SunAxis / SDist) * (SunAxis / SDist)*/* cos(latN * rad) *cos(latN * rad)*cos(Sbeta)*cos(Sbeta)/*(cos(2 * (ST - Sa + 2 * lonN * rad)))*/ << endl;
		//cout << MD / gc /* (MoonAxis / MDist) * (MoonAxis / MDist) * (MoonAxis / MDist)*/ * cos(latN * rad)*cos(latN * rad)*cos(Mbeta)*cos(Mbeta) << endl;
		//pause();
	}

	if (windForcing == REAL_FORCE)
	{
		if (Time_elapsed > timeWind)
		{
			if (FWindX.eof() || FWindY.eof())
			{
				cout << "END OF FILE OF WIND VELOCITIES! Continue with last values" << endl;
			}
			else
			{
				FWindX.read(reinterpret_cast<char*> (xWind), sizeof(double) * Nx * Ny);
				FWindY.read(reinterpret_cast<char*> (yWind), sizeof(double) * Nx * Ny);
			}
			timeWind += timeWindPeriod;
		}
	}

#pragma omp parallel for
	for (int k = 0; k < Nx*Ny; k++)
	{
		double Omx = 0.0;
		double Omy = 0.0;

		if (TideForcing)
		{
			double H_Sun = ST - Sa + Lon[k] * rad;
			double H_Moon = ST - Ma + Lon[k] * rad;

			Omx = -SD * (SunAxis / SDist) * (SunAxis / SDist) * (SunAxis / SDist) * (sin(2 * Lat[k] * rad) * sin(2 * Sbeta) * sin(H_Sun) + 2 * cos(Lat[k] * rad)*cos(Lat[k] * rad) * cos(Sbeta)*cos(Sbeta) * sin(2 * H_Sun));
			Omx += -MD * (MoonAxis / MDist)* (MoonAxis / MDist) * (MoonAxis / MDist) * (sin(2 * Lat[k] * rad) * sin(2 * Mbeta) * sin(H_Moon) + 2 * cos(Lat[k] * rad)*cos(Lat[k] * rad) * cos(Mbeta)*cos(Mbeta) * sin(2 * H_Moon));
			Omx *= hlon * rad / hx;

			Omy = SD * (SunAxis / SDist) * (SunAxis / SDist) * (SunAxis / SDist) * (sin(2 * Lat[k] * rad) * (3 * sin(Sbeta)*sin(Sbeta) - 1) + 2 * cos(2 * Lat[k] * rad) * sin(2 * Sbeta) * cos(H_Sun) - sin(2 * Lat[k] * rad) * cos(Sbeta)*cos(Sbeta) * cos(2 * H_Sun));
			Omy += MD * (MoonAxis / MDist)* (MoonAxis / MDist) * (MoonAxis / MDist) * (sin(2 * Lat[k] * rad) * (3 * sin(Mbeta)*sin(Mbeta) - 1) + 2 * cos(2 * Lat[k] * rad) * sin(2 * Mbeta) * cos(H_Moon) - sin(2 * Lat[k] * rad) * cos(Mbeta)*cos(Mbeta) * cos(2 * H_Moon));
			Omy *= hlat * rad / hy;

			/*Omx = -2 * SD * (SunAxis/SDist) * (SunAxis / SDist) * (SunAxis / SDist) * (sin(2 * Lat[k] * rad)*sin(2 * Sbeta)*(sin(ST - Sa + 2 * Lon[k] * rad)) + 2 * cos(Lat[k] * rad)*cos(Lat[k] * rad)*cos(Sbeta)*cos(Sbeta)*(sin(2 * (ST - Sa + 2 * Lon[k] * rad))));
			Omx += -2 * MD * (MoonAxis / MDist)* (MoonAxis / MDist) * (MoonAxis / MDist) * (sin(2 * Lat[k] * rad)*sin(2 * Mbeta)*(sin(ST - Ma + 2 * Lon[k] * rad)) + 2 * cos(Lat[k] * rad)*cos(Lat[k] * rad)*cos(Mbeta)*cos(Mbeta)*(sin(2 * (ST - Ma + 2 * Lon[k] * rad))));
			Omx *= hlon * rad / hx;

			Omy = SD * (SunAxis / SDist) * (SunAxis / SDist) * (SunAxis / SDist) * (sin(2 * Lat[k] * rad) * (3 * sin(Sbeta)*sin(Sbeta) - 1) + 2 * cos(2 * Lat[k] * rad)*sin(2 * Sbeta)*cos(ST - Sa + 2 * Lon[k] * rad) - sin(2 * Lat[k] * rad)*cos(Sbeta)*cos(Sbeta)*cos(2 * (ST - Sa + 2 * Lon[k] * rad)));
			Omy += MD * (MoonAxis / MDist)* (MoonAxis / MDist) * (MoonAxis / MDist) * (sin(2 * Lat[k] * rad) * (3 * sin(Mbeta)*sin(Mbeta) - 1) + 2 * cos(2 * Lat[k] * rad)*sin(2 * Mbeta)*cos(ST - Ma + 2 * Lon[k] * rad) - sin(2 * Lat[k] * rad)*cos(Mbeta)*cos(Mbeta)*cos(2 * (ST - Ma + 2 * Lon[k] * rad)));
			Omy *= hlat * rad / hy;*/
			// 0.000145842 = 2 * 7.2921 / 100000
		}
		ForceX[k] = TideForcing * Omx + fc * 0.000145842 * sin(rad*Lat[k]) * yU[k];
		ForceY[k] = TideForcing * Omy + -fc * 0.000145842 * sin(rad*Lat[k]) * xU[k];
		PhiX[k] = -mu * sqrt(xU[k] * xU[k] + yU[k] * yU[k]) * xU[k];
		PhiY[k] = -mu * sqrt(xU[k] * xU[k] + yU[k] * yU[k]) * yU[k];

		if (windForcing)
		{
			double gamma = 0.001268 *(1.1 + 0.04 * sqrt(xWind[k] * xWind[k] + yWind[k] * yWind[k])) * 0.001;
			PhiX[k] += gamma * xWind[k];
			PhiY[k] += gamma * yWind[k];
		}
	}
};