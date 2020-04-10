/*
Description file of main class "Raschet" and its functions
*/

#define _USE_MATH_DEFINES
#define _XOPEN_SOURCE 600

#include "Raschet.h"
#include "technical.h"
#include "Constants.h"

#include <iostream>
#include <fstream>
#include <omp.h>
#include <stdarg.h>

#include <cmath>

//#include <time.h>
//#include <stdio.h>
//#include <string>
#include <algorithm>

using namespace std;

// the constructor of the class "Raschet" is inherited from the "Problem_Definition" class
Raschet::Raschet(string Test_name,
	string Postscript,
	double x0,
	double xN,
	double y0,
	double yN,
	double lat0,
	double latN,
	double lon0,
	double lonN,
	int Nx,
	int Ny,
	double T_begin,
	double T_end,
	double alpha,
	double beta,
	double mu,
	int fc,
	double NS,
	double Bmin,
	double eps,
	double B[],
	double H[],
	double xU[],
	double yU[],
	double ForceX[],
	double ForceY[],
	double PhiX[],
	double PhiY[],
	double t_step,
	bool Visualization_to_techplot_flag,
	double sea_level
) : Problem_Definition(x0,
	xN,
	y0,
	yN,
	lat0,
	latN,
	lon0,
	lonN,
	Nx,
	Ny,
	T_begin,
	T_end,
	B,
	H,
	xU,
	yU,
	ForceX,
	ForceY,
	PhiX,
	PhiY
)
{
	Raschet::Test_name = Test_name;
	Raschet::Postscript = Postscript;
	Raschet::alpha = alpha;
	Raschet::beta = beta;
	Raschet::NS = NS;
	Raschet::mu = mu;
	Raschet::fc = fc;
	Raschet::Bmin = Bmin;
	Raschet::eps = eps;
	Raschet::t_step = t_step;
	Raschet::Visualization_to_techplot_flag = Visualization_to_techplot_flag;

	int N = Nx*Ny;
	Raschet::RaschetTime = -1;

	Raschet::Polar = 0;

	Raschet::xUt = new double[N]();
	Raschet::yUt = new double[N]();
	Raschet::Ht = new double[N]();
	//Raschet::PhiXt = new double[N]();
	//Raschet::PhiYt = new double[N]();
	Raschet::S = new int[N]();
	Raschet::sea_level = sea_level;
	//Raschet::Time_elapsed = T_begin;

	Raschet::TransportProblemFlag = false;

	Raschet::path = Data_folder + "/" + Test_name + Postscript;

	Raschet::BoundaryConditionsFromFile = false;
	Raschet::InternalWallsFlag = false;

	Raschet::restart = false;
	Raschet::windForcing = false;
	
	Raschet::SetOpenBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT, RT_CORNER, LT_CORNER, RB_CORNER, LB_CORNER);
	Raschet::SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
}

Raschet::Raschet(string Test_name,
	string Postscript,
	double lat0,
	double latN,
	double lon0,
	double lonN,
	int Nx,
	int Ny,
	double T_begin,
	double T_end,
	double alpha,
	double beta,
	double mu,
	int fc,
	double NS,
	double Bmin,
	double eps,
	double B[],
	double H[],
	double xU[],
	double yU[],
	double ForceX[],
	double ForceY[],
	double PhiX[],
	double PhiY[],
	double t_step,
	bool Visualization_to_techplot_flag,
	double sea_level
) : Problem_Definition(0.0,
	polar_to_decart_x(lonN - lon0, 0.5*(lat0 + latN)),
	0.0,
	polar_to_decart_y(latN - lat0),
	lat0,
	latN,
	lon0,
	lonN,
	Nx,
	Ny,
	T_begin,
	T_end,
	B,
	H,
	xU,
	yU,
	ForceX,
	ForceY,
	PhiX,
	PhiY
)
{
	Raschet::Test_name = Test_name;
	Raschet::Postscript = Postscript;

	Raschet::Polar = 1;

	Raschet::alpha = alpha;
	Raschet::beta = beta;
	Raschet::NS = NS;
	Raschet::mu = mu;
	Raschet::fc = fc;
	Raschet::Bmin = Bmin;
	Raschet::eps = eps;
	Raschet::t_step = t_step;
	Raschet::Visualization_to_techplot_flag = Visualization_to_techplot_flag;

	int N = Nx*Ny;

	Raschet::RaschetTime = -1;

	Raschet::xUt = new double[N](); // Initialization by zero
	Raschet::yUt = new double[N]();
	Raschet::Ht = new double[N]();
	//Raschet::PhiXt = new double[N]();
	//Raschet::PhiYt = new double[N]();
	Raschet::S = new int[N]();
	Raschet::sea_level = sea_level;
	//Raschet::Time_elapsed = T_begin;

	Raschet::TransportProblemFlag = false;

	Raschet::path = Data_folder + "/" + Test_name + Postscript;

	Raschet::BoundaryConditionsFromFile = false;
	Raschet::InternalWallsFlag = false;

	Raschet::restart = false;
	Raschet::windForcing = false;

	Raschet::SetOpenBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT, RT_CORNER, LT_CORNER, RB_CORNER, LB_CORNER);
	Raschet::SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
}


void Raschet::Exec_Raschet()
{
	Raschet::Prepare_Folder("Data"); // Checking and creating "Data" folder

	Raschet::Prepare_Folder(path, false || restart); // Creating Test_name folder
	Raschet::Prepare_Raschet(); // подготовка расчёта:
	cout << "Prepare_Raschet" << endl;
	Raschet::Perform_Calculations(); // выполнение расчёта
}

void Raschet::Perform_Calculations()
{		
	HourMark = ((int)T_begin / 3600);
	dT = ((hx + hy)*0.5*beta) / sqrt(gc*(Hmax + 10));
	
	write_extra_inf(cout,0); // output metadata to the screen
	
	write_extra_inf_to_file(0); // output metadata to the special file extra_inf.txt

	double scheme_time = omp_get_wtime(); // time of calculations 
	double visualization_time = omp_get_wtime(); // time of visualisation 

	outputInputs();
	
	visualization_time = omp_get_wtime() - visualization_time;
	int est_time = 1;
	double estimation_time = omp_get_wtime();

	while (T_end>Time_elapsed && !Stop_Raschet_Flag)
	{
		Numerical_scheme_time_step_parallel();

		if (est_time)
		{
			std::cout << "Estimation of time = " << (omp_get_wtime() - estimation_time)*(T_end - T_begin)/dT + 2 * visualization_time * ((T_end - T_begin)/t_step + 1)<< " seconds." << endl;
			est_time = 0;
		}

		Time_elapsed = Time_elapsed + dT;

		if (Time_elapsed > 3600 * (HourMark + 1))
			HourMark++;

		if (Time_elapsed >= t_graph_export || Stop_Raschet_Flag)
		{
			outputResults();
			t_graph_export = t_graph_export + t_step;
		}
	}

	if (Stop_Raschet_Flag)
		outputResults();

	double Time_of_work = (omp_get_wtime() - scheme_time); 
	std::cout << "Time of work = " << Time_of_work << " seconds." << endl;

	std::ofstream fLog(path + "/timeLog.dat", std::ios::out | std::ios::app);
	fLog << "Time of work = " << Time_of_work << " seconds." << endl;
	fLog.close();
	write_extra_inf_to_file(Time_of_work);
}

Raschet::~Raschet()
{
	/*delete[] B;
	delete[] H;
	delete[] xU;
	delete[] yU;
	delete[] ForceX;
	delete[] ForceY;
	delete[] PhiX;
	delete[] PhiY;*/

	delete[] Ht;
	delete[] xUt;
	delete[] yUt;

	for (int PType = BOTTOM; PType <= LEFT; PType++)
	{
		for (int VType = HEIGHT; VType <= VELOCITY_Y; VType++)
		{
			if (border[VType][PType] == FROM_FILE)
			{
				fclose(FV[VType][PType]);
				delete[] lin_b[VType][PType];
				delete[] lin_k[VType][PType];
			}
		}
		if (TransportProblemFlag)
		{
			if (border[CONCENTRATION][PType] == FROM_FILE)
			{
				fclose(FV[CONCENTRATION][PType]);
				delete[] lin_b[CONCENTRATION][PType];
				delete[] lin_k[CONCENTRATION][PType];
			}
		}
	}

	if (windForcing)
	{
		FWindX.close();
		FWindY.close();

		delete[] xWind;
		delete[] yWind;
	}
}

void Raschet::SetWindSpeed(double WindFrictionCoefficient, double period)
{
	windFrictionCoef = WindFrictionCoefficient;
	timeWind = T_begin - 1;
	timeWindPeriod = period;
	windForcing = true;

	FWindX.open("windSpeed/xWind.dat", std::ios::binary);
	FWindY.open("windSpeed/yWind.dat", std::ios::binary);

	if (!FWindX.is_open() || !FWindY.is_open())
		std::cout << "Error: windSpeed/x(y)Wind was not found!" << endl;
	else
	{
		xWind = new double[Nx*Ny]();
		yWind = new double[Nx*Ny]();
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
		//system("pause");

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
			//system("pause");
	}

	if (windForcing)
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
		PhiX[k] = - mu * sqrt(xU[k] * xU[k] + yU[k] * yU[k]) * xU[k];
		PhiY[k] = - mu * sqrt(xU[k] * xU[k] + yU[k] * yU[k]) * yU[k];

		if (windForcing)
		{
			double gamma = 0.001268 *(1.1 + 0.04 * sqrt(xWind[k] * xWind[k] + yWind[k] * yWind[k])) * 0.001;
			PhiX[k] += gamma * xWind[k];
			PhiY[k] += gamma * yWind[k];
		}
	}
};




