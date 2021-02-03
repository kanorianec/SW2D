#include "Raschet.h"
#include "Constants.h"

#include <iostream>
#include "technical.h"

#include <cmath>
//#include <time.h>
//#include <stdio.h>
//#include <string>
//#include <fstream>
//#include <algorithm>

using namespace std;

void Raschet::Prepare_Raschet()
{
	Ephemeris::setLocationOnEarth((float)(lat0 + latN)*0.5, (float)(lon0 + lonN)*0.5); // for ephemeris calculation

	Time_elapsed = T_begin;

	//**************************************************************************
	//	Ўаг 0-> «адание основных параметров расчета
	//**************************************************************************
	hx = (xN - x0) / (Nx - 1); // grid step by x
	hy = (yN - y0) / (Ny - 1); // grid step by y

	hlat = (latN - lat0) / (Ny - 1); // grid step by latitude
	hlon = (lonN - lon0) / (Nx - 1); // grid step by longitude

	Stop_Raschet_Flag = 0; // 

	//**************************************************************************
	//	Ўаг 1-> «адание матрицы S(i,j) описывающей тип точки, массивов X и Y
	//**************************************************************************

	Raschet::X = new double[Nx*Ny]; // array of coordinate x
	Raschet::Y = new double[Nx*Ny]; // array of coordinate y

	Raschet::Lat = new double[Nx*Ny]; // array of coordinate latitude
	Raschet::Lon = new double[Nx*Ny]; // array of coordinate longitude

	Hmax = H[0]; // зануление массивов дл€ поиска максимумов
	Bmax = B[0]; //

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			X[k] = hx*i + x0;
			Y[k] = hy*j + y0;
			Lat[k] = hlat*j + lat0;
			Lon[k] = hlon*i + lon0;
				

			if (S[k] <INTERNALWALL)
			{
				S[k] = INTERNAL;
				if (j == 0)
				{
					S[k] = BOTTOM;
				}
				if (i == Nx - 1)
				{
					S[k] = RIGHT;
				}
				if (j == Ny - 1)
				{
					S[k] = TOP;
				}
				if (i == 0)
				{
					S[k] = LEFT;
				}

				if (H[k] == -1.0)
				{
					S[k] = EXCLUDED;
				}
			}	


			/* ћаркировка точек отражени€ */
			// 24 - отражение вверх
			// 15 - отражение вправо
			// 18 - отражение вниз
			// 21 - отражение влево

			if (H[i*Ny + j] > Hmax) { Hmax = H[i*Ny + j]; } // поиск максимумов H и B
			if (B[i*Ny + j] > Bmax) { Bmax = B[i*Ny + j]; } //
			//if (B[i*Ny + j] + Bmin > sea_level) { S[i*Ny + j] = EXCLUDED; }
		}
	}

	// угловые точки:
	S[0 * Ny + 0] = LB_CORNER;
	S[(Nx - 1)*Ny + 0] = RB_CORNER;
	S[(Nx - 1)*Ny + Ny - 1] = RT_CORNER;
	S[0 * Ny + Ny - 1] = LT_CORNER;

	//**************************************************************************
	//	Ўаг 2->¬еличина epsilon определ€етс€ в каждой точке (i,j)
	//**************************************************************************

	epsilon = new int[Nx*Ny]();
	tau = new double[Nx*Ny];

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			int m = i*Ny + j;
			if (H[m] <= eps)
				epsilon[m] += 1;

			if (i != Nx - 1) {
				if (H[m] < eps)
					dryFacesX[m] += (int)(H[(i + 1)*Ny + j] + B[(i + 1)*Ny + j] < B[m] + eps || H[(i + 1)*Ny + j] < eps);
				else if (H[(i + 1)*Ny + j] < eps)
					dryFacesX[m] += (int)(eps + B[(i + 1)*Ny + j] > B[m] + H[m]);

				epsilon[m] += (int)(H[(i + 1)*Ny + j] < eps);
			}
			if (i != 0) {
				epsilon[m] += (int)(H[(i - 1)*Ny + j] < eps);
			}
			if (j != Ny - 1) {
				if (H[m] < eps)
					dryFacesY[m] += (int)(H[i*Ny + j + 1] + B[i*Ny + j + 1] < B[m] + eps || H[i*Ny + j + 1] < eps);
				else if (H[i*Ny + j + 1] < eps)
					dryFacesY[m] += (int)(eps + B[i*Ny + j + 1] > B[m] + H[m]);

				epsilon[m] += (int)(H[i*Ny + j + 1] < eps);
			}

			if (j != 0) {
				epsilon[m] += (int)(H[i*Ny + j - 1] < eps);
			}
			/*
			if (i != Nx - 1) {
				if (H[m] < eps)
					dryFacesX[m] += (int)(H[(i + 1)*Ny + j] + B[(i + 1)*Ny + j] < B[m] + eps);
				if (H[(i + 1)*Ny + j] < eps)
					dryFacesX[m] += (int)(eps + B[(i + 1)*Ny + j] > B[m] + H[m]);
			}
			if (j != Ny - 1) {
				if (H[m] < eps)
					dryFacesY[m] += (int)(H[i*Ny + j + 1] + B[i*Ny + j + 1] < B[m] + eps);
				if (H[i*Ny + j + 1] < eps)
					dryFacesY[m] += (int)(eps + B[i*Ny + j + 1] > B[m] + H[m]);
			}
			*/
		}
	}

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			if (H[i*Ny + j] > eps)
			{
				tau[i*Ny + j] = alpha*sqrt(hx*hy) / sqrt(gc*H[i*Ny + j]);
			}
			else
			{
				//tau[i*Ny + j] = alpha*sqrt(hx*hy) / sqrt(gc*eps/2);
				tau[i*Ny + j] = 0.0;
			}
		}
	}
}
