#include "Raschet.h"
#include "Constants.h"

#include <iostream>
#include "technical.h"

//#include <math.h>
//#include <time.h>
//#include <stdio.h>
//#include <string>
//#include <fstream>
//#include <algorithm>
//#include "Variable.h" // Заголовочный файл описания класса сеточных значений для разностной схемы, НЕ ИСПОЛЬЗУЕТСЯ 

using namespace std;

void Raschet::Prepare_Raschet()
{
	Ephemeris::setLocationOnEarth((float)(lat0 + latN)*0.5, (float)(lon0 + lonN)*0.5); // for ephemeris calculation

	//**************************************************************************
	//	Шаг 0-> Задание основных параметров расчета
	//**************************************************************************
	hx = (xN - x0) / (Nx - 1); // grid step by x
	hy = (yN - y0) / (Ny - 1); // grid step by y

	hlat = (latN - lat0) / (Ny - 1); // grid step by latitude
	hlon = (lonN - lon0) / (Nx - 1); // grid step by longitude

	Stop_Raschet_Flag = 0; // 

	//**************************************************************************
	//	Шаг 1-> Задание матрицы S(i,j) описывающей тип точки, массивов X и Y
	//**************************************************************************

	Raschet::X = new double[Nx*Ny]; // array of coordinate x
	Raschet::Y = new double[Nx*Ny]; // array of coordinate y

	Raschet::Lat = new double[Nx*Ny]; // array of coordinate latitude
	Raschet::Lon = new double[Nx*Ny]; // array of coordinate longitude

	Hmax = H[0]; // зануление массивов для поиска максимумов
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


			/* Маркировка точек отражения */
			// 24 - отражение вверх
			// 15 - отражение вправо
			// 18 - отражение вниз
			// 21 - отражение влево

			if (H[i*Ny + j] > Hmax) { Hmax = H[i*Ny + j]; } // поиск максимумов H и B
			if (B[i*Ny + j] > Bmax) { Bmax = B[i*Ny + j]; } //
			if (B[i*Ny + j] + Bmin > sea_level) { S[i*Ny + j] = EXCLUDED; }
		}
	}

	cout << "Nx = " << Nx << "; Ny = " << Ny << ";" << endl;

	cout << "dx = " << hx << "; dy =  " << hy << ";" << endl;
	cout << "x0 = " << X[0] << "; xN = " << X[Nx*Ny - 1] << ";" << endl;
	cout << "y0 = " << Y[0] << "; yN = " << Y[Nx*Ny - 1] << ";" << endl;
	
	if (hlon + hlat > 0)
		cout << "dlat = " << hlat << "; dlon =  " << hlon << ";" << endl;
	if (lonN + latN > 0)
	{
		cout << "lat0 = " << Lat[0] << "; latN = " << Lat[Ny - 1] << ";" << endl;
		cout << "lon0 = " << Lon[0] << "; lonN = " << Lon[Nx*Ny - 1] << ";" << endl;
	}	

	// угловые точки:
	S[0 * Ny + 0] = LB_CORNER;
	S[(Nx - 1)*Ny + 0] = RB_CORNER;
	S[(Nx - 1)*Ny + Ny - 1] = RT_CORNER;
	S[0 * Ny + Ny - 1] = LT_CORNER;

	//**************************************************************************
	//	Шаг 2->Величина epsilon определяется в каждой точке (i,j)
	//**************************************************************************

	epsilon = new int[Nx*Ny]();
	tau = new double[Nx*Ny];

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			int m = i*Ny + j;
			if (H[m] < eps) {
				epsilon[m] += 1;

				if (H[m] < 0)
					H[m] = pow(10, -9);

				//if (H[m] <= pow(10, -9)) {
				if (i != Nx - 1) {
					epsilon[(i + 1)*Ny + j] += (int)(H[(i + 1)*Ny + j] + B[(i + 1)*Ny + j] < B[m]);
				}
				if (i != 0) {
					epsilon[(i - 1)*Ny + j] += (int)(H[(i - 1)*Ny + j] + B[(i - 1)*Ny + j] < B[m]);
				}
				if (j != Ny - 1) {
					epsilon[i*Ny + j + 1] += (int)(H[i*Ny + j + 1] + B[i*Ny + j + 1] < B[m]);
				}

				if (j != 0) {
					epsilon[i*Ny + j - 1] += (int)(H[i*Ny + j - 1] + B[i*Ny + j - 1] < B[m]);
				}

				if (i != Nx - 1 && j != Ny - 1) {
					epsilon[(i + 1)*Ny + j + 1] += (int)(H[(i + 1)*Ny + j + 1] + B[(i + 1)*Ny + j + 1] < B[m]);
				}

				if (i != Nx - 1 && j != 0) {
					epsilon[(i + 1)*Ny + j - 1] += (int)(H[(i + 1)*Ny + j - 1] + B[(i + 1)*Ny + j - 1] < B[m]);
				}

				if (i != 0 && j != Ny - 1) {
					epsilon[(i - 1)*Ny + j + 1] += (int)(H[(i - 1)*Ny + j + 1] + B[(i - 1)*Ny + j + 1] < B[m]);
				}

				if (i != 0 && j != 0) {
					epsilon[(i - 1)*Ny + j - 1] += (int)(H[(i - 1)*Ny + j - 1] + B[(i - 1)*Ny + j - 1] < B[m]);
				}


				//}		

			}
		}
	}

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			if (!epsilon[i*Ny + j])
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
	/*
	for (int i = 0; i < 8; i++)
	{
		cout << HEIGHT << " " << i << " " << border[HEIGHT][i] << " " << border_C[HEIGHT][i] << endl;
		cout << VELOCITY_X << " " << i << " " << border[VELOCITY_X][i] << " " << border_C[VELOCITY_X][i] << endl;
		cout << VELOCITY_Y << " " <<  << " " << border[VELOCITY_Y][i] << " " << border_C[VELOCITY_Y][i] << endl << endl;
	}*/

	cout << "======================= BOUNDARY CONDITIONS =========================" << endl;
	cout << "    LEFT-TOP:           TOP:                RIGHT-TOP:" << endl;
	cout << GetConditionName(HEIGHT, border[HEIGHT][LT_CORNER], border_C[HEIGHT][LT_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][TOP] , border_C[HEIGHT][TOP]) << GetConditionName(HEIGHT, border[HEIGHT][RT_CORNER], border_C[HEIGHT][RT_CORNER]) << endl;
	cout << GetConditionName(VELOCITY_X, border[VELOCITY_X][LT_CORNER], border_C[VELOCITY_X][LT_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][TOP], border_C[VELOCITY_X][TOP]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RT_CORNER], border_C[VELOCITY_X][RT_CORNER]) << endl;
	cout << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LT_CORNER], border_C[VELOCITY_Y][LT_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][TOP], border_C[VELOCITY_Y][TOP]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RT_CORNER], border_C[VELOCITY_Y][RT_CORNER]) << endl;
	if (TransportProblemFlag) cout << GetConditionName(CONCENTRATION, border[CONCENTRATION][LT_CORNER], border_C[CONCENTRATION][LT_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][TOP], border_C[CONCENTRATION][TOP]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RT_CORNER], border_C[CONCENTRATION][RT_CORNER]) << endl;
	cout << "                   \\ ___________________ /" << endl;
	cout << "                    |                   |" << endl;
	cout << "                    |                   |" << endl;
	cout << "    LEFT:           |                   |    RIGHT:" << endl;
	cout << GetConditionName(HEIGHT, border[HEIGHT][LEFT], border_C[HEIGHT][LEFT]) << "|                   |" << GetConditionName(HEIGHT, border[HEIGHT][RIGHT], border_C[HEIGHT][RIGHT]) << endl;
	cout << GetConditionName(VELOCITY_X, border[VELOCITY_X][LEFT], border_C[VELOCITY_X][LEFT]) << "|                   |" << GetConditionName(VELOCITY_X, border[VELOCITY_X][RIGHT], border_C[VELOCITY_X][RIGHT]) << endl;
	cout << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LEFT], border_C[VELOCITY_Y][LEFT]) << "|                   |" << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RIGHT], border_C[VELOCITY_Y][RIGHT]) << endl;
	if (TransportProblemFlag)
		cout << GetConditionName(CONCENTRATION, border[CONCENTRATION][LEFT], border_C[CONCENTRATION][LEFT]) << "|                   |" << GetConditionName(CONCENTRATION, border[CONCENTRATION][RIGHT], border_C[CONCENTRATION][RIGHT]) << endl;
	cout << "                    |                   |" << endl;
	cout << "                    |                   |" << endl;
	cout << "                    |                   |" << endl;
	cout << "                    |___________________|" << endl;
	cout << "                   /                     \\" << endl;
	cout << "    LEFT-BOTTOM:        BOTTOM:             RIGHT-BOTTOM:" << endl;
	cout << GetConditionName(HEIGHT, border[HEIGHT][LB_CORNER], border_C[HEIGHT][LB_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][BOTTOM], border_C[HEIGHT][BOTTOM]) << GetConditionName(HEIGHT, border[HEIGHT][RB_CORNER], border_C[HEIGHT][RB_CORNER]) << endl;
	cout << GetConditionName(VELOCITY_X, border[VELOCITY_X][LB_CORNER], border_C[VELOCITY_X][LB_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][BOTTOM], border_C[VELOCITY_X][BOTTOM]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RB_CORNER], border_C[VELOCITY_X][RB_CORNER]) << endl;
	cout << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LB_CORNER], border_C[VELOCITY_Y][LB_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][BOTTOM], border_C[VELOCITY_Y][BOTTOM]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RB_CORNER], border_C[VELOCITY_Y][RB_CORNER]) << endl;
	if (TransportProblemFlag) cout << GetConditionName(CONCENTRATION, border[CONCENTRATION][LB_CORNER], border_C[CONCENTRATION][LB_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][BOTTOM], border_C[CONCENTRATION][BOTTOM]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RB_CORNER], border_C[CONCENTRATION][RB_CORNER]) << endl;
	cout << "=====================================================================" << endl;
}
