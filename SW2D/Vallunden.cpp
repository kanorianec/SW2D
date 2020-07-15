#define _CRT_SECURE_NO_WARNINGS

#define _USE_MATH_DEFINES
#define _XOPEN_SOURCE 600

#include <cmath>
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "MainCore/Problem_Definition.h"
#include "MainCore/Raschet.h"
#include "MainCore/technical.h"

#include <ctime>   
using namespace std;

double Bmin;

int main() {	

	/* === PROBLEM STATEMENT === */
 
	// T_begin, T_end - start and end time respectively, in seconds х 
	double T_begin = 0;
	double T_end = T_begin + 24 * 3600;

	// start date and time of current problem (UTC)

	int year = 2013;
	int month = 7; // [1,12]
	int day = 1; // [1,31]
	int hour = 12; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] 

	// boundary values of rectangle area
	double x0 = 0;
	double xN = 1000;
	double y0 = 0;
	double yN = 1250;

	double lon0 = 0.0;
	double lonN = 77.866667;
	double lat0 = 77.866667;
	double latN = 0.0;

	// physical parameters 
	double mu = 0; // 0.0026; // bottom friction coefficient 
	int fc = 1;    // use Coriolis force (1) or not (0)

	double D = 0.0; // diffusion

	/* === COEFFICIENTS FOR NUMERICAL SOLUTION  === */

	double beta = 0.2; // CFL number (0; 1)
	double alpha = 0.3; // tuning parameter
	double eps = 0.01; // dry zone parameter

	double NS = 1.0; // coefficient for the Navier-Stokes tensor

	int Nx = 247; // 
	int Ny = 234; // 

	if (!parallelOpenMP)
	{
		omp_set_num_threads(1);
	}

	int threadsNumber;
	#pragma omp parallel
	{
		threadsNumber = omp_get_num_threads();
	}

	/* === TECHNICAL PARAMETERS === */

	// folder name in \Data looks like: Test_namePostscript
	string Test_name = "Vallunden_NEW"; // Test name
	string Postscript = "_" + to_str((T_end - T_begin) / 3600) + "h_NSC_01_" + to_string(Nx) + "x" + to_string(Ny); // short test description

																			  // параметры отрисовки
	int Visualization_to_techplot_flag = 0; // вывод результатов дл€ визуализации в Tecplot
	double t_step = (T_end - T_begin) / 24; // интервал вывода данных в файл 

	double sea_level = 11;  // высота берега, котора€ входит в расчЄт (бќльшие высоты программа игнорирует и не считает)

	// выделение пам€ти дл€ массивов данных
	double* B = new double[Nx*Ny](); // дно
	double* H = new double[Nx*Ny](); // высота сло€
	double* C = new double[Nx*Ny](); // концентраци€ примеси
	double* xU = new double[Nx*Ny](); // скорость по x
	double* yU = new double[Nx*Ny](); // скорость по y

	double* ForceX = new double[Nx*Ny](); // внешние силы по x
	double* ForceY = new double[Nx*Ny](); // внешние силы по y
	double* PhiX = new double[Nx*Ny]();// !!! сила ветра по x
	double* PhiY = new double[Nx*Ny]();// !!! сила ветра по y

	//int sea_points = 0;  //! что это? 
	//int level_points = 0;  //! что это? 
	//double sea_level = 1000;  // высота берега, котора€ входит в расчЄт 
	
	/* === INITIALIZATION === */
	
	FILE *F = fopen("bathymetry/LakeValunden_247x234.dat", "r");
	if (F == NULL)
		cout << "Can't open bathymetry file!" << endl;

	double height = 10.0; // средн€€ глубина 10-12 метров

						  // либо внутри цикла	
	for (int j = Ny - 1; j >= 0; j--)
	{
		double y = y0 + j*(yN - y0) / (Ny - 1);
		for (int i = 0; i < Nx; i++)
		{
			double x = x0 + i*(xN - x0) / (Nx - 1);
			int k = i*Ny + j;

			fscanf(F, "%lf", &B[k]);

			H[k] = pow(10, -9);
			C[k] = 0.0;
			if (B[k] < height)
				H[k] = height - B[k];
			/*
			if (H[k] < 2.0 && H[k] > eps)
			xU[k] = 1.0;*/
		}
	}

	fclose(F);
	
	Raschet *R1 = new Raschet(Test_name,
		Postscript, 
		x0,
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
		alpha,
		beta,
		mu,
		fc,
		NS,		
		Bmin,
		eps,
		B,
		H,
		xU,
		yU,
		ForceX,
		ForceY,
		PhiX,
		PhiY,
		t_step,
		Visualization_to_techplot_flag,
		sea_level
	);
	R1->Initialize_Transport_Problem(C, D);
	R1->SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);

	//R1->Read_Data_from_file("21600-64800");
	

	R1->Exec_Raschet(); // выполнение расчЄта
	//R1->Visualization_to_techplot();

	// очистка выделенной пам€ти
	delete[] B;
	delete[] H;
	delete[] C;
	delete[] xU;
	delete[] yU;
	delete[] ForceX;
	delete[] ForceY;
	delete[] PhiX;
	delete[] PhiY;

	//PlaySound(TEXT("Alarm06.wav"), NULL, SND_FILENAME | SND_LOOP | SND_ASYNC);
	system("pause");
	return 0;
}