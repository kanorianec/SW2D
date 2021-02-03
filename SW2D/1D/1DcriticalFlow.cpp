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
	//initConfiguration("testConfig.txt");
	//system("pause");

	/* === PROBLEM STATEMENT === */
 
	// T_begin, T_end - start and end time respectively, in seconds х 
	double T_begin = 0;
	double T_end = 100;

	// start date and time of current problem (UTC)

	int year = 2013;
	int month = 7; // [1,12]
	int day = 1; // [1,31]
	int hour = 12; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] 

	// boundary values of rectangle area
	double x0 = 0.5;
	double xN = 27.5;
	double y0 = 0;
	double yN = 2;

	double lon0 = 0.0;
	double lonN = 0.0;
	double lat0 = 0.0;// 77.866667;
	double latN = 0.0;// 77.866667;

	// physical parameters 
	double mu = 0; // 0.0026; // bottom friction coefficient 
	int fc = 0;    // use Coriolis force (1) or not (0)

	double D = 0.0; // diffusion

	/* === COEFFICIENTS FOR NUMERICAL SOLUTION  === */

	double beta = 0.1; // CFL number (0; 1)
	double alpha = 0.5; // tuning parameter
	double eps = 0.01; // dry zone parameter

	double NS = 0.0; // coefficient for the Navier-Stokes tensor

	int Nx = 28; // 
	int Ny = 3; // 

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

	/* === TECHNICAL PARAMETERS === */

	// folder name in \Data looks like: Test_namePostscript
	string Test_name; // Test name
	string Postscript; // short test description

	// параметры отрисовки
	int Visualization_to_techplot_flag = 0; // вывод результатов дл€ визуализации в Tecplot
	double t_step = 2.0; // (T_end - T_begin) / 20; // интервал вывода данных в файл 

	double sea_level = 10000;  // высота берега, котора€ входит в расчЄт (бќльшие высоты программа игнорирует и не считает)
	
	/* ### SYMMETRY TEST ### */

	Test_name = "1D_criticFlow"; // Test name
	Postscript = "_" + to_string(Nx) + "x" + to_string(Ny); // short test description

	// damBreak
	/*
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			H[k] = 0.1;
			if (x <= 1000.0)
			{
				H[k] = 10.0;
			}				
		}
	}
	*/
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + (i - 1)*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			if (x > 8.0 && x < 12.0)
			{
				B[k] = 0.2 - 0.05*(x - 10.0)*(x - 10.0);
			}
			H[k] = 0.33 - B[k];
		}
		cout << i << " " << x << " " << B[i*Ny + 1] << " " << H[i*Ny + 1] << endl;
	}

	
	/* === INITIALIZATION === */
	/*
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
			
			//if (H[k] < 2.0 && H[k] > eps)
			//xU[k] = 1.0;
		}
	}

	fclose(F);*/
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
	//R1->Initialize_Transport_Problem(C, D);
	R1->SetFixedBoundaryConditions(FLOW, LEFT, 0.18);
	R1->SetFixedBoundaryConditions(HEIGHT, RIGHT, 0.33);
	R1->SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
	//R1->Read_Data_from_file("21600-64800");
	

	R1->Exec_Raschet(); // выполнение расчЄта

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