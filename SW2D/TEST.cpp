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
	initConfiguration("testConfig.txt");
	//system("pause");

	/* === PROBLEM STATEMENT === */
 
	// T_begin, T_end - start and end time respectively, in seconds х 
	double T_begin = 0;
	double T_end = 10;// 0.03;

	// start date and time of current problem (UTC)

	int year = 2013;
	int month = 7; // [1,12]
	int day = 1; // [1,31]
	int hour = 12; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] 

	// boundary values of rectangle area
	double x0 = 0;
	double xN = 10;
	double y0 = 0;
	double yN = 10;

	double lon0 = 0.0;
	double lonN = 0.0;
	double lat0 = 0.0;// 77.866667;
	double latN = 0.0;// 77.866667;

	// physical parameters 
	double mu = 0; // 0.0026; // bottom friction coefficient 
	int fc = 0;    // use Coriolis force (1) or not (0)

	double D = 0.0; // diffusion

	/* === COEFFICIENTS FOR NUMERICAL SOLUTION  === */

	double beta = 0.2; // CFL number (0; 1)
	double alpha = 0.3; // tuning parameter
	double eps = 0.01; // dry zone parameter

	double NS = 1.0; // coefficient for the Navier-Stokes tensor

	int Nx = 10; // 
	int Ny = 10; // 

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
	double t_step = (T_end - T_begin) / 20; // интервал вывода данных в файл 

	double sea_level = 10000;  // высота берега, котора€ входит в расчЄт (бќльшие высоты программа игнорирует и не считает)
	
	/* ### SYMMETRY TEST ### */

	Test_name = "TEST_SYMMETRY"; // Test name
	Postscript = "_" + to_str((T_end - T_begin) / 3600) + "h_NSC_01_" + to_string(Nx) + "x" + to_string(Ny); // short test description

	for (int i = 0; i < Nx; i++)
	{
		//double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			H[k] = 1.0;
			C[k] = 0.1;
			//double y = y0 + j*(yN - y0) / (Ny - 1);
			if ((i == 4 || i == 5) && (j == 4 || j == 5))
			{
				H[k] = 2.0;
				C[k] = 1.0;
			}				
		}
	}
	
	B[2 * Ny + 3] = 1;
	B[2 * Ny + 6] = 1;
	B[7 * Ny + 3] = 1;
	B[7 * Ny + 6] = 1;
	printArray(std::cout, B, Nx, Ny, "B");
	printArray(std::cout, B, Nx, Ny, "H");
	pause();
	//checkSymmetry(H, Nx, Ny, "H0");
	
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
	parallelOpenMP = true;
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

	checkSymmetry(R1->H, Nx, Ny, "H");
	checkSymmetry(R1->xU, Nx, Ny, "xU");
	checkSymmetry(R1->yU, Nx, Ny, "yU");
	checkSymmetry(R1->C, Nx, Ny, "C");

	/*
	Raschet* R2(R1);
	R2->Time_elapsed = 0;
	
	parallelOpenMP = false;
	R2->Exec_Raschet();

	if (checkEquality(R1->H, R2->H, Nx, Ny))
		cout << "OpenMP for H is good!" << endl;
	else
		cout << "OpenMP for H has errors!" << endl;

	if (checkEquality(R1->xU, R2->xU, Nx, Ny))
		cout << "OpenMP for xU is good!" << endl;
	else
		cout << "OpenMP for xU has errors!" << endl;

	if (checkEquality(R1->yU, R2->yU, Nx, Ny))
		cout << "OpenMP for yU is good!" << endl;
	else
		cout << "OpenMP for yU has errors!" << endl;

	if (checkEquality(R1->C, R2->C, Nx, Ny))
		cout << "OpenMP for C is good!" << endl;
	else
		cout << "OpenMP for C has errors!" << endl;

	checkSymmetry(R2->H, Nx, Ny, "H");
	checkSymmetry(R2->xU, Nx, Ny, "xU");
	checkSymmetry(R2->yU, Nx, Ny, "yU");
	checkSymmetry(R2->C, Nx, Ny, "C");
	*/
	
	Raschet* R3(R1);
	R3->Time_elapsed = 0;

	R3->SetOpenBoundaryConditions(LEFT, RIGHT, TOP, BOTTOM);
	for (int i = 0; i < Nx; i++)
	{
		R3->H[i*Ny + 1] = 3.0;
		R3->H[i*Ny + Ny - 2] = 3.0;
	}
	for (int j = 0; j < Ny; j++)
	{
		R3->H[1*Ny + j] = 3.0;
		R3->H[(Nx - 2)*Ny + j] = 3.0;
	}
	
	R3->B[0 * Ny + (int)Ny / 2] = 4.0;
	R3->B[(Nx - 1) * Ny + (int)Ny / 2] = 4.0;
	R3->B[((int)Nx / 2) * Ny + 0] = 4.0;
	R3->B[((int)Nx / 2) * Ny + Ny - 1] = 4.0;
	//R3->B[(5) * Ny + 4] = 4.0;
	
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			if (R3->H[i*Ny + j] < R3->B[i*Ny + j])
				R3->H[i*Ny + j] = 0.0;
		}
	printArray(std::cout, R3->H, Nx, Ny, "H");
	printArray(std::cout, R3->B, Nx, Ny, "B");
	//R3->Numerical_scheme_time_step_parallel();
	R3->Exec_Raschet(); // выполнение расчЄта
	int leftOpen = 0;
	int rightOpen = 0; 
	int topOpen = 0;
	int botOpen = 0;	

	for (int i = 0; i < Nx; i++)
	{
		int k = i*Ny + 0;
		int k1 = i*Ny + 1;
		botOpen += (R3->H[k] != R3->H[k1]) + (R3->xU[k] != R3->xU[k1]) + (R3->yU[k] != R3->yU[k1]);
		k = i*Ny + Ny - 1;
		k1 = i*Ny + Ny - 2;
		topOpen += (R3->H[k] != R3->H[k1]) + (R3->xU[k] != R3->xU[k1]) + (R3->yU[k] != R3->yU[k1]);
	}
	for (int j = 0; j < Ny; j++)
	{
		int k = 0 * Ny + j;
		int k1 = 1 * Ny + j;
		leftOpen += (R3->H[k] != R3->H[k1]) + (R3->xU[k] != R3->xU[k1]) + (R3->yU[k] != R3->yU[k1]);
		k = (Nx - 1)*Ny + j;
		k1 = (Nx - 2)*Ny + j;
		rightOpen += (R3->H[k] != R3->H[k1]) + (R3->xU[k] != R3->xU[k1]) + (R3->yU[k] != R3->yU[k1]);
	}

	cout << "botOpen = " << botOpen << endl;
	cout << "topOpen = " << topOpen << endl;
	cout << "leftOpen = " << leftOpen << endl;
	cout << "rightOpen = " << rightOpen << endl;

	//printArray(R3->H, Nx, Ny, "H");
	//printArray(R3->xU, Nx, Ny, "xU");
	//printArray(R3->yU, Nx, Ny, "yU");
	//printArray(R3->epsilon, Nx, Ny, "epsilon");
	//printArray(R3->xJ, Nx, Ny, "xJ");
	//printArray(R3->yJ, Nx, Ny, "yJ");
	//printFlux(R3->xJ, R3->yJ, Nx, Ny, "J");
	//printFlux(R3->dryFacesX, R3->dryFacesY, Nx, Ny, "dryFaces");
	//printArray(R3->tau, Nx, Ny, "tau");

	system("pause");
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