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

#include "../MainCore/Problem_Definition.h"
#include "../MainCore/Raschet.h"
#include "../MainCore/technical.h"

#include <ctime>   
using namespace std;

double Bmin;


int main() {
	/* === Problem Definition === */

	// T_begin, T_end - start and end time respectively, in seconds 
	double T_begin = 0;
	double T_end = T_begin + 20;
	int num_of_output_data = 50;// 50;

	int Visualization_to_techplot_flag = 1; //
	double t_step = (T_end - T_begin) / num_of_output_data; 

	// start date and time of current problem (UTC)

	int year = 2010;
	int month = 1; // [1,12]
	int day = 30; // [1,31]
	int hour = 15; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] 

	// boundary values of rectangle area
	double x0 = 0;
	double xN = 100;
	double y0 = 0;
	double yN = 100;

	double lat0 = 0;// 55.0001;
	double latN = 0;// 55.0002;
	double lon0 = 0;// 36.0001;
	double lonN = 0;// 36.0002;

	// физические параметры 
	double mu = 0.0; // bottom friction coffitient
	int fc = 0; // use Coriolis force (1) or not (0) 

	/* === ѕќƒЅ»–ј≈ћџ≈  ќЁ‘‘»÷»≈Ќ“џ „»—Ћ≈ЌЌќ√ќ –≈Ў≈Ќ»я === */

	double beta = 0.1; // CFL number (0; 1)
	double alpha = 0.5; //
	double eps = 1e-2; //

	double NS = 1.0; // коээфициент при тензоре Ќавье-—токса

	int Nx = 100; // 
	int Ny = 100; //

	if (!parallelOpenMP)
	{
		omp_set_num_threads(1);
	}

	int threadsNumber;
	#pragma omp parallel
	{
		threadsNumber = omp_get_num_threads();
	}

	// название папки в \Data формируетс€ как: Test_namePostscript
	string Test_name = "Dry_Zone_Problem"; // название теста 
	string Postscript = "_v0_" + to_string(Nx) + "x" + to_string(Ny); // „ем уникален тест, дл€ отличи€ от остальных


	double sea_level = 10000;  // высота берега, котора€ входит в расчЄт (бќльшие высоты программа игнорирует и не считает)

	double* B = new double[Nx*Ny](); // дно
	double* H = new double[Nx*Ny](); // высота сло€
	double* xU = new double[Nx*Ny](); // скорость по x
	double* yU = new double[Nx*Ny](); // скорость по y

	double* ForceX = new double[Nx*Ny](); // внешние объЄмные силы по x
	double* ForceY = new double[Nx*Ny](); // внешние объЄмные силы по y
	double* PhiX = new double[Nx*Ny]();// внешние поверхностные силы по x
	double* PhiY = new double[Nx*Ny]();// внешние поверхностные силы по y

	//int sea_points = 0;  //!!! что это? 
	//int level_points = 0;  //!!! что это? 

	/* === »Ќ»÷»јЋ»«ј÷»я ƒјЌЌџ’ === */
	double RR = 10.0;
	double xc = 50.0;
	double yc = 50.0;
	double HGT = 30.0;

	// либо внутри цикла	
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			double y = y0 + j*(yN - y0) / (Ny - 1);
			B[k] = 10.0*(x - 50.0)*(x - 50.0) / (50.0*30.0) + 10.0*(y - 50.0)*(y - 50.0) / (50.0*30.0);
			if (y <= 35 && y >= 30)
				B[k] += 5;
			H[k] = 5.0 - B[k];

			if ((x - xc)*(x - xc) + (y - yc)*(y - yc) < RR * RR)
				H[k] = HGT;
			if (H[k] < 0)
				H[k] = 0.0;
			//H[k] = 5.0;
		}
	}

	
	/* === ¬џѕќЋЌ≈Ќ»≈ –ј—„≈“ј === */
	Raschet *R = new Raschet(Test_name,
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
	
	R->SetStartTime(year, month, day, hour, minute, second);
	R->SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
	
	//R->SetFixedBoundaryConditions(VELOCITY_X, RIGHT, -5);
	//R->SetWallBoundaryConditions(TOP, BOTTOM);
	//R->SetFixedBoundaryConditions(HEIGHT, RIGHT, 0.5);
	//R->SetFixedBoundaryConditions(VELOCITY_X, RIGHT, -5);

	// boundary conditions setting, OPEN BOUNDARY is default
	// LOL
	//R->SetWallBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT);

	//pause();
	/*
	Raschet *R2 = new Raschet(Test_name,
	Postscript,
	x0,
	xN,
	y0,
	yN,
	0.0, //lat0
	0.0, //latN
	0.0, //lon0
	0.0, //lonN
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
	t_graph_export,
	Visualization_to_techplot_flag,
	sea_level
	);
	*/

	//R1->Read_Data_from_file("21600-64800");


	R->Exec_Raschet(); // выполнение расчЄта

					   /* === ќ„»—“ ј ѕјћя“» === */

	delete B;
	delete H;
	delete xU;
	delete yU;
	delete ForceX;
	delete ForceY;
	delete PhiX;
	delete PhiY;

	pause();
	return 0;
}