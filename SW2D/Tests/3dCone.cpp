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
#include <algorithm>
#include <cmath>

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
	double T_end = T_begin + 100;// 00;
	int num_of_output_data = 100;// 50;

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
	double x0 = -37.5;
	double xN = 37.5;
	double y0 = -15;
	double yN = 15;

	double lat0 = 0;// 55.0001;
	double latN = 0;// 55.0002;
	double lon0 = 0;// 36.0001;
	double lonN = 0;// 36.0002;

	// физические параметры 
	double mu = 0.0;// 0.002; // bottom friction coffitient
	int fc = 0; // use Coriolis force (1) or not (0) 

	/* === ѕќƒЅ»–ј≈ћџ≈  ќЁ‘‘»÷»≈Ќ“џ „»—Ћ≈ЌЌќ√ќ –≈Ў≈Ќ»я === */

	double beta = 0.1; // CFL number (0; 1)
	double alpha = 0.5; //
	double eps = 1e-3; //

	double NS = 2.0; // коээфициент при тензоре Ќавье-—токса

	int Nx = 76; // 
	int Ny = 31; //

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
	string Test_name = "3Cone_wellB"; // название теста 
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
	double RR = 2.0;
	double xc = 10.0;
	double yc = 10.0;
	double HGT = 4.0;	
	
	//	#funkySetFields - create - field h - expression 'pos().x <= -21.5 ? 1.875 : 0.0' - dimension '[0 1 0 0 0 0 0]' - time 0
	//	funkySetFields - create - field h - expression 'max(10.0 - b, 0.0)' - dimension '[0 1 0 0 0 0 0]' - time 0
	//
	//	#funkySetFields - create - field b - expression '1' - condition 'pos().x > 10.0 && pos().x < 15 && pos().y > 10 && pos().y <15' - dimension '[0 1 0 0 0 0 0]' - time 0
	//	#funkySetFields - create - field h - expression '(pos().x)*(pos().x) + (pos().y)*(pos().y) <= 4  ? 4 : 0' - dimension '[0 1 0 0 0 0 0]' - time 0(base)ivanov@iva:~/ ISP / RSWE / tutori

	// либо внутри цикла	
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			int k = i*Ny + j;
			double y = y0 + j*(yN - y0) / (Ny - 1);
			B[k] = max(max(max(0.0, 3.0 - 0.3*sqrt(sqr(x - 10.0) + sqr(y))), 1.0 - sqrt(sqr(x + 7.5) + sqr(y - 9.0)) / 8.0), 1.0 - sqrt(sqr(x + 7.5) + sqr(y + 9.0)) / 8.0);

			//B[k] = 0.0;// (x - 10.0)*(x - 10.0) / 40.0 + (y - 10.0)*(y - 10.0) / 40.0;
			//if (y <= 35 && y >= 30)
			//	B[k] += 5;
			//H[k] = 0.5 - B[k];
			
			if (x <= -21.5)
				H[k] = 1.875;
			//if (H[k] < 0)
			//H[k] += 1.0;
		}
	}
	//checkSymmetry(H, Nx, Ny, "H");
	
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
	R->SetWallBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT);

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

	//checkSymmetry(H, Nx, Ny, "H");
	//printArray(H, Nx, Ny, "H");

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