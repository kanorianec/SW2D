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

//#include <mmsystem.h>

#include <set>

#include <ctime>    
using namespace std;

double Bmin;

int main() {	
	// T_begin, T_end - start and end time respectively, in seconds 
	double T_begin = 0;
	double T_end = T_begin + 1;// 50;
	int num_of_output_data =  10;
							   // параметры отрисовки
	bool Visualization_to_techplot_flag = true; // вывод результатов дл¤ визуализации в Tecplot
	double t_step = (T_end - T_begin) / num_of_output_data; // интервал вывода данных в файл 

	 // start date and time of current problem (UTC)

	/*int year = 2013;
	int month = 7; // [1,12]
	int day = 1; // [1,31]
	int hour = 12; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] */

	// boundary values of rectangle area
	double x0 = 0;
	double xN = 20;
	double y0 = 0;
	double yN = 2;

	/*double lon0 = 51.0208;
	double lonN = 62;
	double lat0 = 68;
	double latN = 71.9917;*/

	// physical parameters
	double mu = 0.0; // 0.0026; // bottom friction coffitient
	int fc = 0; // use Coriolis force (1) or not (0) 

	double beta = 0.2; // число  уранта - подбор(0; 1)
	double alpha = 0.9; // числовой коэффициент, выбираемый из условий точности и устойчивости счета(экспериментально)
	double eps = 0.01; // величина сухого дна, эпсилон 1 сантиметр	
	double D = 0;// 0.01;

	double NS = 1.0; // коээфициент при тензоре Навье-—токса
	
	// параметры численного решения
	int Nx = 401; // разбиение по координате x
	int Ny = 41; // разбиение по координате y

	// название папки в \Data формируется как: Test_namePostscript
	string Test_name = "RNF_perenos"; // название теста 
	string Postscript = "_PROD_D" + to_str(D) + "_alpha_c" + to_str(alpha_c) + "_" + to_str(Nx) + "x" + to_str(Ny);

	double sea_level = 10000;  // высота берега, которая входит в расчёт 

	// выделение памяти для массивов данных
	double* B = new double[Nx*Ny](); // bathymetry
	double* H = new double[Nx*Ny](); // water height
	double* C = new double[Nx*Ny](); // pollution concentration
	double* xU = new double[Nx*Ny](); // velocity along x axis
	double* yU = new double[Nx*Ny](); // velocity along y axis

	double* ForceX = new double[Nx*Ny](); // volume external forces x
	double* ForceY = new double[Nx*Ny](); // volume external forces y
	double* PhiX = new double[Nx*Ny]();// surface external forces x
	double* PhiY = new double[Nx*Ny]();// surface external forces y

	
	int sea_points = 0;  //! что это? 
	//int level_points = 0;  //! что это? 

	int ust_i_end = (Nx - 1) / 20;
	
	
	
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			double y = y0 + j*(yN - y0) / (Ny - 1);

			int k = i*Ny + j;
			H[k] = 1.0;

			//if (j < (Ny-1)/2 && i < ust_i_end)
			//	H[k] = -1.0;
			//else
			//	H[k] = 1;
			C[k] = 0.0;
			xU[k] = 0.0;
			yU[k] = 0.0;
			ForceX[k] = 0.0;
			ForceY[k] = 0.0;
			PhiX[k] = 0.0;
			PhiY[k] = 0.0;
		}
	}
	
	//omp_set_num_threads(1);

	Raschet *R = new Raschet(Test_name,
		Postscript,
		x0,
		xN,
		y0,
		yN,
		0.0, // lat0
		0.0, // latN
		0.0, // lon0
		0.0, // lonN
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

	//R->Restart_from_time_moment(1.000443);

	R->Initialize_Transport_Problem(C, D);
	R->SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);

	R->SetFixedBoundaryConditions(VELOCITY_X, LEFT, 1.0);
	R->SetZeroDerivativeConditions(VELOCITY_X, LT_CORNER, LB_CORNER);
	R->SetWallBoundaryConditions(TOP, BOTTOM);
	R->SetInternalWall(LEFT, ust_i_end, 0, (Ny - 1) / 2);
	R->SetInternalWall(RIGHT, ust_i_end-1, 0, (Ny - 1) / 2);
	R->SetInternalWall(BOTTOM, (Ny - 1) / 2, 0, ust_i_end);
	R->SetInternalWall(TOP, (Ny - 1) / 2 - 1, 0, ust_i_end - 1);
	
	R->SetInternalWall(LEFT, 0, 0, (Ny - 1) / 2-1);

	R->SetFixedBoundaryConditions(CONCENTRATION, LEFT, 1.0);
	
	/*R->SetInternalWall(LEFT, Nx / 2, 0, 2 * Ny / 5 - 1);
	R->SetInternalWall(RIGHT, Nx / 2 - 1, 0, 2 * Ny/ 5 - 1);

	R->SetInternalWall(LEFT, Nx / 2, 3 * Ny / 5, Ny - 1);
	R->SetInternalWall(RIGHT, Nx / 2 - 1, 3* Ny / 5, Ny - 1);*/
	//R1->Read_Data_from_file("21600-64800");
	
	R->Exec_Raschet(); // выполнение расчёта
	//R1->Visualization_to_techplot();

	// очистка выделенной памяти
	delete B;
	delete H;
	delete C;
	delete xU;
	delete yU;
	delete ForceX;
	delete ForceY;
	delete PhiX;
	delete PhiY;

	system("pause");
	return 0;
}