
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
	// T_begin, T_end - start and end time respectively, in seconds 
	double T_begin = 0;
	double T_end = 2 * 24 * 3600;
	int num_of_output_data = 1;// 50;

	int Visualization_to_techplot_flag = 0; //
	double t_step = 3600.0;// (T_end - T_begin) / num_of_output_data; 

	// start date and time of current problem (UTC)

	int year = 2013;
	int month = 7; // [1,12]
	int day = 1; // [1,31]
	int hour = 12; // [0,23]
	int minute = 0; // [0, 59]
	int second = 0; // [0, 59] 

	// boundary values of rectangle area
	double x0 = 0;
	double xN = 0;
	double y0 = 0;
	double yN = 0;

	double lon0 = 51.4250;
	double lonN = 69.2750;
	double lat0 = 68.0083;
	double latN = 71.8;

	double mu = 0.0026; // bottom friction coffitient
	int fc = 1; // use Coriolis force (1) or not (0) 


	double beta = 0.1; // CFL number (0; 1)
	double alpha = 0.5; //
	double eps = 0.1; //

	double NS = 0.0;//1.0; //
	
	int Nx = 715;//1320; // 3240; // 
	int Ny = 456;// 960; // 1080;
	
	int threadsNumber;
	#pragma omp parallel
	{
		threadsNumber = omp_get_num_threads();
	}
	
	//omp_set_num_threads(2);

	string Test_name = "KaraGate_tides_wind_OMP_" + to_str(threadsNumber);//_noForce"; 
	string Postscript = "_" + to_str((T_end - T_begin)/3600) + "h_" + to_string(Nx) + "x" + to_string(Ny); 

	//double t_graph_export = T_begin;  

	double sea_level = 10000;  


	double* B = new double[Nx*Ny](); //
	double* Bn = new double[Nx*Ny](); // 
	double* H = new double[Nx*Ny](); // 
	double* xU = new double[Nx*Ny](); //
	double* yU = new double[Nx*Ny](); //

	double* ForceX = new double[Nx*Ny](); //
	double* ForceY = new double[Nx*Ny](); //
	double* PhiX = new double[Nx*Ny]();//
	double* PhiY = new double[Nx*Ny]();//

	//int sea_points = 0;  //
	//int level_points = 0;  //

	/* === »Õ»÷»¿À»«¿÷»ﬂ ƒ¿ÕÕ€’ ===  GEBCO_2019_30.0_70.0_45.0_63.0_ESRIASCII.asc*/

	FILE *F = fopen("bathymetry/Kara_Gate_extend_715x456.dat", "r");
	if (F == NULL)
	    cout<<"Can't open bathymetry file!"<<endl;
	FILE *FH = fopen("initial/H_in.dat", "r");
	if (FH == NULL) 
	    cout<<"Can't open FH file"<<endl;
	FILE *FU = fopen("initial/U_in.dat", "r");
	FILE *FV = fopen("initial/V_in.dat", "r");
	Bmin = 0;
	// ÎË·Ó ‚ÌÛÚË ˆËÍÎ‡	
	for (int j = Ny - 1; j >= 0; j--)
	{
		for (int i = 0; i < Nx; i++)
		{
			int k = i*Ny + j;

			fscanf(F, "%lf ", &B[k]);
			fscanf(FH, "%lf ", &H[k]);
			fscanf(FU, "%lf ", &xU[k]);
			fscanf(FV, "%lf ", &yU[k]);

			if (B[k] < Bmin)
				Bmin = B[k];
		}
	}

	fclose(F);
	fclose(FH);
	fclose(FU);
	fclose(FV);

	double height_0 = 0.0;

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			int k = i*Ny + j;
			H[k] += -B[k];
			
			if (H[k] < 0)
			{
				H[k] = 0.0;// pow(10, -6);
				xU[k] = 0.0;
				yU[k] = 0.0;
			}
				
			B[k] -= Bmin;
			
		}
	}

	Raschet *R = new Raschet(Test_name,
		Postscript,
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
	//R->SetVisualizationProperties(0, T_end, Nx - 580, Ny - 150, Nx - 1, Ny - 1);
	R->SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
//	R->SetWallBoundaryConditions(TOP,RIGHT);
	R->SetFileBoundaryConditions(VELOCITY_X, /*RIGHT,*/ LEFT, TOP/*, BOTTOM*/);
	R->SetFileBoundaryConditions(VELOCITY_Y, /*RIGHT,*/ LEFT, TOP/*, BOTTOM*/);

	R->SetWindSpeed(0.0, 3600);// 3600);
	//R->SetFileBoundaryConditions(HEIGHT, /*RIGHT,*/ LEFT, TOP/*, BOTTOM*/);
	/*
	R->SetFileBoundaryConditions(VELOCITY_X, LEFT);//, TOP, BOTTOM);
	R->SetFileBoundaryConditions(VELOCITY_Y, LEFT);
	R->SetFileBoundaryConditions(HEIGHT, LEFT);
	*/
	
	//R->SetWallBoundaryConditions(RIGHT, LEFT, TOP, BOTTOM);
	//R->SetFixedBoundaryConditions(VELOCITY_X, RIGHT, -0.01);

	R->Exec_Raschet();

	delete[] B;
	delete[] H;
	delete[] xU;
	delete[] yU;
	delete[] ForceX;
	delete[] ForceY;
	delete[] PhiX;
	delete[] PhiY;
	
	return 0;
}