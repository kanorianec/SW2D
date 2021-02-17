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

	Raschet::mu = mu;
	Raschet::fc = fc;
	Raschet::windForcing = NO_FORCE;
	Raschet::tidesHarmonics = false;

	Raschet::restart = false;
	
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

	Raschet::mu = mu;
	Raschet::fc = fc;
	Raschet::windForcing = NO_FORCE;
	Raschet::tidesHarmonics = false;

	Raschet::restart = false;

	Raschet::SetOpenBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT, RT_CORNER, LT_CORNER, RB_CORNER, LB_CORNER);
	Raschet::SetVisualizationProperties(T_begin, T_end, 0, 0, Nx - 1, Ny - 1);
}


void Raschet::Exec_Raschet()
{
	if (!parallelOpenMP)
	{
		omp_set_num_threads(1);
	}
	#pragma omp parallel
	{
		Raschet::OMP_THREADS_NUMBER = omp_get_num_threads();
	}
	Raschet::Prepare_Folder("Data"); // Checking and creating "Data" folder

	Raschet::Prepare_Folder(path, ignore_warning || restart); // Creating Test_name folder
	Raschet::Prepare_Raschet(); // подготовка расчёта:
	cout << "Prepare_Raschet" << endl;
	Raschet::Perform_Calculations(); // выполнение расчёта
}

void Raschet::Perform_Calculations()
{		
	HourMark = ((int)T_begin / 3600);
	dT = ((hx + hy)*0.5*beta) / sqrt(gc*(Hmax)); // + 10
	dT = 0.005;
	write_extra_inf(cout,0); // output metadata to the screen
	
	write_extra_inf_to_file(0); // output metadata to the special file extra_inf.txt

	double scheme_time = omp_get_wtime(); // time of calculations 
	double visualization_time = omp_get_wtime(); // time of visualisation 

	outputInputs();
	
	visualization_time = omp_get_wtime() - visualization_time;
	int est_time = 1;
	double estimation_time = omp_get_wtime();

	while (T_end>=Time_elapsed && !Stop_Raschet_Flag)
	{
		Numerical_scheme_time_step_parallel();
		//int Sym = 1;
		//Sym *= (int)checkSymmetry(H, Nx, Ny, "H");
		//Sym *= (int)checkSymmetry(xU, Nx, Ny, "xU");
		//Sym *= (int)checkSymmetry(yU, Nx, Ny, "yU");
		////printFlux(xJ, yJ, Nx, Ny, "J");
		////printFlux(dryFacesX, dryFacesY, Nx, Ny, "dryFaces");
		//if (!Sym)
		//{
		//	printArray(H, Nx, Ny, "H");
		//	printArray(xU, Nx, Ny, "xU");
		//	printArray(yU, Nx, Ny, "yU");
		//	//printFlux(xJ, yJ, Nx, Ny, "J");
		//	//printFlux(dryFacesX, dryFacesY, Nx, Ny, "dryFaces");
		//	printArray(epsilon, Nx, Ny, "epsilon");
		//	pause();
		//}
			
		//printArray(H, Nx, Ny, "H");
		//printArray(xU, Nx, Ny, "xU");
		//printArray(yU, Nx, Ny, "yU");
		//printFlux(xJ, yJ, Nx, Ny, "J");
		//printFlux(dryFacesX, dryFacesY, Nx, Ny, "dryFaces");
		//printArray(epsilon, Nx, Ny, "epsilon");
		//pause();
		Time_elapsed = Time_elapsed + dT;

		if (est_time)
		{
			std::cout << "Estimation of time = " << (omp_get_wtime() - estimation_time)*(T_end - T_begin)/dT + 2 * visualization_time * ((T_end - T_begin)/t_step + 1)<< " seconds." << endl;
			est_time = 0;
		}		

		if (Time_elapsed > 3600 * (HourMark + 1))
			HourMark++;

		if (Time_elapsed >= t_graph_export || Time_elapsed + dT > t_graph_export || Stop_Raschet_Flag)
		{
			outputResults();
			t_graph_export = t_graph_export + t_step;
		}
	}
	/*
	for (int i = 0; i < Nx; i++)
	{
		cout << i << " " << Ht[i*Ny + 1] << endl;
	}*/

	//if (Stop_Raschet_Flag)
	//	outputResults();

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

	if (tidesHarmonics)
	{
		for (int PType = BOTTOM; PType <= LEFT; PType++)
		{
			if (tideSide[PType])
			{
				for (int h = 0; h < tideNum; h++)
				{
					delete[] tideA[PType][h];
				}
			}
		}
	}
}







