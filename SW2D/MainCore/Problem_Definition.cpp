#include <cmath>
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "Problem_Definition.h"
#include "Raschet.h"
using namespace std;

// описание конструктора 
Problem_Definition::Problem_Definition(
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
	double B[],
	double H[],
	double xU[],
	double yU[],
	double ForceX[],
	double ForceY[],
	double PhiX[],
	double PhiY[]
)
{
	Problem_Definition::Nx = Nx;
	Problem_Definition::Ny = Ny;
	Problem_Definition::x0 = x0;
	Problem_Definition::xN = xN;
	Problem_Definition::y0 = y0;
	Problem_Definition::yN = yN;
	Problem_Definition::lat0 = lat0;
	Problem_Definition::latN = latN;
	Problem_Definition::lon0 = lon0;
	Problem_Definition::lonN = lonN;
	Problem_Definition::T_begin = T_begin;
	Problem_Definition::T_end = T_end;

	Problem_Definition::B = new double[Nx*Ny];
	Problem_Definition::H = new double[Nx*Ny];
	Problem_Definition::xU = new double[Nx*Ny];
	Problem_Definition::yU = new double[Nx*Ny];
	Problem_Definition::ForceX = new double[Nx*Ny];
	Problem_Definition::ForceY = new double[Nx*Ny];
	Problem_Definition::PhiX = new double[Nx*Ny];
	Problem_Definition::PhiY = new double[Nx*Ny];
	Problem_Definition::xJ = new double[Nx*Ny]();
	Problem_Definition::yJ = new double[Nx*Ny]();
	Problem_Definition::dryFacesX = new int[Nx*Ny]();
	Problem_Definition::dryFacesY = new int[Nx*Ny]();
	//Problem_Definition::dT_ = new double[Nx*Ny]();
	if (!parallelOpenMP)
	{
		omp_set_num_threads(1);
	}
	#pragma omp parallel
	{
		Problem_Definition::OMP_THREADS_NUMBER = omp_get_num_threads();
	}	

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			Problem_Definition::B[i*Ny + j] = B[i*Ny + j];
			Problem_Definition::H[i*Ny + j] = H[i*Ny + j];
			Problem_Definition::xU[i*Ny + j] = xU[i*Ny + j];
			Problem_Definition::yU[i*Ny + j] = yU[i*Ny + j];
			Problem_Definition::ForceX[i*Ny + j] = ForceX[i*Ny + j];
			Problem_Definition::ForceY[i*Ny + j] = ForceY[i*Ny + j];
			Problem_Definition::PhiX[i*Ny + j] = PhiX[i*Ny + j];
			Problem_Definition::PhiY[i*Ny + j] = PhiY[i*Ny + j];
		}
	}

}

// пустой конструктор
Problem_Definition::Problem_Definition(){}

// деструктор
Problem_Definition::~Problem_Definition()
{
	delete[] B;
	delete[] H;
	delete[] xU;
	delete[] yU;
	delete[] ForceX;
	delete[] ForceY;
	delete[] PhiX;
	delete[] PhiY;

	delete[] xJ;
	delete[] yJ;

	delete[] dryFacesX;
	delete[] dryFacesY;
	//delete[] dT_;
}

// !!! функция считывания данных из текстового файла
void Problem_Definition::Continue_from(char Filename[80])
{
	ifstream  FILE;
	FILE.open(Filename, ios::out);

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			FILE >> X[i*Ny + j];
			FILE >> Y[i*Ny + j];
			FILE >> H[i*Ny + j];
			FILE >> xU[i*Ny + j];
			FILE >> yU[i*Ny + j];
			FILE >> ForceX[i*Ny + j];
			FILE >> ForceY[i*Ny + j];
			FILE >> B[i*Ny + j];
			//FILE >> epsilon.V[i*Ny+j];
		}
	}
	FILE.close();
}