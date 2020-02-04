#include <cmath>
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
//#include "Variable.h" // Заголовочный файл описания класса сеточных значений для разностной схемы, НЕ ИСПОЛЬЗУЕТСЯ 
#include "Problem_Defenition.h"
#include "Raschet.h"
using namespace std;

// описание конструктора 
Problem_Defenition::Problem_Defenition(
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
	Problem_Defenition::Nx = Nx;
	Problem_Defenition::Ny = Ny;
	Problem_Defenition::x0 = x0;
	Problem_Defenition::xN = xN;
	Problem_Defenition::y0 = y0;
	Problem_Defenition::yN = yN;
	Problem_Defenition::lat0 = lat0;
	Problem_Defenition::latN = latN;
	Problem_Defenition::lon0 = lon0;
	Problem_Defenition::lonN = lonN;
	Problem_Defenition::T_begin = T_begin;
	Problem_Defenition::T_end = T_end;

	Problem_Defenition::B = new double[Nx*Ny];
	Problem_Defenition::H = new double[Nx*Ny];
	Problem_Defenition::xU = new double[Nx*Ny];
	Problem_Defenition::yU = new double[Nx*Ny];
	Problem_Defenition::ForceX = new double[Nx*Ny];
	Problem_Defenition::ForceY = new double[Nx*Ny];
	Problem_Defenition::PhiX = new double[Nx*Ny];
	Problem_Defenition::PhiY = new double[Nx*Ny];

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			Problem_Defenition::B[i*Ny + j] = B[i*Ny + j];
			Problem_Defenition::H[i*Ny + j] = H[i*Ny + j];
			Problem_Defenition::xU[i*Ny + j] = xU[i*Ny + j];
			Problem_Defenition::yU[i*Ny + j] = yU[i*Ny + j];
			Problem_Defenition::ForceX[i*Ny + j] = ForceX[i*Ny + j];
			Problem_Defenition::ForceY[i*Ny + j] = ForceY[i*Ny + j];
			Problem_Defenition::PhiX[i*Ny + j] = PhiX[i*Ny + j];
			Problem_Defenition::PhiY[i*Ny + j] = PhiY[i*Ny + j];
		}
	}

}

// пустой конструктор
Problem_Defenition::Problem_Defenition(){}

// деструктор
Problem_Defenition::~Problem_Defenition()
{
	delete B;
	delete H;
	delete xU;
	delete yU;
	delete ForceX;
	delete ForceY;
	delete PhiX;
	delete PhiY;
}

// !!! функция считывания данных из текстового файла
void Problem_Defenition::Continue_from(char Filename[80])
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