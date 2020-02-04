/*
�������� ���� �������� ������.
*/

#define _USE_MATH_DEFINES
#define _XOPEN_SOURCE 600

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

//#include "Variable.h" // ������������ ���� �������� ������ �������� �������� ��� ���������� �����, �� ������������ 
#include "MainCore/Problem_Defenition.h"
#include "MainCore/Raschet.h"
#include "MainCore/technical.h"

#include <Windows.h>
//#include <mmsystem.h>

//#include <set>

#include <ctime>    
using namespace std;

double Bmin;



int main() {
	/* === Problem Definition === */

	// T_begin, T_end - start and end time respectively, in seconds 
	double T_begin = 0;
	double T_end = T_begin + 20;

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

					// ���������� ��������� 
	double mu = 0;// 0.0026; // ����������� ������ 
	double fc = 0;// 2 * 7.2921 / (100000) * sin((lat0 + latN)* 0.5 * M_PI / 180); // ����������� ���������

				  /* === ����������� ������������ ���������� ������� === */

	double beta = 0.1; // ����� ������� - ������(0; 1)
	double alpha = 0.1; // �������� �����������, ���������� �� ������� �������� � ������������ �����(����������������)
	double eps = 0.05;// pow(10, -6);// 0.0001; // �������� ������ ���, ������� 1 ���������	

	double NS = 1.0; // ����������� ��� ������� �����-������

	int Nx = 5; // ��������� �� ���������� x
	int Ny = 5; // ��������� �� ���������� y

				 /* === ����������� ��������� === */

				 // �������� ����� � \Data ����������� ���: Test_namePostscript
	string Test_name = "Dry_Zone_Problem"; // �������� ����� 
	string Postscript = "_v0_" + to_string(Nx) + "x" + to_string(Ny); // ��� �������� ����, ��� ������� �� ���������

																	  // ��������� ���������
	int Visualization_to_techplot_flag = 1; // ����� ����������� ��� ������������ � Tecplot
	double t_step = (T_end - T_begin) / 2; // �������� ������ ������ � ���� 
	double t_graph_export = T_begin;  // ������ �������, � �������� ������� ����� � ����

	double sea_level = 1000;  // ������ ������, ������� ������ � ������ (������� ������ ��������� ���������� � �� �������)

							  /* === ��������� ������ === */

	double* B = new double[Nx*Ny](); // ���
	double* H = new double[Nx*Ny](); // ������ ����
	double* xU = new double[Nx*Ny](); // �������� �� x
	double* yU = new double[Nx*Ny](); // �������� �� y

	double* ForceX = new double[Nx*Ny](); // ������� �������� ���� �� x
	double* ForceY = new double[Nx*Ny](); // ������� �������� ���� �� y
	double* PhiX = new double[Nx*Ny]();// ������� ������������� ���� �� x
	double* PhiY = new double[Nx*Ny]();// ������� ������������� ���� �� y

	//int sea_points = 0;  //!!! ��� ���? 
	//int level_points = 0;  //!!! ��� ���? 

	/* === ������������� ������ === */



	// ���� ������ �����	
	for (int i = 0; i < Nx; i++)
	{
		double x = x0 + i*(xN - x0) / (Nx - 1);
		for (int j = 0; j < Ny; j++)
		{
			double y = y0 + j*(yN - y0) / (Ny - 1);
			int k = i*Ny + j;
			B[k] = 10.0*(x - 50.0)*(x - 50.0) / (50.0*30.0) + 10.0*(y - 50.0)*(y - 50.0) / (50.0*30.0);
			if (y <= 35 && y >= 30)
				B[k] += 5;
			H[k] = 5.0 - B[k];

			if ((x - 50.0)*(x - 50.0) + (y - 50.0)*(y - 50.0) < 10.0 * 10.0)
				H[k] = 20.0;
			if (H[k] < 0)
				H[k] = 0.0;
		}
	}


	// ���� �� �����
	// !!!
	
	/* === ���������� ������� === */
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
	
	//R->SetFixedBoundaryConditions(VELOCITY_X, RIGHT, -5);
	//R->SetWallBoundaryConditions(TOP, BOTTOM);
	//R->SetFixedBoundaryConditions(HEIGHT, RIGHT, 0.5);
	//R->SetFixedBoundaryConditions(VELOCITY_X, RIGHT, -5);

	// boundary conditions setting, OPEN BOUNDARY is default
	// LOL
	//R->SetWallBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT);

	//system("pause");
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


	R->Exec_Raschet(); // ���������� �������

					   /* === ������� ������ === */

	delete B;
	delete H;
	delete xU;
	delete yU;
	delete ForceX;
	delete ForceY;
	delete PhiX;
	delete PhiY;

	system("pause");
	return 0;
}