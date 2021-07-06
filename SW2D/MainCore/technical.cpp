/*
Файл описания технических функций
*/
#pragma once

#define _USE_MATH_DEFINES
#define _XOPEN_SOURCE 600
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>
#include <iomanip>

#ifdef __linux__
#include <sys/types.h>
#include <sys/stat.h>
#endif // linux

#if defined(_WIN64) || defined(_WIN32)
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif // _WIN64) || defined(_WIN32)

using namespace std;

#include "Constants.h"
#include "Raschet.h"
#include <omp.h>
#include "technical.h"


//#include <stdio.h>
//#include <stdlib.h>

//#include <algorithm>

// вычисление x-координаты из географических координат 
double polar_to_decart_x(double dlon, double lat) {
	double x = length_ekv * cos(lat*rad) * fabs(dlon) / 360;
	return x;
}
// вычисление y-координаты из географических координат 
double polar_to_decart_y(double dlat) {
	double y = length_mer * fabs(dlat) / 360;
	return y;
}

bool folderNotExists(string folder_path)
{
	struct stat buffer;
	return stat(folder_path.c_str(), &buffer);
}

// Creating folder
void Raschet::Prepare_Folder(string folder_path, bool ignore_warning)
{
	struct stat buffer;
	if (folderNotExists(folder_path))
	{
			if (mkdir(folder_path.c_str(), 0755))
				cout << "Error: can't create folder " << folder_path << endl;
	}
	else
	{
		if (!ignore_warning)
		{
			cout << "Warning, " << folder_path << " already exists! Would you like to replace it? yes [y], no [n]:" << endl;
			char c;
			cin >> c;
			while (c != 'y' && c != 'n')
			{
				cin >> c;
				//cout << "!" << c << endl;
			}
			if (c == 'y')
				Prepare_Folder(folder_path);
			else
			{
				cin.ignore(1024, '\n');
				cout << "Press Enter to stop the programm." << endl;
				cin.get();
				exit(0);
			}
		}
	}
}

// 
void Raschet::Print_info_about_point(string name, int index) {
	int  i = int(index / Ny);
	int j = index % Ny;
	if (i <= Nx - 1 && j <= Ny - 1 && i >= 0 && j >= 0)
	{
		cout << "=============================================" << endl;
		cout << "Info about point " << name << endl;
		cout << "index = " << index << "; i = " << i << ", j = " << j << endl;
		cout << "X = " << X[index] << ", Y = " << Y[index] << ";" << endl;
		cout << "eps = " << epsilon[index] << ";" << endl;
		cout << "High = " << H[index] << "; B =  " << B[index] + Bmin << ";" << endl;
		cout << "Ux = " << xU[index] << "; Uy = " << yU[index] << ";" << endl;
		cout << "Utx = " << xUt[index] << "; Uty = " << yUt[index] << ";" << endl;
		cout << "Ht = " << Ht[index] << ";" << endl;
		if (TransportProblemFlag)
		{
			cout << "C = " << C[index] << "; Ct = " << Ct[index] << ";" << endl;
		}
		cout << "tau = " << tau[index] << "; " << endl;
		cout << "Type of point: " << S[index] << endl;
		cout << "Time = " << Time_elapsed + dT << "; Time index = " << (int)((Time_elapsed + dT)/dT)<< endl;
	}	
};

// 
void Raschet::Write_point_to_file(int index, double X_cord, double Y_Cord, string file_name)
{
	file_name = path + "/" + file_name;
	int  i = int(index / Ny);
	int j = index % Ny;
	if ((i*hx <= X_cord) && ((i + 1)*hx > X_cord) && (j*hy <= Y_Cord) && ((j + 1)*hy > Y_Cord))
	{
		FILE *F = fopen(file_name.c_str(),"a");
		fprintf(F, "%lf %lf %lf %lf %lf %lf\n", Time_elapsed, H[i*Ny + j] + B[i*Ny + j] - Bmax, xU[i*Ny + j], yU[i*Ny + j], PhiX[i*Ny + j], PhiY[i*Ny + j]);
		fclose(F);
	}
}

void Raschet::Write_point_to_file(int index, string file_name)
{
	file_name = path + "/" + file_name;
	FILE *F = fopen(file_name.c_str(), "a");
	fprintf(F, "%lf %lf %lf %lf\n", Time_elapsed, H[index] + B[index], xU[index], yU[index]);
	fclose(F);
}


// write all metadata of the solving problem
void Raschet::write_extra_inf_to_file(double Time_of_work /*, double tt*/)
{
	string file_name = "extra_inf.txt";
	ofstream file;
	file.precision(10);
	file.open(path + "/" + file_name);
	if (file.is_open())
		cout << path + "/" + file_name << endl;
	write_extra_inf(file,Time_of_work);
	file.close();
}

// write all metadata of the solving problem
void Raschet::write_extra_inf(ostream &out, double Time_of_work /*, double tt*/)
{
	out << GetTimeStamp();

	out << endl << "====== COEFFICIENTS ======" << endl;
	out << "Nx = " << Nx << "; Ny = " << Ny << ";" << endl;
	out << "alpha = " << alpha << "; beta = " << beta << ";" << endl;
	out  << "NS = " << NS << "; eps0 = " << eps << ";" << endl;

	out << endl << "====== FORCING ======" << endl;
	out << "fc = " << fc << "; mu = " << mu << ";" << endl;
	out << "Tide forcing = " << TideForcing << ";" << endl;
	out << "Wind forcing = " << windForcing << ";" << endl;
	if (TransportProblemFlag) out << "Diffusion coefficient: D = " << D << ";" << endl;

	out << endl << "======== AREA & TIME ========" << endl;
	out << "dx = " << hx << "; dy = " << hy << ";" << endl;
	out << "x0 = " << X[0] << ";" << endl << "xN = " << X[Nx*Ny - 1] << ";" << endl;
	out << "y0 = " << Y[0] << ";" << endl << "yN = " << Y[Ny - 1] << ";" << endl;
	if (hlon + hlat > 0)
		out << "dlat = " << hlat << "; dlon =  " << hlon << ";" << endl;
	if (lonN + latN > 0)
	{
		out << "lat0 = " << Lat[0] << "; latN = " << Lat[Ny - 1] << ";" << endl;
		out << "lon0 = " << Lon[0] << "; lonN = " << Lon[Nx*Ny - 1] << ";" << endl;
	}
	out << "Bmin = " << Bmin << endl;
	out << "dt = " << dT << "; Number of iterations: " << (int)((T_end - T_begin)/dT) << ";" << endl;
	out << "Time_begin = " << T_begin << ";" << endl << "Time_End = " << T_end << ";" << endl;
	out << "Start Date and time of problem solving: " << asctime(gmtime(&RaschetTime));

	out << "======================= BOUNDARY CONDITIONS =========================" << endl;
	out << "Tide harmonics boundary conditions: " << tidesHarmonics << ";" << endl;
	out << "Combine file and free boundary conditions: " << COMBINE_FILE_AND_FREE_BOUNDARY_CONDITIONS << ";"<<endl;
	out << "    LEFT-TOP:           TOP:                RIGHT-TOP:" << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LT_CORNER], border_C[HEIGHT][LT_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][TOP], border_C[HEIGHT][TOP]) << GetConditionName(HEIGHT, border[HEIGHT][RT_CORNER], border_C[HEIGHT][RT_CORNER]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LT_CORNER], border_C[VELOCITY_X][LT_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][TOP], border_C[VELOCITY_X][TOP]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RT_CORNER], border_C[VELOCITY_X][RT_CORNER]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LT_CORNER], border_C[VELOCITY_Y][LT_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][TOP], border_C[VELOCITY_Y][TOP]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RT_CORNER], border_C[VELOCITY_Y][RT_CORNER]) << endl;
	if (TransportProblemFlag) out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LT_CORNER], border_C[CONCENTRATION][LT_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][TOP], border_C[CONCENTRATION][TOP]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RT_CORNER], border_C[CONCENTRATION][RT_CORNER]) << endl;
	out << "                     " << GetConditionName(TIDE, int(tideSide[TOP]), 0) << endl;
	out << "                   \\ ___________________ /" << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "    LEFT:           |                   |    RIGHT:" << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LEFT], border_C[HEIGHT][LEFT]) << "|                   |" << GetConditionName(HEIGHT, border[HEIGHT][RIGHT], border_C[HEIGHT][RIGHT]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LEFT], border_C[VELOCITY_X][LEFT]) << "|                   |" << GetConditionName(VELOCITY_X, border[VELOCITY_X][RIGHT], border_C[VELOCITY_X][RIGHT]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LEFT], border_C[VELOCITY_Y][LEFT]) << "|                   |" << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RIGHT], border_C[VELOCITY_Y][RIGHT]) << endl;
	if (TransportProblemFlag)
		out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LEFT], border_C[CONCENTRATION][LEFT]) << "|                   |" << GetConditionName(CONCENTRATION, border[CONCENTRATION][RIGHT], border_C[CONCENTRATION][RIGHT]) << endl;
	out << GetConditionName(TIDE, int(tideSide[LEFT]), 0) << "|                   |" << GetConditionName(TIDE, int(tideSide[RIGHT]), 0) << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "                    |___________________|" << endl;
	out << "                   /                     \\" << endl;
	out << "    LEFT-BOTTOM:        BOTTOM:             RIGHT-BOTTOM:" << endl;
	out << "                     " << GetConditionName(TIDE, int(tideSide[BOTTOM]), 0) << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LB_CORNER], border_C[HEIGHT][LB_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][BOTTOM], border_C[HEIGHT][BOTTOM]) << GetConditionName(HEIGHT, border[HEIGHT][RB_CORNER], border_C[HEIGHT][RB_CORNER]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LB_CORNER], border_C[VELOCITY_X][LB_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][BOTTOM], border_C[VELOCITY_X][BOTTOM]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RB_CORNER], border_C[VELOCITY_X][RB_CORNER]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LB_CORNER], border_C[VELOCITY_Y][LB_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][BOTTOM], border_C[VELOCITY_Y][BOTTOM]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RB_CORNER], border_C[VELOCITY_Y][RB_CORNER]) << endl;
	if (TransportProblemFlag) out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LB_CORNER], border_C[CONCENTRATION][LB_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][BOTTOM], border_C[CONCENTRATION][BOTTOM]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RB_CORNER], border_C[CONCENTRATION][RB_CORNER]) << endl;

	out << endl << "======== TECHNICAL ========" << endl;
	out << "massFluxCorrection flag = " << massFluxCorrection << endl;
	int thrnum = 0;
	#pragma omp parallel
	{
		thrnum = omp_get_num_threads();
	}
	if (thrnum == 1)
		out << "OpenMP is off. " << endl;
	else
		out << "OpenMP is on. Program uses " << thrnum << " threads." << endl;
	if (Time_of_work > 0.0)
	    out << "Time of work = " << Time_of_work << " seconds;" << endl;
	out << "include forcing to regularization: F_reg = " << F_reg << "; Phi_reg = " << Phi_reg << endl;
	if (TransportProblemFlag)
	{
		out << "type of transport equation regularization: alpha_c = 1 - normal, alpha_c = 0 - simplified, also could be between (0,1). alpha_c = " << alpha_c << ";" << endl;
		out << "Coefficient of viscosity in the transport equation, basic = 0.0, for special cases = 1.0/gc; NSC = " << NSC << ";" << endl;
	}
}

inline char* GetTimeStamp()
{
	time_t rawtime = time(NULL);
	struct tm * timeinfo = localtime(&rawtime);
	return asctime(timeinfo);
}

std::string GetVariableName(TypeOfVariable VType)
{
	std::string Var;

	switch (VType)
	{
	case HEIGHT:
		Var = "H";
		break;
	case VELOCITY_X:
		Var = "Ux";
		break;
	case VELOCITY_Y:
		Var = "Uy";
		break;
	case CONCENTRATION:
		Var = "C";
		break;
	case TIDE:
		Var = "Tide";
	}

	return Var;
}

std::string GetConditionName(TypeOfVariable VType, int border, double border_C)
{
	std::string Var;
	std::string Name = "   ";

	const int StrLen = 20;

	Var = GetVariableName(VType);

	std::ostringstream streamObj;
	streamObj << border_C;

	if (border == 1 && border_C == 0.0)
		Name += "d" + Var + "/dn = 0";
	else if (border == -1)
		Name += Var + " = " + streamObj.str();
	else if (border == FROM_FILE)
		Name += Var + ": file " + GetVariableName(VType) + "_?.dat";
	else
		Name += Var + ": unknown" + to_string(border) + " " +  streamObj.str();

	if (VType == TIDE)
	{
		Name = "   ";
		if (border)
			Name += "~~~Tide~~~";
	}
	
	
	int CurLen = Name.length();

	for (int i = 1; i <= StrLen - CurLen; i++)
		Name += " ";
	//cout << Name.length() << endl;
	//pause();
	return Name;
}

std::string to_str(double num, int p)
{
	if (p < 0)
	{
		std::string res = to_string(num);

		while (res.back() == '0' && res != "0")
		{
			res.pop_back();
		}
		if (res.back() == '.')
			res.pop_back();
		return res;
	}
	else
	{
		std::ostringstream out;
		out.precision(p);
		out << std::fixed << num;
		return out.str();
	}
}

std::string to_str(int num)
{
	return to_string(num);
}

void pause()
{
	cout << "Pause. Press Enter to continue." << endl;
	cin.get();
}

bool checkSymmetry(double* A, int Nx, int Ny, string name)
{
	bool Sym = true;
	// horizontal symmetry
	bool hSym = true;
	bool ahSym = true;
	bool vSym = true;
	bool avSym = true;

	cout << " == Array \"" << name << "\": ==" << endl;
	cout.precision(10);

	double difZero = 0;// 1e-16;

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{		
			if (i <= Nx / 2)
			{
				hSym = hSym && (fabs(A[i*Ny + j] - A[(Nx - i - 1)*Ny + j]) <= difZero); //(A[i*Ny + j] == A[(Nx - i - 1)*Ny + j]);
				ahSym = ahSym && (fabs(A[i*Ny + j] + A[(Nx - i - 1)*Ny + j]) <= difZero); //(A[i*Ny + j] == - A[(Nx - i - 1)*Ny + j]);
				if (!(hSym || ahSym))
				{
					cout << "h";
					//DEBUG
					//cout << fabs(A[i*Ny + j] - A[(Nx - i - 1)*Ny + j]) << " " << fabs(A[i*Ny + j] + A[(Nx - i - 1)*Ny + j]) << endl;
					//system("pause");
				}
					
			}
				
			if (j <= Ny / 2)
			{
				vSym = vSym && (fabs(A[i*Ny + j] - A[i*Ny + Ny - j - 1]) <= difZero); //(A[i*Ny + j] == A[i*Ny + Ny - j - 1]);
				avSym = avSym && (fabs(A[i*Ny + j] + A[i*Ny + Ny - j - 1]) <= difZero); //(A[i*Ny + j] == - A[i*Ny + Ny - j - 1]);
				if (!(vSym || avSym))
				{
					cout << "v";
					//DEBUG
					//cout << fabs(A[i*Ny + j] - A[i*Ny + Ny - j - 1]) << " " << fabs(A[i*Ny + j] + A[i*Ny + Ny - j - 1]) << endl;
					//system("pause");
				}
					
			}
			cout << A[i*Ny + j] << " ";
		}
		cout << endl;
	}

	cout << "== ==" << endl;

	if (hSym)
		cout << "Array \"" << name << "\" is symmetrical horizontally: >-<" << endl;
	else if (ahSym)
		cout << "Array \"" << name << "\" is symmetrical horizontally up to a sign: +>-<-" << endl;
	else
		cout << "Array \"" << name << "\" is NOT symmetrical horizontally!" << endl;
	if (vSym)
		cout << "Array \"" << name << "\" is symmetrical vertically: >|< " << endl;
	else if (avSym)
		cout << "Array \"" << name << "\" is symmetrical vertically up to a sign: +>|<- " << endl;
	else
		cout << "Array \"" << name << "\" is NOT symmetrical vertically!" << endl;
	Sym = Sym && (hSym || ahSym) && (vSym || avSym);
	return Sym;
}
bool checkEquality(double* A1, double* A2, int Nx, int Ny)
{
	bool eq = true;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			eq = eq && (A1[i*Ny + j] == A2[i*Ny + j]);
			if (!eq)
				cout << "Equality error in (" << i << ", " << j << ");" << endl;
		}
	}
	return eq;
}

void printArray(std::ostream& os, double* A, int Nx, int Ny, string name)
{
	printTArray(os, A, Nx, Ny, name);
}
void printArray(std::ostream& os, int* A, int Nx, int Ny, string name)
{
	printTArray(os, A, Nx, Ny, name);
}

template <typename Temp>
void printTArray(std::ostream& os, Temp* A, int Nx, int Ny, string name)
{
	os << " == Array \"" << name << "\": ==" << endl;
	os.precision(4);
	for (int j = 0; j < Ny; j++)
	{
		if (j == 0)
		{
			os << left << setw(10) << " /";
			for (int i = 0; i < Nx; i++)
			{
				if (j != Ny - 1)
					os << setfill('-');
				os << left << setw(10)  << to_str(i) + " " << setfill(' ');
			}
			os << endl << endl;
			setfill(' ');
		}

		for (int i = 0; i < Nx; i++)
		{
			if (i == 0)
				os << left << setw(10) << " " + to_str(j);
			os << left << setw(10) << to_str(A[i*Ny + j]);
		}

		os << endl << endl;
	}
	os << endl;
}


template <typename Temp>
void printTFlux(std::ostream& os, Temp* Ax, Temp* Ay, int Nx, int Ny, string name)
{
	os << " == Flux \"" << name << "\": ==" << endl;

	int p = 4;
	//int interval = 5;

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			//os << left << setw(10) << A[i*Ny + j] << " ";
			string Xdirection = " ";
			double xAval = 0.0; 

			if (j!=0)
				xAval = Ax[(i - 1)*Ny + j];

			string Axstr = to_str(fabs(xAval), p);

			if (xAval == 0.0)
				Axstr = string((p + 2) / 2-1, ' ') + "o" + string((p + 2) / 2, ' ');
			if (xAval < 0)
				Xdirection = "<";
			else if (xAval > 0)
				Xdirection = ">";
			
			if (i != 0 && j != 0 && j != Ny - 1)
				os << string((p + 2) / 2 + 1, ' ') << Xdirection + Axstr << string((p + 2) / 2 + 1, ' ') + "x";
			else if (i == 0)
				os << "  " + to_str(j);
			else
				os << " " + string(2 * (p + 3) - 1, '-') + " " + to_str(i);
		}
		os << endl << endl;
		if (j != Ny - 1)
			for (int i = 0; i < Nx; i++)
			{
				string Ydirection = " ";
				double yAval = 0.0;

				if (i != 0)
					yAval = Ay[i*Ny + j];

				string Aystr = to_str(fabs(yAval), p);

				if (yAval == 0.0)
					Aystr = string((p + 2) / 2 - 1, ' ') + "o" + string((p + 2) / 2, ' ');
				if (yAval < 0)
					Ydirection = "^";
				else if (yAval > 0)
					Ydirection = "v";
				//if (j == 0)
				//	os << left << setw(interval * 2) << " |";
				//else
				//	os << internal << setw(interval * 2) << "  |  ";

				if (i != Nx - 1 && i != 0)
					os << string((p + 2) / 2 + 1, ' ') << "+" << string((p + 2) / 2 + 1, ' ') + Ydirection + Aystr;
				//cout << string((p + 2) / 2, ' ') + "x" + string((p + 2) / 2, ' ') << Ydirection + Aystr;			
				else if (i == Nx - 1)
					os << string((p + 2) / 2 + 1, ' ') << "+" << string((p + 2) + 1, ' ') + "|";
				else if (i == 0)
					os << "  " + (string)"|" + string((p + 2) / 2, ' ');
			}
		os << endl << endl;
	}
	os << endl;
}

void printFlux(std::ostream& os, double* Ax, double* Ay, int Nx, int Ny, string name)
{
	printTFlux(os, Ax, Ay, Nx, Ny, name);
}

void printFlux(std::ostream& os, int* Ax, int* Ay, int Nx, int Ny, string name)
{
	printTFlux(os, Ax, Ay, Nx, Ny, name);
}

double sqr(double x)
{
	return x*x;
}