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
#include "technical.h"
#include "Raschet.h"

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

bool folderExists(string folder_path)
{
	struct stat buffer;
	return stat(folder_path.c_str(), &buffer);
}

// Creating folder
void Raschet::Prepare_Folder(string folder_path, bool ignore_warning)
{
	struct stat buffer;
	if (folderExists(folder_path))
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
/*
// функция сохранения текущих данных
void Raschet::Save_Data(double Time_of_work) {
	string save_path = path + "/" + to_string((int)T_begin) + "-" + to_string((int)T_end); // путь папки сохраненных данных

	Prepare_Folder(save_path);

	string name_h = save_path + "/H.dat";
	string name_xU = save_path + "/xU.dat";
	string name_yU = save_path + "/yU.dat";
	string name_C = save_path + "/C.dat";

	FILE *fH = fopen(name_h.c_str(), "w");
	FILE *fxU = fopen(name_xU.c_str(), "w");
	FILE *fyU = fopen(name_yU.c_str(), "w");
	FILE *fC;

	if (TransportProblemFlag)
	{
		fC = fopen(name_C.c_str(), "w");
	}

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			int k = i*Ny + j;
			fprintf(fH, "%lf ", H[k]);
			fprintf(fxU, "%lf ", xU[k]);
			fprintf(fyU, "%lf ", yU[k]);
			if (TransportProblemFlag)
			{
				fprintf(fC, "%lf ", C[k]);
			}
		}
	}

	fclose(fH);
	fclose(fxU);
	fclose(fyU);
	if (TransportProblemFlag)
	{
		fclose(fC);
	}

	//write_extra_inf_to_file(Time_of_work);
}

// функция загрузки величин из файла
void Raschet::Read_Data_from_file(string path_name) {

	string load_path = "Data/" + path_name;

	string name_h = load_path + "/H.dat";
	string name_xU = load_path + "/xU.dat";
	string name_yU = load_path + "/yU.dat";
	string name_C = load_path + "/C.dat";

	FILE *fH = fopen(name_h.c_str(), "r");
	FILE *fxU = fopen(name_xU.c_str(), "r");
	FILE *fyU = fopen(name_yU.c_str(), "r");
	FILE *fC;

	if (TransportProblemFlag)
	{
		fC = fopen(name_C.c_str(), "r");
	}

	if (fH && fxU && fyU)
	{
		for (int i = 0; i<Nx; i++)
		{
			for (int j = 0; j<Ny; j++)
			{
				int k = i*Ny + j;
				fscanf(fH, "%lf", &H[k]);
				fscanf(fxU, "%lf", &xU[k]);
				fscanf(fyU, "%lf", &yU[k]);
				if (TransportProblemFlag)
				{
					fscanf(fC, "%lf ", &C[k]);
				}
			}
		}
	}
	else {
		cout << "FILE " << path_name << " is not found!" << endl;
	}

	fclose(fH);
	fclose(fxU);
	fclose(fyU);
	if (TransportProblemFlag)
	{
		fclose(fC);
	}
};
*/
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
		cout << "Time = " << Time_elapsed + dT << endl;
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
		fprintf(F, "%lf %lf %lf %lf %lf %lf\n", Time_elapsed, H[i*Ny + j] + B[i*Ny + j] - Bmax, xU[i*Ny + j], yU[i*Ny + j], PhiXt[i*Ny + j], PhiYt[i*Ny + j]);
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
	out << "Combine file and free boundary conditions: " << COMBINE_FILE_AND_FREE_BOUNDARY_CONDITIONS << ";"<<endl;
	out << "    LEFT-TOP:           TOP:                RIGHT-TOP:" << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LT_CORNER], border_C[HEIGHT][LT_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][TOP], border_C[HEIGHT][TOP]) << GetConditionName(HEIGHT, border[HEIGHT][RT_CORNER], border_C[HEIGHT][RT_CORNER]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LT_CORNER], border_C[VELOCITY_X][LT_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][TOP], border_C[VELOCITY_X][TOP]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RT_CORNER], border_C[VELOCITY_X][RT_CORNER]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LT_CORNER], border_C[VELOCITY_Y][LT_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][TOP], border_C[VELOCITY_Y][TOP]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RT_CORNER], border_C[VELOCITY_Y][RT_CORNER]) << endl;
	if (TransportProblemFlag) out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LT_CORNER], border_C[CONCENTRATION][LT_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][TOP], border_C[CONCENTRATION][TOP]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RT_CORNER], border_C[CONCENTRATION][RT_CORNER]) << endl;
	out << "                   \\ ___________________ /" << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "    LEFT:           |                   |    RIGHT:" << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LEFT], border_C[HEIGHT][LEFT]) << "|                   |" << GetConditionName(HEIGHT, border[HEIGHT][RIGHT], border_C[HEIGHT][RIGHT]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LEFT], border_C[VELOCITY_X][LEFT]) << "|                   |" << GetConditionName(VELOCITY_X, border[VELOCITY_X][RIGHT], border_C[VELOCITY_X][RIGHT]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LEFT], border_C[VELOCITY_Y][LEFT]) << "|                   |" << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RIGHT], border_C[VELOCITY_Y][RIGHT]) << endl;
	if (TransportProblemFlag)
		out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LEFT], border_C[CONCENTRATION][LEFT]) << "|                   |" << GetConditionName(CONCENTRATION, border[CONCENTRATION][RIGHT], border_C[CONCENTRATION][RIGHT]) << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "                    |                   |" << endl;
	out << "                    |___________________|" << endl;
	out << "                   /                     \\" << endl;
	out << "    LEFT-BOTTOM:        BOTTOM:             RIGHT-BOTTOM:" << endl;
	out << GetConditionName(HEIGHT, border[HEIGHT][LB_CORNER], border_C[HEIGHT][LB_CORNER]) << GetConditionName(HEIGHT, border[HEIGHT][BOTTOM], border_C[HEIGHT][BOTTOM]) << GetConditionName(HEIGHT, border[HEIGHT][RB_CORNER], border_C[HEIGHT][RB_CORNER]) << endl;
	out << GetConditionName(VELOCITY_X, border[VELOCITY_X][LB_CORNER], border_C[VELOCITY_X][LB_CORNER]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][BOTTOM], border_C[VELOCITY_X][BOTTOM]) << GetConditionName(VELOCITY_X, border[VELOCITY_X][RB_CORNER], border_C[VELOCITY_X][RB_CORNER]) << endl;
	out << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][LB_CORNER], border_C[VELOCITY_Y][LB_CORNER]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][BOTTOM], border_C[VELOCITY_Y][BOTTOM]) << GetConditionName(VELOCITY_Y, border[VELOCITY_Y][RB_CORNER], border_C[VELOCITY_Y][RB_CORNER]) << endl;
	if (TransportProblemFlag) out << GetConditionName(CONCENTRATION, border[CONCENTRATION][LB_CORNER], border_C[CONCENTRATION][LB_CORNER]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][BOTTOM], border_C[CONCENTRATION][BOTTOM]) << GetConditionName(CONCENTRATION, border[CONCENTRATION][RB_CORNER], border_C[CONCENTRATION][RB_CORNER]) << endl;

	out << endl << "======== TECHNICAL ========" << endl;
	if (OMP_THREADS_NUMBER == 1)
		out << "OpenMP is off. " << endl;
	else
		out << "OpenMP is on. Program uses " << OMP_THREADS_NUMBER << " threads." << endl;
	if (Time_of_work > 0.0)
	    out << "Time of work = " << Time_of_work << " seconds;" << endl;
	out << "include forcing to regularization: F_reg = " << F_reg << "; Phi_reg" << Phi_reg << endl;
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
	
	int CurLen = Name.length();

	for (int i = 1; i <= StrLen - CurLen; i++)
		Name += " ";
	//cout << Name.length() << endl;
	//system("pause");
	return Name;
}

std::string to_str(double num)
{
	std::string res = to_string(num);

	while (res.back() == '0')
	{
		res.pop_back();
	}
	if (res.back() == '.')
		res.pop_back();
	return res;
}

std::string to_str(int num)
{
	return to_string(num);
}