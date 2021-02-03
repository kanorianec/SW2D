#pragma once

double polar_to_decart_x(double dlon, double lat);
double polar_to_decart_y(double dlat);
char* GetTimeStamp();
std::string GetConditionName(TypeOfVariable VType, int border, double border_C);
std::string GetVariableName(TypeOfVariable VType);
std::string to_str(int num);
std::string to_str(double num, int p = -1);
bool folderNotExists(std::string folder_path);

void pause();
void checkSymmetry(double* A, int Nx, int Ny, std::string name = "noname");
bool checkEquality(double* A1, double* A2, int Nx, int Ny);

//int check_dot(Stations *S1, Stations *S2, int i, int j);