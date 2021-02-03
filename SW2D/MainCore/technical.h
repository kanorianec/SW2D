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
bool checkSymmetry(double* A, int Nx, int Ny, std::string name = "noname");

template <typename Temp>
void printTArray(Temp* A, int Nx, int Ny, std::string name);
void printArray(double* A, int Nx, int Ny, std::string name = "noname");
void printArray(int* A, int Nx, int Ny, std::string name = "noname");

template <typename Temp>
void printTFlux(Temp* Ax, Temp* Ay, int Nx, int Ny, std::string name);
void printFlux(double* Ax, double* Ay, int Nx, int Ny, std::string name = "noname");
void printFlux(int* Ax, int* Ay, int Nx, int Ny, std::string name = "noname");

bool checkEquality(double* A1, double* A2, int Nx, int Ny);


//int check_dot(Stations *S1, Stations *S2, int i, int j);