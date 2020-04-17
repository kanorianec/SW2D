#pragma once /* Защита от двойного подключения заголовочного файла */

// !!!заголовочный файл описания класса, работающего с задачей, но непонятен выбор величин для работы с ним

#include <string>

using namespace std;

class Problem_Definition {

public:

	// Basic variables
	double* H; // water height (depth)
	double* xU; // horizontal velocity in x direction
	double* yU; // horizontal velocity in y direction
	double* B; // batymetry	
	double* ForceX; // external volume forces in x direction
	double* ForceY; // external volume forces in y direction
	double* PhiX; // external surface forces in x direction
	double* PhiY; // external surface forces in x direction

	// mass fluxes
	double* xJ;
	double* yJ;
	//double* dT_;

	// Additional variables
	double* C; // concentration of pollutant in transport equation

	// Area describing: X = LATITUDE, Y = LONGITUDE

	// area in Cartesian system
	double x0; // left x value
	double xN; // right x value

	double y0; // bottom y value
	double yN; // top y value

	// area in polar system in terms of Lat/Lon
	double lat0; // left latitude value
	double latN; // right latitude value

	double lon0; // bottom longitude value
	double lonN; // top longitude value

	int Nx; // Number of grid partitions in X direction
 	int Ny; // Number of grid partitions in Y direction

	double T_begin; // Time of start (in seconds)
	double T_end; // Time of end (in seconds)

	double* X; // array of X coordinates
	double* Y; // array of Y coordinates

	double* Lat; // array of Latitude coordinates
	double* Lon; // array of Longitude coordinates

	int PARALLEL_FLAG;
	int OMP_THREADS_NUMBER;

	// конструктор
	Problem_Definition();
	Problem_Definition(double x0, 
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
	);
	// Деструктор 
	~Problem_Definition();

	// !!! функция считывания данных из текстового файла
	void Continue_from(char[80]);

};