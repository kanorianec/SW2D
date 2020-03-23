#pragma once
#include <string>

#include "Problem_Definition.h"
#include "../ephemeris/Ephemeris.h"
//#include <time.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include "Constants.h"

class Raschet : public Problem_Definition
{
public:
	std::string Test_name; //  Name of problem
	std::string Postscript; // additional postscript
	std::string path; // path to the test directory

	int Polar; // flag of Polar(1) or Cartesian(0) is main coordinate system

	double hx; // grid step by x
	double hy; // grid step by y

	double hlat; // grid step by latitude
	double hlon; // grid step by longitude

	double alpha; // tuning parameter in range (0,1)
	double beta; // CFL number (0, 1)
	double NS; // coefficient at the Navier-Stokes tensor
	
	// Forcing:
	double mu; // bottom friction coefficient
	double windFrictionCoef; // wind friction coefficient
	int fc; // use Coriolis force (1) or not (0)

	// Transport problem:
	double D; // Diffusion coefficient
	

	double Hmax; // maximal depth
	double Bmax; // maximal value of bathymetry
	double Bmin; // minimal value of batymetry
	double dT; // time step
	double Time_elapsed; // elapsed time

	double eps; // dry-zones parameter
	int* epsilon; // flags for dry zones detecting

	// arrays of data in the next time step k + 1
	double* tau;
	double* xUt;
	double* yUt;
	double* Ht;
	//double* PhiXt;
	//double* PhiYt;

	double* xWind;
	double* yWind;

	double* Ct;
	bool TransportProblemFlag;

	// type of border conditions, coeff: 1 (free), -1 (wall or constant), 0 (from_file); coef_c - constant value, only with
	int border[4][8];
	double border_C[4][8];

	// boundary conditions from file
	bool BoundaryConditionsFromFile;
	bool InternalWallsFlag;
	double* lin_b[3][4];
	double* lin_k[3][4];
	double t1_bound; // first time moment in seconds
	double t2_bound;
	FILE* FV[3][4]; // 3 vars x 4 walls: arrray of pointers to boundary conditions file 
	FILE* FT; // arrray of pointers to time shifts file

	bool windForcing;
	std::ifstream FWindX; // xU arrray of pointers to wind speed file
	std::ifstream FWindY; // yU arrray of pointers to wind speed file
	double timeWind;
	double timeWindPeriod;

	// ephemeris for tides calculation
	double Sbeta; // declination of Sun in rad
	double Sa; // right ascension of Sun in rad
	double ST; // apparent Sidereal Time in rad
	double SDist; // Distance between Earth and Sun
	double MDist; // Distance between Eath and Moon
	double Mbeta; // declination of Moon in rad
	double Ma; // right ascension of Moon in rad

	int* S;//TypeOfPoint* S; // type of the pount
	//double X_spot, Y_spot;

	int Stop_Raschet_Flag; //stop

	
	int HourMark;
	double sea_level; // shore height, which is included in the calculation

	bool Visualization_to_techplot_flag; // output flag to Tecplot 

	int i_minVis;
	int i_maxVis;
	int j_minVis;
	int j_maxVis;

	double t_step;
	double t_graph_export;
	double t_end_graph_export;
	bool first_visualization;
	int output_per_file_counter;
	double current_file_sizeMB;

	bool restart;

	time_t RaschetTime; // Date and time (UTC) in seconds from 1900 year of problem start 

	//Raschet();

	// Constructor for ordinary Catesian area
	Raschet(std::string Test_name,
		std::string Postscript,
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
	);

	// Constructor for Polar area
	Raschet(std::string Test_name,
		std::string Postscript,
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
	);

	void Prepare_Raschet();	// Preparing calculations
	void Exec_Raschet(); // Calculations process execution 
	void Perform_Calculations(); // Preforming calculations

	void Initialize_Transport_Problem(double InitialC[],
		double D
	);
	
	void outputResults();
	void outputInputs();
	void Visualization_to_techplot_input(); // визуализация данных в формате, используемом для построения в Techplot
	void Visualization_to_techplot_result();
	
	//void Visualization_to_techplot_TECIO(int start=0);
	// void Visualization_to_techplot(); // старая реализация была сделана через fstream. Оказалась медленной, удалена. 

	void Numerical_scheme_time_step_parallel(); // параллельная реализация(через openMP) функции численного счета
	//void Numerical_scheme_time_step(); // убрана, поскольку используется параллельная реализация

	void Prepare_Folder(std::string folder_path, bool ignore_warning = true); // Function to create folder

	void Save_Data(); // Save all variable at time moment "Time_elapsed"
	void Restart_from_time_moment(double Time_moment); // Reastart: reading all variables from time moment "Time_moment"
	void Save_Grid(); // function for saving grid data: B, X, Y, ...

	void Print_info_about_point(std::string name, int index); // вывод информации о расчетной точке
	void Write_point_to_file(int index, double X_cord, double Y_Cord, std::string file_name); // вывод информации о точке в файл, ЛУЧШЕ ОПТИМИЗИРОВАТЬ: ВЫНЕСТИ УСЛОВИЕ В БЛОК ПРОГРАММЫ
	void Write_point_to_file(int index, std::string file_name);
	void write_extra_inf_to_file(double Time_of_work); // Writes test description to file "extra_inf.txt"
	void write_extra_inf(std::ostream &out, double Time_of_work); // Writes test description to file "extra_inf.txt"
	
	void SetVisualizationProperties(double Vis_T_start, double Vis_T_end, int iMin, int jMin, int iMax, int jMax);
	void SetStartTime(int Year, int Month, int Day, int Hour, int Minute, int Second); // Function to set start time (UTC) in terms of time_t variable
	
	void SetFixedBoundaryConditions(TypeOfVariable VType, TypeOfPoint PType, double Value);
	
	// Set up boundary conditions of d/dn = 0
	void SetZeroDerivativeConditions(TypeOfVariable VType, TypeOfPoint PType1, TypeOfPoint PType2 = EXCLUDED, TypeOfPoint PType3 = EXCLUDED, TypeOfPoint PType4 = EXCLUDED, TypeOfPoint PType5 = EXCLUDED, TypeOfPoint PType6 = EXCLUDED, TypeOfPoint PType7 = EXCLUDED, TypeOfPoint PType8 = EXCLUDED);

	// Wall boundary conditions: u_n = 0, du_tau/dn = 0, dh/dn = 0;
	void SetWallBoundaryConditions(TypeOfPoint PType1, TypeOfPoint PType2 = EXCLUDED, TypeOfPoint PType3 = EXCLUDED, TypeOfPoint PType4 = EXCLUDED);
	
	void SetInternalWall(TypeOfPoint PType, int solid_ind, int start_ind, int end_ind);
	// Open boundary conditions: du_n/dn = du_tau/dn = dh/dn = 0;
	void SetOpenBoundaryConditions(TypeOfPoint PType1, TypeOfPoint PType2 = EXCLUDED, TypeOfPoint PType3 = EXCLUDED, TypeOfPoint PType4 = EXCLUDED, TypeOfPoint PType5 = EXCLUDED, TypeOfPoint PType6 = EXCLUDED, TypeOfPoint PType7 = EXCLUDED, TypeOfPoint PType8 = EXCLUDED);
	// Boundary conditions from file
	void SetFileBoundaryConditions(TypeOfVariable VType, TypeOfPoint PType1, TypeOfPoint PType2 = EXCLUDED, TypeOfPoint PType3 = EXCLUDED, TypeOfPoint PType4 = EXCLUDED);
	void RecalcFileBoundaryConditions();

	// Setting wind speed to initialize wind friction force
	void SetWindSpeed(double WindFrictionCoefficient, double period);

	void Recalc_forces_parallel(); // Forces recalculation with OpenMP

	//void Recalc_forces_parallel();
	//void Operate_City_data(int i, int j, float X_cord, float Y_Cord, char FName[50]); // старая версия Write_point_to_file

	// Destructor
	~Raschet();

};

