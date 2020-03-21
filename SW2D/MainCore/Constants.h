#pragma once
//================= Technical Constants ====================
// include forcing to regularization
const int F_reg = 0;
const int Phi_reg = 0;

const bool COMBINE_FILE_AND_FREE_BOUNDARY_CONDITIONS = false; //true;

// include forcing in bpundary conditions
const int F_bound = 0;

// type of transport equation regularization: alpha_c = 1 - normal, alpha_c = 0 - simplified, also could be between (0,1).
const double alpha_c = 1.0;
// Coefficient of viscosity in the transport equation, basic = 0.0, for special cases = 1.0/gc
const double NSC = 0.0;
//=========================================================

//const float pi = 3.1415926535;
//const int Nmax = 1555000; // !? 
const double gc = 9.81;
//const double G = 0.0000000000667408; // gravitational constant
const double CriticalVal = 10000.0;

const double SD = 1.2035; // Dudson's constant for Sun
const double MD = 2.6206; // Dudson's constant for Moon
const int TideForcing = 1; //
const int CoriolisForcing = 1;

const double length_mer = 40007.86 * 1000; // Earth 2*Pi*R in meridian direction
const double length_ekv = 40075.017 * 1000; // Earth 2*Pi*R in equator direction
const double pi = 3.141592653589793238462643; // Pi value
const double rad = pi / 180.0; // degrees to radian factor
const double Earth_R = 6371200;

const double SunAxis = 1.0; // Semi-major axis between Earth and Sun in astronomical unit
const double MoonAxis = 0.00257188153; // Semi-major axis between Earth and Moon in astronomical unit

const double outputMaxSizeMB = 800; // max size of one output file in MB

const int border_WALL[4][8] = {
	//BOTTOM, RIGHT, TOP, LEFT, RB_CORNER, LB_CORNER, RT_CORNER, LT_CORNER
	{1, 1, 1, 1, 1, 1, 1, 1 }, // HEIGHT
	{ 1, -1, 1, -1, -1, -1, -1, -1}, // VELOCITY_X
	{-1, 1, -1, 1, -1, -1, -1, -1 }, // VELOCITY_Y
	{ 1, 1, 1, 1, 1, 1, 1, 1 } // CONCENTRATION
};

const std::string Data_folder = "Data"; // main directory that contains all data folders

enum TypeOfPoint {
	BOTTOM = 0,
	RIGHT = 1,
	TOP = 2,
	LEFT = 3,
	RB_CORNER = 4,
	LB_CORNER = 5,
	RT_CORNER = 6,
	LT_CORNER = 7,
	INTERNAL = 8,	
	INTERNALWALL = 9,
	EXCLUDED = -9
};

enum TypeOfVariable {
	HEIGHT = 0,
	VELOCITY_X = 1,
	VELOCITY_Y = 2,
	CONCENTRATION = 3
};

enum BoundaryConditions {
	OPEN_BORDER = 1,
	WALL = -1,
	CONSTANT_VALUE = 0,
	FROM_FILE = 2
};

//char *name_of_boundary_condition[] = {"Open border", "Wall", "Fixed boundary" };
//const float normirovka = 1000.0000000000; //! что это? 