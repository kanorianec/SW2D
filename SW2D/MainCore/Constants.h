#ifndef CONSTANTS_H
#define CONSTANTS_H
#pragma once
#include <string>

// =========================================================
// Physical and Mathematical constants

const double gc = 10.0; // 9.81;
//const double G = 0.0000000000667408; // gravitational constant

const double SD = 1.2035; // Dudson's constant for Sun
const double MD = 2.6206; // Dudson's constant for Moon
						  //const int CoriolisForcing = 1;

const double length_mer = 40007.86 * 1000; // Earth 2*Pi*R in meridian direction
const double length_ekv = 40075.017 * 1000; // Earth 2*Pi*R in equator direction
const double pi = 3.141592653589793238462643; // Pi value
const double rad = pi / 180.0; // degrees to radian factor
const double Earth_R = 6371200;

const double n = 0.018;

const double SunAxis = 1.0; // Semi-major axis between Earth and Sun in astronomical unit
const double MoonAxis = 0.00257188153; // Semi-major axis between Earth and Moon in astronomical unit

//================= Technical Constants and parameters ====================
// use OpenMP
extern bool parallelOpenMP;
// include forcing to regularization
extern int F_reg;
extern int Phi_reg;

// mass fluxes correction for dry zone condition 
extern bool massFluxCorrection; // flag
extern double epsFlux;

const bool wellBalancedScheme = true;

// ignore warning of existing case folder
extern bool ignore_warning; 

const bool COMBINE_FILE_AND_FREE_BOUNDARY_CONDITIONS = false; //true;

// include forcing in boundary conditions
const int F_bound = 0;

// type of transport equation regularization: alpha_c = 1 - normal, alpha_c = 0 - simplified, also could be between (0,1).
extern double alpha_c;
// Coefficient of viscosity in the transport equation, basic = 0.0, for special cases = 1.0/gc
extern double NSC;


const bool binaryOutputFlag = true; // output to binary flag
const bool TXToutputFlag = true; // output to text flag

const double CriticalVal = 10000.0;

const int TideForcing = 0; //

const double outputMaxSizeMB = 800; // max size of one output file in MB

const int border_WALL[4][8] = {
	//BOTTOM, RIGHT, TOP, LEFT, RB_CORNER, LB_CORNER, RT_CORNER, LT_CORNER
	{ 1, 1, 1, 1, 1, 1, 1, 1 }, // HEIGHT
	{ 1, -1, 1, -1, -1, -1, -1, -1 }, // VELOCITY_X
	{ -1, 1, -1, 1, -1, -1, -1, -1 }, // VELOCITY_Y
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
	CONCENTRATION = 3,
	FLOW = 4,
	TIDE = 777
};

enum BoundaryConditions {
	OPEN_BORDER = 1,
	WALL = -1,
	CONSTANT_VALUE = 0,
	FROM_FILE = 2,
	ONLY_TIDE = 3
};

const int tideNum = 8; // number of tides harmonics
extern const char* HRM[];
const std::string tideHarmonicsFolder = "tidesHarmonics"; // directory of tide harmonics characteristics
extern bool tideSide[4];

enum tideType {  // tide harmonics
	K1,
	K2,
	M2,
	N2,
	O1,
	P1,
	Q1,
	S2
};

extern const double qTide[tideNum];

//char *name_of_boundary_condition[] = {"Open border", "Wall", "Fixed boundary" };
//const float normirovka = 1000.0000000000; // 

enum forceType {
	NO_FORCE = 0,
	REAL_FORCE = 1,
	CONST_FORCE = 2,	
};

void initConfiguration(std::string configFileName = ""); 
bool checkRange(std::string name, double value, double minVal, double maxVal);

// costil indication
// COSTIL! COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  

#endif