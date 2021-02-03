#include "Constants.h"

#include <iostream>
#include <fstream>
#include <algorithm>

//const std::string Data_folder = "Data";
//const std::string tideHarmonicsFolder = "tideHarmonics";

extern const char* HRM[] = { "K1",
						"K2",
						"M2",
						"N2",
						"O1",
						"P1",
						"Q1",
						"S2"
};

extern bool tideSide[4] = { false, false, false, false };

extern const double qTide[tideNum] = {rad*15.041069, 
										rad*30.082137, 
										rad*28.984104, 
										rad*28.439730, 
										rad*13.943036, 
										rad*13.398661,
										rad*14.958931,
										rad * 30};
// === DEFAULT VALUES OF TECHNICAL PARAMETERS ===

// use OpenMP
bool parallelOpenMP = true;
// include forcing to regularization
int F_reg = 0;
int Phi_reg = 0;

// mass fluxes correction for dry zone condition 
bool massFluxCorrection = true; // flag
double epsFlux = 1e-4;

bool ignore_warning = false; // ignore warning of existing case folder

// type of transport equation regularization: alpha_c = 1 - normal, alpha_c = 0 - simplified, also could be between (0,1).
double alpha_c = 1.0;
// Coefficient of viscosity in the transport equation, basic = 0.0, for special cases = 1.0/gc
double NSC = 0.1;


void initConfiguration(std::string configFileName)
{
	if (configFileName != "")
	{
		// std::ifstream is RAII, i.e. no need to call close
		std::ifstream cFile(configFileName);
		if (cFile.is_open())
		{
			std::cout << "Reading configuration file \"" << configFileName << "\":" << std::endl;
			std::string line;
			while (getline(cFile, line)) {
				line.erase(std::remove_if(line.begin(), line.end(), isspace),
					line.end());
				if (line[0] == '#' || line.empty())
					continue;
				auto delimiterPos = line.find("=");
				auto name = line.substr(0, delimiterPos);
				auto value = line.substr(delimiterPos + 1);
				std::cout << name << " = " << std::stod(value) << std::endl;
				double dvalue = std::stod(value);
				if (name == "F_reg")
				{
					if (checkRange(name, dvalue, 0, 1))
						F_reg = std::stoi(value);
				}
				if (name == "Phi_reg")
				{
					if (checkRange(name, dvalue, 0, 1))
						Phi_reg = std::stoi(value);
				}
				if (name == "parallelOpenMP")
				{
					if (checkRange(name, dvalue, 0, 1))
						parallelOpenMP = std::stoi(value);
				}
				if (name == "massFluxCorrection")
				{
					if (checkRange(name, dvalue, 0, 1))
						massFluxCorrection = std::stoi(value);
				}
				if (name == "epsFlux")
				{
					epsFlux = std::stod(value);
				}
				if (name == "ignore_warning")
				{
					if (checkRange(name, dvalue, 0, 1))
						ignore_warning = std::stoi(value);
				}
				if (name == "alpha_c")
				{
					if (checkRange(name, dvalue, 0, 1))
						alpha_c = std::stod(value);
				}
				if (name == "NSC")
				{
					if (checkRange(name, dvalue, 0, 1))
						NSC = std::stod(value);
				}
			}

		}
		else {
			std::cerr << "Couldn't open config file \"" << configFileName << "\" for reading.\n";
		}
	}	
}

bool checkRange(std::string name, double value, double minVal, double maxVal)
{
	bool result = value >= minVal && value <= maxVal;
	if (!result)
	{
		std::cerr << "WARNING! Variable " << name << " is out of range: [" << minVal << ", " << maxVal << "]." << std::endl;
	}
	return result;
}
