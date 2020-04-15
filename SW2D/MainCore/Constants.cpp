#include "Constants.h"

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

