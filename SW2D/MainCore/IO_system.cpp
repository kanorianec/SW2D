#include "Raschet.h"
#include "technical.h"
//#include <iostream>
#include <fstream>


// Save all variable at time moment "Time_of_work"
void Raschet::Save_Grid() {
	string save_path = path + "/Grid"; // full path name to the preparing folder

	Prepare_Folder(save_path);

	string name_X = save_path + "/X.dat";
	string name_Y = save_path + "/Y.dat";
	string name_B = save_path + "/B.dat";


	std::ofstream fX(name_X, std::ios::out | std::ios::binary);
	std::ofstream fY(name_Y, std::ios::out | std::ios::binary);
	std::ofstream fB(name_B, std::ios::out | std::ios::binary);

	fX.write(reinterpret_cast<const char*> (X), sizeof(double) * Nx * Ny);
	fY.write(reinterpret_cast<const char*> (Y), sizeof(double) * Nx * Ny);
	fB.write(reinterpret_cast<const char*> (B), sizeof(double) * Nx * Ny);

	fX.close();
	fY.close();
	fB.close();
}

// Save all variable at time moment "Time_of_work"
void Raschet::Save_Data() {
	string save_path = path + "/" + to_str(Time_elapsed); // full path name to the preparing folder

	Prepare_Folder(save_path);

	string name_h = save_path + "/H.dat";
	string name_xU = save_path + "/xU.dat";
	string name_yU = save_path + "/yU.dat";
	string name_C = save_path + "/C.dat";

	std::ofstream fH (name_h, std::ios::out | std::ios::binary);
	std::ofstream fxU(name_xU, std::ios::out | std::ios::binary);
	std::ofstream fyU(name_yU, std::ios::out | std::ios::binary);
	std::ofstream fC;
	
	fH.write(reinterpret_cast<const char*> (H), sizeof(double) * Nx * Ny);
	fxU.write(reinterpret_cast<const char*> (xU), sizeof(double) * Nx * Ny);
	fyU.write(reinterpret_cast<const char*> (yU), sizeof(double) * Nx * Ny);

	fH.close();
	fxU.close();
	fyU.close();
	
	if (TransportProblemFlag)
	{
		fC.open(name_C, std::ios::out | std::ios::binary);
		fC.write(reinterpret_cast<const char*> (C), sizeof(double) * Nx * Ny);
		fC.close();
	}
}



// Reastart: reading all variables from time moment "Time_moment"
void Raschet::Restart_from_time_moment(double Time_moment) {

	string load_path = path + "/Data/" + to_str(Time_moment);

	cout << "Restarting from time moment " << Time_moment << endl;

	T_begin = Time_moment;
	restart = true;

	if (folderExists(load_path))
	{
		string name_h = load_path + "/H.dat";
		string name_xU = load_path + "/xU.dat";
		string name_yU = load_path + "/yU.dat";
		string name_C = load_path + "/C.dat";

		std::ifstream fH(name_h, std::ios::binary);
		std::ifstream fxU(name_xU, std::ios::binary);
		std::ifstream fyU(name_yU, std::ios::binary);
		std::ifstream fC;

		fH.read(reinterpret_cast<char*> (H), sizeof(double) * Nx * Ny);
		fxU.read(reinterpret_cast<char*> (H), sizeof(double) * Nx * Ny);
		fyU.read(reinterpret_cast<char*> (H), sizeof(double) * Nx * Ny);

		fH.close();
		fxU.close();
		fyU.close();

		if (TransportProblemFlag)
		{
			fC.open(name_C, std::ios::out | std::ios::binary);
			fC.read(reinterpret_cast<char*> (H), sizeof(double) * Nx * Ny);
			fC.close();
		}
	}
	else
	{
		cout << "Error: folder " << load_path << " does not exist!" << endl;
		cin.ignore(1024, '\n');
		cout << "Press Enter to stop the programm." << endl;
		cin.get();
		exit(0);
	}
			
};