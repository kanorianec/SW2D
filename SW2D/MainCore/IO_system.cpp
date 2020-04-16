#include "Raschet.h"
#include "technical.h"
//#include <iostream>
#include <fstream>


// output input data procedure
void Raschet::outputInputs()
{
	if (!restart)
	{
		cout << "Visualization of initial data; " << GetTimeStamp();
		std::ofstream fLog(path + "/timeLog.dat", std::ios::out);
		fLog << "Visualization of initial data; " << GetTimeStamp();
		fLog.close();

		if (binaryOutputFlag)
		{
			Save_Grid();
			Save_Data();			
		}

		if (Visualization_to_techplot_flag) 
		{
			Visualization_to_techplot_input();
			Visualization_to_techplot_result();
		}
	}
	t_graph_export = t_graph_export + t_step;
}

// output results procedure
void Raschet::outputResults()
{
	std::ofstream fLog(path + "/timeLog.dat", std::ios::out | std::ios::app);
	fLog << "Visualization time moment = " << t_graph_export << "; time =  ";
	

	cout << "Visualization time moment = " << t_graph_export << "; time =  ";
	if (t_step >= 3600.0)
	{
		cout << Time_elapsed / 3600.0 << " hours; " << GetTimeStamp() << flush;
		fLog << Time_elapsed / 3600.0 << " hours; " << GetTimeStamp();
	}		
	else
	{
		cout << Time_elapsed << "; " << GetTimeStamp() << flush;
		fLog << Time_elapsed << "; " << GetTimeStamp();
	}
		
	fLog.close();

	if (binaryOutputFlag)
		Save_Data();

	if (Visualization_to_techplot_flag)
		Visualization_to_techplot_result();
}

// Save all variable at time moment "Time_of_work"
void Raschet::Save_Grid() {
	string save_path = path + "/Grid"; // full path name to the preparing folder

	Prepare_Folder(save_path);

	string name_X = save_path + "/X.dat";
	string name_Y = save_path + "/Y.dat";
	
	string name_Lon = save_path + "/Lon.dat";
	string name_Lat = save_path + "/Lat.dat";
	
	string name_B = save_path + "/B.dat";
	string name_t = save_path + "/t.dat";

	std::ofstream fX(name_X, std::ios::out | std::ios::binary);
	std::ofstream fY(name_Y, std::ios::out | std::ios::binary);
	
	std::ofstream fB(name_B, std::ios::out | std::ios::binary);

	std::ofstream fNN(save_path + "/Nx_Ny.dat", std::ios::out);
	fNN << Nx << " " << Ny;

	fX.write(reinterpret_cast<const char*> (X), sizeof(double) * Nx * Ny);
	fY.write(reinterpret_cast<const char*> (Y), sizeof(double) * Nx * Ny);
	
	if (hlat + hlon > 0)
	{
		std::ofstream fLon(name_Lon, std::ios::out | std::ios::binary);
		std::ofstream fLat(name_Lat, std::ios::out | std::ios::binary);

	    fLon.write(reinterpret_cast<const char*> (Lon), sizeof(double) * Nx * Ny);
	    fLat.write(reinterpret_cast<const char*> (Lat), sizeof(double) * Nx * Ny);

		fLon.close();
		fLat.close();
	}
	
	fB.write(reinterpret_cast<const char*> (B), sizeof(double) * Nx * Ny);

	fNN.close();
	fX.close();
	fY.close();	
	fB.close();
	
	std::ofstream ft(name_t, std::ios::out);
	ft.close();
}

// Save all variable at time moment "Time_of_work"
void Raschet::Save_Data() {
	string save_path = path + "/" + to_str(Time_elapsed); // full path name to the preparing folder

	Prepare_Folder(save_path);

	string name_h = save_path + "/H.dat";
	string name_xU = save_path + "/xU.dat";
	string name_yU = save_path + "/yU.dat";
	string name_C = save_path + "/C.dat";

	string name_Fx = save_path + "/ForceX.dat";
	string name_Fy = save_path + "/ForceY.dat";
	
	string name_t = path + "/Grid/t.dat";

	std::ofstream fH (name_h, std::ios::out | std::ios::binary);
	std::ofstream fxU(name_xU, std::ios::out | std::ios::binary);
	std::ofstream fyU(name_yU, std::ios::out | std::ios::binary);
	//std::ofstream fFx(name_Fx, std::ios::out | std::ios::binary);
	//std::ofstream fFy(name_Fy, std::ios::out | std::ios::binary);
	
	std::ofstream fC;
	
	std::ofstream ft(name_t, std::ios::out | std::ios::app);
	
	ft << to_str(Time_elapsed) << " ";
	
	fH.write(reinterpret_cast<const char*> (H), sizeof(double) * Nx * Ny);
	fxU.write(reinterpret_cast<const char*> (xU), sizeof(double) * Nx * Ny);
	fyU.write(reinterpret_cast<const char*> (yU), sizeof(double) * Nx * Ny);
	//fFx.write(reinterpret_cast<const char*> (xWind), sizeof(double) * Nx * Ny); //ForceX
	//fFy.write(reinterpret_cast<const char*> (yWind), sizeof(double) * Nx * Ny); //ForceY

	fH.close();
	fxU.close();
	fyU.close();
	//fFx.close();
	//fFy.close();
	ft.close();
	
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

	if (t_graph_export < T_begin)
		t_graph_export = T_begin;

	restart = true;

	if (!folderNotExists(load_path))
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
		fxU.read(reinterpret_cast<char*> (xU), sizeof(double) * Nx * Ny);
		fyU.read(reinterpret_cast<char*> (yU), sizeof(double) * Nx * Ny);

		fH.close();
		fxU.close();
		fyU.close();

		if (TransportProblemFlag)
		{
			fC.open(name_C, std::ios::out | std::ios::binary);
			fC.read(reinterpret_cast<char*> (C), sizeof(double) * Nx * Ny);
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