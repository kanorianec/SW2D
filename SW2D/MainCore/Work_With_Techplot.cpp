﻿/*
Visualization function description file. The main way of visualization is using "Tecplot 360 EX 2016 R2"
*/
#define _CRT_SECURE_NO_WARNINGS
#include "Raschet.h"
#include "Constants.h"
#include "technical.h"

#include <iostream>
#include <algorithm>
#include <fstream>

//#include <stdio.h>
//#include <string>
//#include <math.h>
//#include <time.h>
//#include <fstream>

using namespace std;

/*
format of output tecplot file:
TITLE= Test_name
VARIABLES = "X" "Y" "H" "B" "Ux" "Uy" "S" // ordered list of variables
ZONE T="Time 0.000000" DATAPACKING=POINT, NODES=10000, ELEMENTS=9801, ZONETYPE=FEQUADRILATERAL // Data descriptor

Next, the data is displayed in columns 

After that, the indices of the points are displayed in a specific form (as it necessary for the Tecplot)
*/

void Raschet::SetVisualizationProperties(double Vis_T_start, double Vis_T_end, int iMin, int jMin, int iMax, int jMax)
{
	Raschet::t_graph_export = max(Vis_T_start,T_begin);
	Raschet::t_end_graph_export = min(Vis_T_end,T_end);
	Raschet::i_maxVis = iMax;
	Raschet::i_minVis = iMin;
	Raschet::j_maxVis = jMax;
	Raschet::j_minVis = jMin;	
	Raschet::first_visualization = true;
	Raschet::output_per_file_counter = 0;
	Raschet::current_file_sizeMB = 0.0;

}

void Raschet::Visualization_to_techplot_input()
{
	string InputDataFileName = path + "/Input.dat";

	FILE *F_in = fopen(InputDataFileName.c_str(), "w");

	fprintf(F_in, "TITLE= \"%s\"\n", Test_name.c_str());
	//fprintf(F, "VARIABLES = \"X\" \"Y\" \"H\" \"B\" \"Ux\" \"Uy\" \"Phix\" \"Phiy\" \"Fx\" \"Fy\"\n");
	fprintf(F_in, "VARIABLES = \"X\" \"Y\" \"Longitude\" \"Latitude\" \"B\" \"H\" \"Ux\" \"Uy\" \"S\" ");
	if (TransportProblemFlag)
		fprintf(F_in, "\"C\" ");
	fprintf(F_in, "\"epsilon\" \n");
	fprintf(F_in, "ZONE T=\"Time %lf\" DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FEQUADRILATERAL\n", t_graph_export, (i_maxVis - i_minVis + 1)*(j_maxVis - j_minVis + 1), (i_maxVis - i_minVis)*(j_maxVis - j_minVis));
	for (int i = i_minVis; i <= i_maxVis; i++)
	{
		for (int j = j_minVis; j <= j_maxVis; j++)
		{
			int k = i*Ny + j;
			fprintf(F_in, "%lf %lf %lf %lf %lf %lf %lf %lf %d ", X[k], Y[k], Lon[k], Lat[k], B[k] + Bmin, H[k] + B[k] + Bmin, xU[k], yU[k], S[k]);
			if (TransportProblemFlag)
				fprintf(F_in, "%lf ",C[k]);
			fprintf(F_in, "%d\n", epsilon[k]);
		}
	}

	for (int i = 0; i < (i_maxVis - i_minVis); i++)
	{
		for (int j = 0; j < (j_maxVis - j_minVis); j++)
		{
			fprintf(F_in, "%d %d %d %d\n", ((j_maxVis - j_minVis + 1)*i + j + 1), ((j_maxVis - j_minVis + 1)*(i + 1) + j + 1), ((j_maxVis - j_minVis + 1)*(i + 1) + j + 2), ((j_maxVis - j_minVis + 1)*i + j + 2));
		}
	}

	fclose(F_in);
}

void Raschet::Visualization_to_techplot_result()
{
	double timeStr;
	if (T_begin > 3600.0)
		timeStr = T_begin / 3600.0;

	if (current_file_sizeMB > outputMaxSizeMB)
	{
		output_per_file_counter++;
		first_visualization = true;
		current_file_sizeMB = 0.0;
	}
	if (first_visualization) // a condition for outputting first part of output data, it is necessary to write a header once
	{
		string RaschetFileName = path + "/Result_t" + to_str(timeStr, 2) + "(" + to_string(output_per_file_counter) + ").dat"; //to_string(output_per_file_counter)
		FILE *F = fopen(RaschetFileName.c_str(), "w");
		
		//cout << "Visualization to techplot time moment = " << t_graph_export << "; time =  "<<Time_elapsed << "; " << RaschetFileName << "; " << GetTimeStamp();

		fprintf(F, "TITLE= \"%s\"\n", Test_name.c_str());
		//fprintf(F, "VARIABLES = \"X\" \"Y\" \"H\" \"B\" \"Ux\" \"Uy\" \"Phix\" \"Phiy\" \"Fx\" \"Fy\"\n");
		fprintf(F, "VARIABLES = \"X\" \"Y\" \"B\" \"Longitude\" \"Latitude\" \"H\" \"Ux\" \"Uy\" ");
		if (TransportProblemFlag)
			fprintf(F, "\"C\" ");
		fprintf(F, "\"Fx\" \"Fy\" \"Phix\" \"Phiy\"\n"); // 
		fprintf(F, "ZONE T=\"Time %lf\" DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FEQUADRILATERAL\n", t_graph_export, (i_maxVis - i_minVis + 1)*(j_maxVis - j_minVis + 1), (i_maxVis - i_minVis)*(j_maxVis - j_minVis));
		for (int i = i_minVis; i <= i_maxVis; i++)
		{
			for (int j = j_minVis; j <= j_maxVis; j++)
			{
				int k = i*Ny + j;		
				fprintf(F, "%lf %lf %lf %lf %lf %lf %lf %lf ", X[k], Y[k], B[k] + Bmin, Lon[k], Lat[k], min(H[k] + B[k] + Bmin, 200.0), xU[k], yU[k]);
				if (TransportProblemFlag)
					fprintf(F, "%lf ", C[k]);
				fprintf(F, "%lf %lf %lf %lf\n", H[k] * ForceX[k], H[k] * ForceY[k], PhiX[k], PhiY[k]); //, H[k] * ForceX[k], H[k] * ForceY[k]
			}
		}

		for (int i = 0; i < (i_maxVis - i_minVis); i++)
		{
			for (int j = 0; j < (j_maxVis - j_minVis); j++)
			{
				fprintf(F, "%d %d %d %d\n", ((j_maxVis - j_minVis + 1)*i + j + 1), ((j_maxVis - j_minVis + 1)*(i + 1) + j + 1), ((j_maxVis - j_minVis + 1)*(i + 1) + j + 2), ((j_maxVis - j_minVis + 1)*i + j + 2));
			}
		}

		fclose(F);
		first_visualization = false;
	}
	else
	{
		string RaschetFileName = path + "/Result_t" + to_str(timeStr, 2) + "(" + to_string(output_per_file_counter) + ").dat"; //to_string(output_per_file_counter)

		FILE *F = fopen(RaschetFileName.c_str(), "ab");

		//fprintf(F, "TITLE= \"%s\"\n", Test_name.c_str());
		//fprintf(F, "VARIABLES = \"X\" \"Y\" \"H\" \"B\" \"Ux\" \"Uy\" \"Phix\" \"Phiy\" \"Fx\" \"Fy\"\n");
		fprintf(F, "VARIABLES = \"X\" \"Y\" \"B\" \"Longitude\" \"Latitude\" \"H\" \"Ux\" \"Uy\" ");
		if (TransportProblemFlag)
			fprintf(F, "\"C\" ");
		fprintf(F, "\"Fx\" \"Fy\" \"Phix\" \"Phiy\"\n"); // 
		fprintf(F, "ZONE T=\"Time %lf\" DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FEQUADRILATERAL, VARSHARELIST = ([1, 2, 3, 4, 5] = 1), CONNECTIVITYSHAREZONE = 1\n", t_graph_export, (i_maxVis - i_minVis + 1)*(j_maxVis - j_minVis + 1), (i_maxVis - i_minVis)*(j_maxVis - j_minVis));
		for (int i = i_minVis; i <= i_maxVis; i++)
		{
			for (int j = j_minVis; j <= j_maxVis; j++)
			{
				int k = i*Ny + j;
				fprintf(F, "%lf %lf %lf ", min(H[k] + B[k] + Bmin, 200.0), xU[k], yU[k]);
				if (TransportProblemFlag)
					fprintf(F, "%lf ", C[k]);
				fprintf(F, "%lf %lf %lf %lf\n",H[k] * ForceX[k], H[k] * ForceY[k], PhiX[k], PhiY[k]); 

			}
		}

		fclose(F);
	}
	current_file_sizeMB += Nx*Ny * 9 * 7 / (1024 * 1024);
	//cout << current_file_sizeMB << endl;
	//pause();
}