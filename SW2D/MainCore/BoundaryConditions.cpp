#include "Raschet.h"
#include "technical.h"
#include <cmath>
#include <algorithm>

// Set up boundary conditions of d/dn = 0
void Raschet::SetZeroDerivativeConditions(TypeOfVariable VType, TypeOfPoint PType1, TypeOfPoint PType2, TypeOfPoint PType3, TypeOfPoint PType4, TypeOfPoint PType5, TypeOfPoint PType6, TypeOfPoint PType7, TypeOfPoint PType8)
{
	const int VarNum = 8;
	TypeOfPoint tempType[VarNum] = { PType1, PType2, PType3, PType4, PType5, PType6, PType7, PType8 };
	for (int i = 0; i < VarNum; i++)
	{
		if (tempType[i] != EXCLUDED)
		{
			if (VType == CONCENTRATION && !TransportProblemFlag)
			{
				cout << "ERROR setting zero derivative conditions for concentration: there is NO TRANSPORT EQUATION in the PROBLEM DEFINITION!" << endl;
				pause();
			}
			border[VType][tempType[i]] = 1;
			border_C[VType][tempType[i]] = 0.0;
		}
	}
}

// Fixed "Value" for variable "Vtype" in the "Ptype" border
void Raschet::SetFixedBoundaryConditions(TypeOfVariable VType, TypeOfPoint PType, double Value)
{
	if (VType == CONCENTRATION && !TransportProblemFlag)
	{
		cout << "ERROR setting fixed boundary conditions for concentration: there is NO TRANSPORT EQUATION in the PROBLEM DEFINITION!" << endl;
		pause();
	}

	if (PType != EXCLUDED)
	{
		border[VType][PType] = -1;
		border_C[VType][PType] = Value;

		switch (PType) // set up adjacent corners
		{
		case LEFT:
			border[VType][LT_CORNER] = -1;
			border_C[VType][LT_CORNER] = Value;

			border[VType][LB_CORNER] = -1;
			border_C[VType][LB_CORNER] = Value;
			break;
		case RIGHT:
			border[VType][RT_CORNER] = -1;
			border_C[VType][RT_CORNER] = Value;

			border[VType][RB_CORNER] = -1;
			border_C[VType][RB_CORNER] = Value;
			break;
		case TOP:
			border[VType][LT_CORNER] = -1;
			border_C[VType][LT_CORNER] = Value;

			border[VType][RT_CORNER] = -1;
			border_C[VType][RT_CORNER] = Value;
			break;
		case BOTTOM:
			border[VType][LB_CORNER] = -1;
			border_C[VType][LB_CORNER] = Value;

			border[VType][RB_CORNER] = -1;
			border_C[VType][RB_CORNER] = Value;
			break;
		}
	}
}

// Open boundary conditions: du_n/dn = du_tau/dn = dh/dn = 0;
void Raschet::SetOpenBoundaryConditions(TypeOfPoint PType1, TypeOfPoint PType2, TypeOfPoint PType3, TypeOfPoint PType4, TypeOfPoint PType5, TypeOfPoint PType6, TypeOfPoint PType7, TypeOfPoint PType8)
{
	SetZeroDerivativeConditions(HEIGHT, PType1, PType2, PType3, PType4, PType5, PType6, PType7, PType8);
	SetZeroDerivativeConditions(VELOCITY_X, PType1, PType2, PType3, PType4, PType5, PType6, PType7, PType8);
	SetZeroDerivativeConditions(VELOCITY_Y, PType1, PType2, PType3, PType4, PType5, PType6, PType7, PType8);
	if (TransportProblemFlag)
	{
		SetZeroDerivativeConditions(CONCENTRATION, PType1, PType2, PType3, PType4, PType5, PType6, PType7, PType8);
	}

	t2_bound = -1; // COSTIL!
}

// Wall boundary conditions: u_n = 0, du_tau/dn = 0, dh/dn = 0;
void Raschet::SetWallBoundaryConditions(TypeOfPoint PType1, TypeOfPoint PType2, TypeOfPoint PType3, TypeOfPoint PType4)
{
	const int VarNum = 4;
	TypeOfPoint tempType[VarNum] = { PType1, PType2, PType3, PType4 };
	for (int i = 0; i < VarNum; i++)
	{
		if (tempType[i] != EXCLUDED)
		{
			SetZeroDerivativeConditions(HEIGHT, tempType[i]);
			if (TransportProblemFlag)
			{
				SetZeroDerivativeConditions(CONCENTRATION, tempType[i]);
			}

			if (tempType[i] == LEFT || tempType[i] == RIGHT)
			{
				SetZeroDerivativeConditions(VELOCITY_Y, tempType[i]);
				SetFixedBoundaryConditions(VELOCITY_X, tempType[i], 0.0);
			}
			else
			{
				SetZeroDerivativeConditions(VELOCITY_X, tempType[i]);
				SetFixedBoundaryConditions(VELOCITY_Y, tempType[i], 0.0);
			}
		}
	}
}



void Raschet::SetFileBoundaryConditions(TypeOfVariable VType, TypeOfPoint PType1, TypeOfPoint PType2, TypeOfPoint PType3, TypeOfPoint PType4)
{
	BoundaryConditionsFromFile = true;
	if (t2_bound < 0)
	{
		FT = fopen("borders/time_shift.dat", "r");
		if (FT == NULL) 
		{
			std::cout << "Error: borders/time_shift.dat was not found!" << endl;
			pause();
		}			
		fscanf(FT, "%lf\n", &t2_bound);		
	}

	const int VarNum = 4;
	TypeOfPoint tempType[VarNum] = { PType1, PType2, PType3, PType4 };
	for (int i = 0; i < VarNum; i++)
	{
		TypeOfPoint PType = tempType[i];
		if (PType != EXCLUDED)
		{
			string Name = "borders/" + GetVariableName(VType) + "_" + to_string(PType) + ".dat";

			FV[VType][PType] = fopen(Name.c_str(), "r");
			if (FV[VType][PType] == NULL)
				std::cout << "Error: " << Name << " was not found!" << endl;

			int N = Nx;
			if (PType != EXCLUDED)
			{
				if (PType == LEFT || PType == RIGHT)
					N = Ny;

				lin_b[VType][PType] = new double[N];
				lin_k[VType][PType] = new double[N];

				border[VType][PType] = FROM_FILE;
				switch (PType) // set up corners
				{
				case LEFT:
					border[VType][LT_CORNER] = FROM_FILE;
					border[VType][LB_CORNER] = FROM_FILE;
					break;
				case RIGHT:
					border[VType][RT_CORNER] = FROM_FILE;
					border[VType][RB_CORNER] = FROM_FILE;
					break;
				case TOP:
					border[VType][LT_CORNER] = FROM_FILE;
					border[VType][RT_CORNER] = FROM_FILE;
					break;
				case BOTTOM:
					border[VType][LB_CORNER] = FROM_FILE;
					border[VType][RB_CORNER] = FROM_FILE;
					break;
				}
			}

			for (int k = 0; k < N; k++)
			{
				double y1;
				fscanf(FV[VType][PType], "%lf ", &y1);
				lin_b[VType][PType][k] = y1;
				lin_k[VType][PType][k] = 0.0;
			}
		}
	}
}

void Raschet::RecalcFileBoundaryConditions()
{
	TypeOfVariable FinalVtype;
	if (TransportProblemFlag)
		FinalVtype = CONCENTRATION;
	else
		FinalVtype = VELOCITY_Y;

	if (t2_bound < Time_elapsed + dT)
	{
		t1_bound = t2_bound;
		fscanf(FT, "%lf\n", &t2_bound);

		for (int PType = BOTTOM; PType <= LEFT; PType++)
			for (int VType = HEIGHT; VType <= FinalVtype; VType++)
				if (border[VType][PType] == FROM_FILE)
				{
					int N = Nx;
					if (PType == LEFT || PType == RIGHT)
						N = Ny;
					for (int k = 0; k < N; k++)
					{
						double y1;
						double y2;
						y1 = lin_k[VType][PType][k] * t1_bound + lin_b[VType][PType][k];
						fscanf(FV[VType][PType], "%lf ", &y2);
						lin_k[VType][PType][k] = (y2 - y1) / (t2_bound - t1_bound);
						lin_b[VType][PType][k] = (y1 * t2_bound - y2 * t1_bound) / (t2_bound - t1_bound);
					}
				}
	}


	for (int PType = BOTTOM; PType <= LEFT; PType++)
	{
		int Na, Nb, N, dn;
		switch (PType) // set up corners
		{
		case LEFT:
			Na = 1;
			Nb = 0;
			N = Ny;
			dn = Ny;
			break;
		case RIGHT:
			Na = 1;
			Nb = (Nx - 1)*Ny;
			N = Ny;
			dn = -Ny;
			break;
		case TOP:
			Na = Ny;
			Nb = Ny - 1;
			N = Nx;
			dn = -1;
			break;
		case BOTTOM:
			Na = Ny;
			Nb = 0;
			N = Nx;
			dn = 1;
			break;
		}

		for (int VType = HEIGHT; VType <= FinalVtype; VType++)
			if (border[VType][PType] == FROM_FILE)
			{
				switch (VType) // set up corners
				{
				case HEIGHT:
				#pragma omp parallel for
					for (int k = 0; k < N; k++)
					{
						Ht[Na*k + Nb] = (Time_elapsed + dT) * lin_k[VType][PType][k] + lin_b[VType][PType][k];
						Ht[Na*k + Nb] -= Bmin + B[Na*k + Nb];
					}
					break;
				case VELOCITY_X:
				#pragma omp parallel for
					for (int k = 0; k < N; k++)
					{
						xUt[Na*k + Nb] = (Time_elapsed + dT) * lin_k[VType][PType][k] + lin_b[VType][PType][k];
					}
					break;
				case VELOCITY_Y:
				#pragma omp parallel for
					for (int k = 0; k < N; k++)
					{
						yUt[Na*k + Nb] = (Time_elapsed + dT) * lin_k[VType][PType][k] + lin_b[VType][PType][k];
					}
					break;
				case CONCENTRATION:
				#pragma omp parallel for
					for (int k = 0; k < N; k++)
					{
						Ct[Na*k + Nb] = (Time_elapsed + dT) * lin_k[VType][PType][k] + lin_b[VType][PType][k];
					}
					break;
				}
			}
		
		if (COMBINE_FILE_AND_FREE_BOUNDARY_CONDITIONS)
		{
			#pragma omp parallel for
			for (int k = 0; k < N; k++)
			{
				int i = Na*k + Nb;
				if (PType == LEFT || PType == RIGHT)
				{
					if (dn*xUt[i] < 0)
					{
						xUt[i] = xUt[i + dn];
						yUt[i] = yUt[i + dn];
						//Ht[i] = Ht[i + dn];
					}
				}
				else
				{
					if (dn*yUt[i] < 0)
					{
						xUt[i] = xUt[i + dn];
						yUt[i] = yUt[i + dn];
						//Ht[i] = Ht[i + dn];
					}
				}
			}
		}
	}

	/*
	for (int k = 0; k < Nx; k++)
	{
	cout << yUt[k*Ny + 0] << " ";
	}
	pause();	*/
}

void Raschet::SetInternalWall(TypeOfPoint PType, int solid_ind, int start_ind, int end_ind)
{
	InternalWallsFlag = true;
	int Na = Ny;
	int Nb = solid_ind;
	if (PType == LEFT || PType == RIGHT)
	{
		Na = 1;
		Nb = solid_ind*Ny;
	}
	for (int k = start_ind; k <= end_ind; k++)
		S[Na*k + Nb] = PType + INTERNALWALL;
	int corner_1;
	int corner_2;
	switch (PType) // set up adjacent corners
	{
	case LEFT:
		corner_1 = LB_CORNER + INTERNALWALL;
		corner_2 = LT_CORNER + INTERNALWALL;
		break;
	case RIGHT:
		corner_1 = RB_CORNER + INTERNALWALL;
		corner_2 = RT_CORNER + INTERNALWALL;
		break;
	case TOP:
		corner_1 = RT_CORNER + INTERNALWALL;
		corner_2 = LT_CORNER + INTERNALWALL;
		break;
	case BOTTOM:
		corner_1 = RB_CORNER + INTERNALWALL;
		corner_2 = LB_CORNER + INTERNALWALL;
		break;
	}
	if (start_ind != 0 && start_ind != Nx - 1 && start_ind != Ny - 1)
		S[Na*start_ind + Nb] = corner_1;
	if (end_ind != 0 && end_ind != Nx - 1 && end_ind != Ny - 1)
		S[Na*end_ind + Nb] = corner_2;
}

void Raschet::SetTidesHarmonicsBoundaryConditions(TypeOfPoint PType1, TypeOfPoint PType2, TypeOfPoint PType3, TypeOfPoint PType4)
{
	tidesHarmonics = true;
	const int VarNum = 4;
	TypeOfPoint tempType[VarNum] = { PType1, PType2, PType3, PType4};

	if (!folderNotExists(tideHarmonicsFolder))
	{
		for (int i = 0; i < VarNum; i++)
		{
			TypeOfPoint PType = tempType[i];
			if (PType != EXCLUDED)
			{
				if (border[HEIGHT][PType] != FROM_FILE && border[HEIGHT][PType] != CONSTANT_VALUE)
					border[HEIGHT][PType] = ONLY_TIDE;
				tideSide[PType] = true;
				int N = Nx;
				if (PType == LEFT || PType == RIGHT)
					N = Ny;

				for (int h = 0; h < tideNum; h++)
				{
					tideA[PType][h] = new double[N];
					tidePh[PType][h] = new double[N];

					string fname = tideHarmonicsFolder + "/" + HRM[h] + "_" + to_str(PType) + "_";
					std::ifstream fAmp(fname + "AMP.dat", std::ios::out | std::ios::binary);
					if (!fAmp.is_open())
						cout << "Can't open " + fname + "AMP.dat" << endl;
					std::ifstream fPh(fname + "PH.dat", std::ios::out | std::ios::binary);
					if (!fPh.is_open())
						cout << "Can't open " + fname + "PH.dat" << endl;

					fAmp.read(reinterpret_cast<char*> (tideA[PType][h]), sizeof(double) *N);
					fPh.read(reinterpret_cast<char*> (tidePh[PType][h]), sizeof(double) *N);
					/*
					for (int o = 0; o < N; o++)
					{
						cout << HRM[h] << " " << PType << " " << tideA[PType][h][o] << endl;
					}
					pause();
					*/

					fAmp.close();
					fPh.close();
				}
			}
		}
	}
	else
	{
		cout << "Error: folder " << tideHarmonicsFolder << " does not exist!" << endl;
		cin.ignore(1024, '\n');
		cout << "Press Enter to stop the programm." << endl;
		cin.get();
		exit(0);
	}
}

void Raschet::addTidesHarmonicsBoundaryConditions()
{
	// convert time to convenient format
	struct tm* tt;
	time_t to_convert = dT + Time_elapsed + RaschetTime;

	tt = gmtime(&to_convert);

	int D = tt->tm_yday;
	int Y = tt->tm_year + 1900;
	int hour = tt->tm_hour;
	double minute = tt->tm_min;
	double second = tt->tm_sec;
	double t = hour + minute / 60 + second / 3600;

	//int month = tt->tm_mon + 1;	
	//int minute = tt->tm_min;
	//int second = tt->tm_sec;

	int L = (int)(0.25*(Y - 1901));
	int N = (259.157 - 19.3282*(Y - 1900) - 0.05295*(D + L))*rad;
	JulianDay jd = Calendar::julianDayForDateAndTime(tt->tm_mday, tt->tm_mon + 1, tt->tm_year + 1900, hour, tt->tm_min, tt->tm_sec);
	double T = (jd.day + jd.time - 2415020.0) / 36525;
	double h = 279.6967 + T*36000.7689 + T*T*0.0003032;
	double s = 270.434164 + 481267.8831*T - 0.001133*T*T + 0.0000019*T*T*T;

	double nu = 12.94*sin(N) - 1.34*sin(2 * N) + 0.19*sin(3 * N);
	double ksi = 11.87*sin(N) - 1.34*sin(2 * N) + 0.19*sin(3 * N);
	double nut = 8.86*sin(N) - 0.68*sin(2 * N) + 0.07*sin(3 * N);
	double nutt = 0.5*(17.74*sin(N) - 0.68*sin(2 * N) + 0.04*sin(3 * N));
	double p = 334.385 + 40.6625*(Y - 1900) + 0.11140*(D + L);

	double f[tideNum];
	double u[tideNum];
	double v0[tideNum];
	double q[tideNum];

	// K1: q[K1] = 15.041069;
	f[K1] = 1.0060 + 0.1150*cos(N) - 0.0088*cos(2 * N) + 0.0006*cos(3 * N);
	u[K1] = -nut;
	v0[K1] = h + 90;
	
	// K2: q['K2'] = 30.082137;
	f[K2] = 1.0241 + 0.2863*cos(N) + 0.0083*cos(2 * N) - 0.0015*cos(3 * N);
	u[K2] = -2 * nutt;
	v0[K2] = 2 * h;

	// M2: q['M2'] = 28.984104;
	f[M2] = 1.0004 - 0.0373*cos(N) + 0.0002*cos(2 * N);
	u[M2] = 2 * ksi - 2 * nu;
	v0[M2] = 2 * h - 2 * s;

	// N2: q['N2'] = 28.439730;
	f[N2] = f[M2];
	u[N2] = u[M2];
	v0[N2] = 2 * h - 3 * s + p;

	// O1: q['O1'] = 13.943036;
	f[O1] = 1.0089 + 0.1871*cos(N) - 0.0147*cos(2 * N) + 0.0014*cos(3 * N);
	u[O1] = 2 * ksi - nu;
	v0[O1] = h - 2 * s - 90;

	// Q1: q[Q1] = 13.398661;
	f[Q1] = f[O1];
	u[Q1] = u[O1];
	v0[Q1] = h - 3 * s + p - 90;
	
	// P1: q['P1'] = 14.958931;
	f[P1] = 1.000;
	u[P1] = 0;
	v0[P1] = -h - 90;

	// S2: q[S2] = 30;
	f[S2] = 1.000;
	u[S2] = 0;
	v0[S2] = 0;
	
	for (int i = 0; i < tideNum; i++)
	{
		u[i] *= rad;
		v0[i] *= rad;
		q[i] *= rad;
	}
	

	for (int PType = BOTTOM; PType <= LEFT; PType++)
	{
		if (tideSide[PType])
		{
			if (border[HEIGHT][PType] != CONSTANT_VALUE && border[HEIGHT][PType] != FROM_FILE && border[HEIGHT][PType] != ONLY_TIDE)
			{
				cout << PType << " boundary conditions has to be const value or from file or ONLY_TIDE!" << endl;
				cin.ignore(1024, '\n');
				Stop_Raschet_Flag = 1;
			}
			else
			{
				int Na, Nb, NN;
				switch (PType) // set up corners
				{
				case LEFT:
					Na = 1;
					Nb = 0;
					NN = Ny;
					break;
				case RIGHT:
					Na = 1;
					Nb = (Nx - 1)*Ny;
					NN = Ny;
					break;
				case TOP:
					Na = Ny;
					Nb = Ny - 1;
					NN = Nx;
					break;
				case BOTTOM:
					Na = Ny;
					Nb = 0;
					NN = Nx;
					break;
				}
				#pragma omp parallel for
				for (int k = 0; k < NN; k++)
				{
					if (Ht[Na*k + Nb] > eps)
					{
						for (int it = 0; it < tideNum; it++)
						{
							if (border[HEIGHT][PType] == ONLY_TIDE)
								Ht[Na*k + Nb] = std::max(-B[Na*k + Nb] - Bmin, eps);
							Ht[Na*k + Nb] += f[it] * 0.01 * tideA[PType][it][k] * (cos(t*qTide[it] + v0[it] + u[it] - tidePh[PType][it][k]));
							if (Ht[Na*k + Nb] < 0)
								Ht[Na*k + Nb] = eps;
						}
					}
					
				}
			}
		}
	}
}
