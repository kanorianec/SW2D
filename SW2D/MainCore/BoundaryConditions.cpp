#include "Raschet.h"
#include "technical.h"

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
				system("pause");
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
		system("pause");
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
		FT = fopen("borders\\time_shift.dat", "r");
		fscanf(FT, "%lf\n", &t2_bound);
		if (FT == NULL)
			std::cout << "Error: borders\\time_shift.dat was not found!" << endl;
	}

	const int VarNum = 4;
	TypeOfPoint tempType[VarNum] = { PType1, PType2, PType3, PType4 };
	for (int i = 0; i < VarNum; i++)
	{
		TypeOfPoint PType = tempType[i];
		if (PType != EXCLUDED)
		{
			string Name = "borders\\" + GetVariableName(VType) + "_" + to_string(PType) + ".dat";

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
		int Na, Nb, N;
		switch (PType) // set up corners
		{
		case LEFT:
			Na = 1;
			Nb = 0;
			N = Ny;
			break;
		case RIGHT:
			Na = 1;
			Nb = (Nx - 1)*Ny;
			N = Ny;
			break;
		case TOP:
			Na = Ny;
			Nb = Ny - 1;
			N = Nx;
			break;
		case BOTTOM:
			Na = Ny;
			Nb = 0;
			N = Nx;
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
	}

	/*
	for (int k = 0; k < Nx; k++)
	{
	cout << yUt[k*Ny + 0] << " ";
	}
	system("pause");	*/
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
