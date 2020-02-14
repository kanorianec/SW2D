#include "Raschet.h"

void Raschet::Initialize_Transport_Problem(double InitialC[],
	double D)
{
	TransportProblemFlag = true;

	Raschet::C = new double[Nx*Ny]();
	Raschet::Ct = new double[Nx*Ny]();

	Raschet::D = D;

	for (int i = 0; i<Nx; i++)
	{
		for (int j = 0; j<Ny; j++)
		{
			Raschet::C[i*Ny + j] = InitialC[i*Ny + j];
		}
	}
	//reset open boundary conditions because of changing value of TransportProblemFlag
	SetOpenBoundaryConditions(TOP, BOTTOM, RIGHT, LEFT, RT_CORNER, LT_CORNER, RB_CORNER, LB_CORNER);
}