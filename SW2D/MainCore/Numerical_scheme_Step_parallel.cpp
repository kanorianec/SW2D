/*
File of calculation module, includes difference scheme
*/
#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <stdio.h>
#include <iostream>
//#include <fstream>

#include "Raschet.h"
#include "Constants.h"
#include <omp.h>
#include "technical.h"
#include <algorithm>

using namespace std;

// parallel realization (using openMP) 
void Raschet::Numerical_scheme_time_step_parallel()
{
	//omp_set_num_threads(4);

	Recalc_forces_parallel(); // forces calculating 

	// mass fluxes xJm, yJm calculation 
	//#pragma omp parallel for
	for (int i = 0; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny - 1; j++)
		{
			double xW_05R, yW_05T;
			//double xJ_05R, yJ_05T;
			double B_c, B_1R, B_1T;
			double H_05R, H_05T, H_c, H_05L05T, H_05R05B, H_05R05T, H_1R, H_1T;
			double xU_05R, xU_c, xU_05L05T, xU_05R05B, xU_05R05T, xU_1R;
			double yU_05T, yU_c, yU_05L05T, yU_05R05B, yU_05R05T, yU_1T;
			double tau_05R, tau_05T;
			double ForceX_05R;
			double ForceY_05T;
			double PhiY_05T;
			double PhiX_05R;

			int k = i*Ny + j;

			// Calculating indices R - right, L - left, T - top, B - bottom, and corners. 
			int k_R = k + Ny;
			int k_L = k - Ny;
			int k_T = k + 1;
			int k_B = k - 1;
			int k_RT = k_R + 1;
			int k_RB = k_R - 1;
			int k_LT = k_L + 1;

			if (i == 0)
			{
				k_L += Ny;
				k_LT += Ny;
			}
			if (j == 0)
			{
				k_B += 1;
				k_RB += 1;
			}

			// Values in the cell centers
			H_05R05T = 0.25*(H[k] + H[k_R] + H[k_T] + H[k_RT]);
			H_05R05B = 0.25*(H[k] + H[k_R] + H[k_B] + H[k_RB]);
			H_05L05T = 0.25*(H[k] + H[k_L] + H[k_T] + H[k_LT]);

			// Values on the cell ribs 
			H_05R = 0.5*(H[k] + H[k_R]);
			H_05T = 0.5*(H[k] + H[k_T]);

			// Right/Left/Top/Bottom value
			H_1R = H[k_R];
			H_1T = H[k_T];

			// Center-value
			H_c = H[k];			

			// Bathymetry B
			B_c = B[k];

			B_1R = B[k_R];
			B_1T = B[k_T];

			// Velocity-x Ux
			xU_c = xU[k];

			xU_05R05T = 0.25*(xU[k] + xU[k_R] + xU[k_T] + xU[k_RT]);
			xU_05R05B = 0.25*(xU[k] + xU[k_R] + xU[k_B] + xU[k_RB]);
			xU_05L05T = 0.25*(xU[k] + xU[k_L] + xU[k_T] + xU[k_LT]);

			xU_05R = 0.5*(xU[k] + xU[k_R]);

			xU_1R = xU[k_R];

			// Velosity-y Uy
			yU_c = yU[k];

			yU_05R05T = 0.25*(yU[k] + yU[k_R] + yU[k_T] + yU[k_RT]);
			yU_05R05B = 0.25*(yU[k] + yU[k_R] + yU[k_B] + yU[k_RB]);
			yU_05L05T = 0.25*(yU[k] + yU[k_L] + yU[k_T] + yU[k_LT]);

			yU_05T = 0.5*(yU[k] + yU[k_T]);

			yU_1T = yU[k_T];

			// tau
			tau_05R = 0.5*(tau[k] + tau[k_R]);
			tau_05T = 0.5*(tau[k] + tau[k_T]);


			ForceX_05R = 0.5*(ForceX[k] + ForceX[k_R]);
			ForceY_05T = 0.5*(ForceY[k] + ForceY[k_T]);


			PhiX_05R = 0.5*(PhiX[k] + PhiX[k_R]);
			PhiY_05T = 0.5*(PhiY[k] + PhiY[k_T]);

			// W
			xW_05R = (tau_05R)*(
				(H_1R * xU_1R * xU_1R - H_c * xU_c * xU_c) / hx
				+ (H_05R05T * xU_05R05T * yU_05R05T - H_05R05B * xU_05R05B * yU_05R05B) / hy
				+ gc * H_05R * ((H_1R + B_1R) - (H_c + B_c)) / hx
				- F_reg * H_05R * ForceX_05R
				- Phi_reg * PhiX_05R
				); //**APPROVED**//
			
			yW_05T = (tau_05T)*(
				(H_05R05T * xU_05R05T * yU_05R05T - H_05L05T * xU_05L05T * yU_05L05T) / hx
				+ (H_1T * yU_1T * yU_1T - H_c * yU_c * yU_c) / hy
				+ gc * H_05T * ((H_1T + B_1T) - (H_c + B_c)) / hy
				- F_reg * H_05T * ForceY_05T
				- Phi_reg * PhiY_05T
				); //**APPROVED**//

			//xJ //yJ_05T = H_05T*yU_05T - yW_05T; //**APPROVED**//
			xJ[k] = (dT / hx) * (H_05R*xU_05R - xW_05R);

			//yJ //yJ_05B = H_05B*yU_05B - yW_05B; //**APPROVED**//
			yJ[k] = (dT / hy) * (H_05T*yU_05T - yW_05T);

			double M = (xJ[k] - xJ[k_L]) + (yJ[k] - yJ[k_B]);
			
			if (massFluxCorrection)
			{
				if (i > 0 && j > 0 && (H_c - M) < 0.0 && !epsilon[k])
				{
					double dT0 = (min(epsFlux, eps/10) - H_c) / M;
					#pragma omp atomic
					xJ[k_L] *= dT0;
					#pragma omp atomic
					yJ[k_B] *= dT0;
					#pragma omp atomic
					xJ[k] *= dT0;
					#pragma omp atomic
					yJ[k] *= dT0;
				}
			}
		}
	}
	
	// calculate internal points
	#pragma omp parallel for
	for (int k = Ny + 1; k < (Nx - 2)*Ny + Ny - 2 + 1; k++)
	{ 
		// Internal points
		if (S[k] == INTERNAL)
		{	
			// Resetting of additional variables

			double xWS_05R, xWS_05L, xWS_05B, xWS_05T;
			double yWS_05R, yWS_05L, yWS_05B, yWS_05T;
			double RS_05R, RS_05L, RS_05B, RS_05T;
			double xxPT_05R, xxPT_05L, xyPT_05R, xyPT_05L, yxPT_05T, yxPT_05B, yyPT_05T, yyPT_05B;
			double xJ_05R, xJ_05L, yJ_05B, yJ_05T;
			double B_05R, B_05L, B_05B, B_05T, B_c, B_05L05B, B_05L05T, B_05R05B, B_05R05T, B_1R, B_1L, B_1B, B_1T;
			double H_05R, H_05L, H_05B, H_05T, H_c, H_05L05B, H_05L05T, H_05R05B, H_05R05T, H_1R, H_1L, H_1B, H_1T;
			double xU_05R, xU_05L, xU_05B, xU_05T, xU_c, xU_05L05B, xU_05L05T, xU_05R05B, xU_05R05T, xU_1R, xU_1L, xU_1B, xU_1T;
			double yU_05R, yU_05L, yU_05B, yU_05T, yU_c, yU_05L05B, yU_05L05T, yU_05R05B, yU_05R05T, yU_1R, yU_1L, yU_1B, yU_1T;

			double ForceX_05R, ForceX_05L, ForceX_05B, ForceX_05T, ForceX_c; 
			double ForceY_05R, ForceY_05L, ForceY_05B, ForceY_05T, ForceY_c; 
			double PhiY_05R, PhiY_05L, PhiY_05B, PhiY_05T, PhiY_c; 
			double PhiX_05R, PhiX_05L, PhiX_05B, PhiX_05T, PhiX_c; 
			double tau_05R, tau_05L, tau_05B, tau_05T, tau_c; 

			// Calculating indices R - right, L - left, T - top, B - bottom, and corners. 
			int k_R = k + Ny;
			int k_L = k - Ny;
			int k_T = k + 1;
			int k_B = k - 1;
			int k_RT = k_R + 1;
			int k_RB = k_R - 1;
			int k_LT = k_L + 1;
			int k_LB = k_L - 1;

			// Values in the cell centers
			H_05R05T = 0.25*(H[k] + H[k_R] + H[k_T] + H[k_RT]);
			H_05R05B = 0.25*(H[k] + H[k_R] + H[k_B] + H[k_RB]);
			H_05L05T = 0.25*(H[k] + H[k_L] + H[k_T] + H[k_LT]);
			H_05L05B = 0.25*(H[k] + H[k_L] + H[k_B] + H[k_LB]);

			// Values on the cell ribs 
			H_05R = 0.5*(H[k] + H[k_R]);
			H_05L = 0.5*(H[k] + H[k_L]);
			H_05T = 0.5*(H[k] + H[k_T]);
			H_05B = 0.5*(H[k_B] + H[k]);
			
			// Right/Left/Top/Bottom value
			H_1R = H[k_R];
			H_1L = H[k_L];
			H_1T = H[k_T];
			H_1B = H[k_B];

			// Center-value
			H_c = H[k];
			B_c = B[k];

			// Bathymetry B
			B_05R05T = 0.25*(B[k] + B[k_R] + B[k_T] + B[k_RT]);
			B_05R05B = 0.25*(B[k] + B[k_R] + B[k_B] + B[k_RB]);
			B_05L05T = 0.25*(B[k] + B[k_L] + B[k_T] + B[k_LT]);
			B_05L05B = 0.25*(B[k] + B[k_L] + B[k_B] + B[k_LB]);

			B_05R = 0.5*(B[k] + B[k_R]);
			B_05L = 0.5*(B[k] + B[k_L]);
			B_05T = 0.5*(B[k] + B[k_T]);
			B_05B = 0.5*(B[k_B] + B[k]);

			B_1R = B[k_R];
			B_1L = B[k_L];
			B_1T = B[k_T];
			B_1B = B[k_B];

			// Velocity-x Ux
			xU_c = xU[k];

			xU_05R05T = 0.25*(xU[k] + xU[k_R] + xU[k_T] + xU[k_RT]);
			xU_05R05B = 0.25*(xU[k] + xU[k_R] + xU[k_B] + xU[k_RB]);
			xU_05L05T = 0.25*(xU[k] + xU[k_L] + xU[k_T] + xU[k_LT]);
			xU_05L05B = 0.25*(xU[k] + xU[k_L] + xU[k_B] + xU[k_LB]);

			xU_05R = 0.5*(xU[k] + xU[k_R]);
			xU_05L = 0.5*(xU[k_L] + xU[k]);
			xU_05T = 0.5*(xU[k] + xU[k_T]);
			xU_05B = 0.5*(xU[k_B] + xU[k]);

			xU_1R = xU[k_R];
			xU_1L = xU[k_L];
			xU_1T = xU[k_T];
			xU_1B = xU[k_B];

			// Velosity-y Uy
			yU_c = yU[k];

			yU_05R05T = 0.25*(yU[k] + yU[k_R] + yU[k_T] + yU[k_RT]);
			yU_05R05B = 0.25*(yU[k] + yU[k_R] + yU[k_B] + yU[k_RB]);
			yU_05L05T = 0.25*(yU[k] + yU[k_L] + yU[k_T] + yU[k_LT]);
			yU_05L05B = 0.25*(yU[k] + yU[k_L] + yU[k_B] + yU[k_LB]);

			yU_05R = 0.5*(yU[k] + yU[k_R]);
			yU_05L = 0.5*(yU[k_L] + yU[k]);
			yU_05T = 0.5*(yU[k] + yU[k_T]);
			yU_05B = 0.5*(yU[k_B] + yU[k]);

			yU_1R = yU[k_R];
			yU_1L = yU[k_L];
			yU_1T = yU[k_T];
			yU_1B = yU[k_B];

			ForceX_c = ForceX[k];

			ForceX_05R = 0.5*(ForceX[k] + ForceX[k_R]);
			ForceX_05L = 0.5*(ForceX[k_L] + ForceX[k]);
			ForceX_05T = 0.5*(ForceX[k] + ForceX[k_T]);
			ForceX_05B = 0.5*(ForceX[k_B] + ForceX[k]);

			ForceY_c = ForceY[k];

			ForceY_05R = 0.5*(ForceY[k] + ForceY[k_R]);
			ForceY_05L = 0.5*(ForceY[k_L] + ForceY[k]);
			ForceY_05T = 0.5*(ForceY[k] + ForceY[k_T]);
			ForceY_05B = 0.5*(ForceY[k_B] + ForceY[k]);

			PhiX_c = PhiX[k];

			PhiX_05R = 0.5*(PhiX[k] + PhiX[k_R]);
			PhiX_05L = 0.5*(PhiX[k_L] + PhiX[k]);
			PhiX_05T = 0.5*(PhiX[k] + PhiX[k_T]);
			PhiX_05B = 0.5*(PhiX[k_B] + PhiX[k]);

			PhiY_c = PhiY[k];

			PhiY_05R = 0.5*(PhiY[k] + PhiY[k_R]);
			PhiY_05L = 0.5*(PhiY[k_L] + PhiY[k]);
			PhiY_05T = 0.5*(PhiY[k] + PhiY[k_T]);
			PhiY_05B = 0.5*(PhiY[k_B] + PhiY[k]);

			// tau
			tau_c = tau[k];

			tau_05R = 0.5*(tau[k] + tau[k_R]);
			tau_05L = 0.5*(tau[k_L] + tau[k]);
			tau_05T = 0.5*(tau[k] + tau[k_T]);
			tau_05B = 0.5*(tau[k_B] + tau[k]);				

			//xJ
			xJ_05R = (hx / dT) * xJ[k]; //H_05R*xU_05R - xW_05R; //**APPROVED**//
			xJ_05L = (hx / dT) *xJ[k_L]; ////H_05L*xU_05L - xW_05L; //**APPROVED**//

			//yJ
			yJ_05T = (hy / dT) *yJ[k];// H_05T*yU_05T - yW_05T; //**APPROVED**//
			yJ_05B = (hy / dT) *yJ[k_B]; // H_05B*yU_05B - yW_05B; //**APPROVED**//

			// Next time-step H

			Ht[k] = H_c - (dT / hx)*(xJ_05R - xJ_05L) - (dT / hy) * (yJ_05T - yJ_05B); //**APPROVED**//
			if (Ht[k] < 0.0 && massFluxCorrection)
				Ht[k] = eps / 2;
					
			// Dry-zone condition
			if (Ht[k]>eps && !epsilon[k])//
			{
				
				//xWS
				xWS_05R = tau_05R * H_05R * (
					xU_05R *(xU_1R - xU_c) / hx 
					+ yU_05R * (xU_05R05T - xU_05R05B) / hy
					+ gc * ((H_1R + B_1R) - (H_c + B_c)) / hx 
					- F_reg * ForceX_05R
					) 
					- Phi_reg * tau_05R * PhiX_05R; //**APPROVED**//
				xWS_05L = tau_05L * H_05L * (
					xU_05L *(xU_c - xU_1L) / hx 
					+ yU_05L * (xU_05L05T - xU_05L05B) / hy
					+ gc * ((H_c + B_c) - (H_1L + B_1L)) / hx 
					- F_reg * ForceX_05L
					) 
					- Phi_reg * tau_05L * PhiX_05L; //**APPROVED**//
				xWS_05T = tau_05T * H_05T * (
					xU_05T *(xU_05R05T - xU_05L05T) / hx 
					+ yU_05T * (xU_1T - xU_c) / hy
					+ gc * ((H_05R05T + B_05R05T) - (H_05L05T + B_05L05T)) / hx 
					- F_reg * ForceX_05T
					) 
					- Phi_reg * tau_05T * PhiX_05T; //**APPROVED**//
				xWS_05B = tau_05B * H_05B * (
					xU_05B *(xU_05R05B - xU_05L05B) / hx 
					+ yU_05B * (xU_c - xU_1B) / hy
					+ gc * ((H_05R05B + B_05R05B) - (H_05L05B + B_05L05B)) / hx 
					- F_reg * ForceX_05B
					) 
					- Phi_reg * tau_05B * PhiX_05B; //**APPROVED**//

				//yWS
				yWS_05R = tau_05R * H_05R * (
					xU_05R *(yU_1R - yU_c) / hx 
					+ yU_05R * (yU_05R05T - yU_05R05B) / hy
					+ gc * ((H_05R05T + B_05R05T) - (H_05R05B + B_05R05B)) / hy 
					- F_reg * ForceY_05R
					) 
					- Phi_reg * tau_05R * PhiY_05R; //**APPROVED**//
				yWS_05L = tau_05L * H_05L * (
					xU_05L *(yU_c - yU_1L) / hx 
					+ yU_05L * (yU_05L05T - yU_05L05B) / hy
					+ gc * ((H_05L05T + B_05L05T) - (H_05L05B + B_05L05B)) / hy 
					- F_reg * ForceY_05L
					) 
					- Phi_reg * tau_05L * PhiY_05L; //**APPROVED**//
				yWS_05T = tau_05T * H_05T * (
					xU_05T *(yU_05R05T - yU_05L05T) / hx 
					+ yU_05T * (yU_1T - yU_c) / hy
					+ gc * ((H_1T + B_1T) - (H_c + B_c)) / hy 
					- F_reg * ForceY_05T
					) 
					- Phi_reg * tau_05T * PhiY_05T; //**APPROVED**//
				yWS_05B = tau_05B * H_05B * (
					xU_05B *(yU_05R05B - yU_05L05B) / hx 
					+ yU_05B * (yU_c - yU_1B) / hy
					+ gc * ((H_c + B_c) - (H_1B + B_1B)) / hy 
					- F_reg * ForceY_05B
					) 
					- Phi_reg * tau_05B * PhiY_05B; //**APPROVED**//

				//RS
				RS_05R = gc * tau_05R * (
					xU_05R * 0.5 * (H_1R * H_1R - H_c * H_c) / hx
					+ yU_05R * 0.5 * (H_05R05T * H_05R05T - H_05R05B * H_05R05B) / hy
					+ H_05R * H_05R * ((xU_1R - xU_c) / hx + (yU_05R05T - yU_05R05B) / hy) 
					); //**APPROVED**//
				RS_05L = gc * tau_05L * (
					xU_05L * 0.5 * (H_c * H_c - H_1L * H_1L) / hx
					+ yU_05L * 0.5 * (H_05L05T * H_05L05T - H_05L05B * H_05L05B) / hy
					+ H_05L * H_05L * ((xU_c - xU_1L) / hx + (yU_05L05T - yU_05L05B) / hy)
					); //**APPROVED**//
				RS_05T = gc * tau_05T * (
					xU_05T * 0.5 * (H_05R05T * H_05R05T - H_05L05T * H_05L05T) / hx
					+ yU_05T * 0.5 * (H_1T * H_1T - H_c * H_c) / hy
					+ H_05T * H_05T * ((xU_05R05T - xU_05L05T) / hx + (yU_1T - yU_c) / hy) 
					); //**APPROVED**//
				RS_05B = gc * tau_05B * (
					xU_05B * 0.5 * (H_05R05B * H_05R05B - H_05L05B * H_05L05B) / hx
					+ yU_05B * 0.5 * (H_c * H_c - H_1B * H_1B) / hy
					+ H_05B * H_05B * ((xU_05R05B - xU_05L05B) / hx + (yU_c - yU_1B) / hy) 
					); //**APPROVED**//
				
				////xxPT
				xxPT_05R = xU_05R * xWS_05R + RS_05R + NS * 2 * gc * tau_05R * H_05R * H_05R * 0.5 * (xU_1R - xU_c) / hx; //**APPROVED**//
				xxPT_05L = xU_05L * xWS_05L + RS_05L + NS * 2 * gc * tau_05L * H_05L * H_05L * 0.5 * (xU_c - xU_1L) / hx; //**APPROVED**//

				////yyPT
				yyPT_05T = yU_05T * yWS_05T + RS_05T + NS * 2 * gc * tau_05T * H_05T * H_05T * 0.5 * (yU_1T - yU_c) / hy; //**APPROVED**//
				yyPT_05B = yU_05B * yWS_05B + RS_05B + NS * 2 * gc * tau_05B * H_05B * H_05B * 0.5 * (yU_c - yU_1B) / hy; //**APPROVED**//

				////yxPT
				yxPT_05T = yU_05T * xWS_05T + NS * gc * tau_05T * H_05T * H_05T * 0.5 * ((xU_1T - xU_c) / hy + (yU_05R05T - yU_05L05T) / hx); //**APPROVED**//
				yxPT_05B = yU_05B * xWS_05B + NS * gc * tau_05B * H_05B * H_05B * 0.5 * ((xU_c - xU_1B) / hy + (yU_05R05B - yU_05L05B) / hx); //**APPROVED**//

				////xyPT
				xyPT_05R = xU_05R * yWS_05R + NS * gc * tau_05R * H_05R * H_05R * 0.5 * ((xU_05R05T - xU_05R05B) / hy + (yU_1R - yU_c) / hx); //**APPROVED**//
				xyPT_05L = xU_05L * yWS_05L + NS * gc * tau_05L * H_05L * H_05L * 0.5 * ((xU_05L05T - xU_05L05B) / hy + (yU_c - yU_1L) / hx); //**APPROVED**//
								
				////xUt	
				//// Well balanced: H[k] _c 0.5*(H2y2+H2y1)
				xUt[k] = (H_c * xU_c
					+ dT * PhiX_c
					+ (dT / hx) * ((xxPT_05R - xxPT_05L) - (xU_05R * xJ_05R - xU_05L * xJ_05L))
					+ (dT / hy) * ((yxPT_05T - yxPT_05B) - (xU_05T * yJ_05T - xU_05B * yJ_05B))
					- 0.5 * gc * (dT / hx) * (H_05R + H_05L) * (
						ForceX_c * hx / (-gc)
						+ 0.5 * ((H_1R + B_1R) - (H_1L + B_1L))
						) 
					- dT * tau_c * ((H_05R * xU_05R - H_05L * xU_05L) / hx
						+ (H_05T * yU_05T - H_05B * yU_05B) / hy
						) * (ForceX_c - gc * (B_05R - B_05L) / hx)
					) / Ht[k]; //**APPROVED**//
				
				////yUt
				//// Well balanced: H[k] _c 0.5*(H2y2+H2y1)
				yUt[k] = (H_c * yU_c
					+ dT * PhiY_c
					+ (dT / hx) * ((xyPT_05R - xyPT_05L) - (yU_05R * xJ_05R - yU_05L * xJ_05L))
					+ (dT / hy) * ((yyPT_05T - yyPT_05B) - (yU_05T * yJ_05T - yU_05B * yJ_05B))
					- 0.5 * gc * (dT / hy) * (H_05T + H_05B) * (
						ForceY_c * hy / (-gc)
						+ 0.5 * ((H_1T + B_1T) - (H_1B + B_1B))
						)
					- dT * tau_c * ((H_05R * xU_05R - H_05L * xU_05L) / hx
						+ (H_05T * yU_05T - H_05B * yU_05B) / hy
						) * (ForceY_c - gc * (B_05T - B_05B) / hy)
					) / Ht[k]; //**APPROVED**//
			}
			else
			{
				xUt[k] = 0.0;
				yUt[k] = 0.0;
			}	
				
			if (TransportProblemFlag)
			{
				double C_05R, C_05L, C_05B, C_05T, C_c, C_05L05B, C_05L05T, C_05R05B, C_05R05T, C_1R, C_1L, C_1B, C_1T;
				C_05R05T = 0.25*(C[k] + C[k_R] + C[k_T] + C[k_RT]);
				C_05R05B = 0.25*(C[k] + C[k_R] + C[k_B] + C[k_RB]);
				C_05L05T = 0.25*(C[k] + C[k_L] + C[k_T] + C[k_LT]);
				C_05L05B = 0.25*(C[k] + C[k_L] + C[k_B] + C[k_LB]);

				C_05R = 0.5*(C[k] + C[k_R]);
				C_05L = 0.5*(C[k_L] + C[k]);
				C_05T = 0.5*(C[k] + C[k_T]);
				C_05B = 0.5*(C[k_B] + C[k]);

				C_1R = C[k_R];
				C_1L = C[k_L];
				C_1T = C[k_T];
				C_1B = C[k_B];

				C_c = C[k];
				// Concentration equation
				if (Ht[k] > eps && !epsilon[k])
				{
					Ct[k] = (C_c*H_c - (dT / hx)*(xJ_05R*C_05R - xJ_05L*C_05L) - (dT / hy)*(yJ_05T*C_05T - yJ_05B*C_05B)
						+ (dT / hx)*(H_05R*((C_1R - C_c) / hx * (D + NSC*gc*tau_05R*H_05R + alpha_c * tau_05R*xU_05R*xU_05R) + alpha_c * tau_05R*xU_05R*yU_05R * (C_05R05T - C_05R05B) / hy)
							- H_05L*((C_c - C_1L) / hx * (D + NSC*gc*tau_05L*H_05L + alpha_c * tau_05L*xU_05L*xU_05L) + alpha_c * tau_05L*xU_05L*yU_05L * (C_05L05T - C_05L05B) / hy))
						+ (dT / hy)*(H_05T*((C_1T - C_c) / hy * (D + NSC*gc*tau_05T*H_05T + alpha_c * tau_05T*yU_05T*yU_05T) + alpha_c * tau_05T*xU_05T*yU_05T * (C_05R05T - C_05L05T) / hx)
							- H_05B*((C_c - C_1B) / hy * (D + NSC*gc*tau_05B*H_05B + alpha_c * tau_05B*yU_05B*yU_05B) + alpha_c * tau_05B*xU_05B*yU_05B * (C_05R05B - C_05L05B) / hx))
						) / Ht[k];
				}
				else
				{
					Ct[k] = C_c;
				}
			}
		} // end if (TypeOfPoint == INTERNAL)
	} //end for (int k = Ny + 1; k < (Nx - 2)*Ny + Ny - 2 + 1; k++)
	
	// Boundary conditions
	if (BoundaryConditionsFromFile)
		RecalcFileBoundaryConditions();

	if (!InternalWallsFlag)
	{
		// outer boundaries
		#pragma omp parallel for
		for (int m = 0; m < 2 * Nx + 2 * Ny; m++)
		{
			int i = 0;
			int j = 0;

			if (m<Ny)
			{
				i = 0;
				j = m;
			}
			else if ((m >= Ny) && (m < Ny + Nx))
			{
				i = m - Ny;
				j = Ny - 1;
			}
			else if ((m >= Ny + Nx) && (m < Ny + 2 * Nx))
			{
				i = m - Ny - Nx;
				j = 0;
			}
			else
			{
				i = Nx - 1;
				j = m - 2 * Nx - Ny;
			}

			int n = i*Ny + j;

			int type = S[n];
			

			int k = n + 1 * (type == BOTTOM) - Ny*(type == RIGHT) - 1 * (type == TOP) + Ny*(type == LEFT) + (Ny + 1)*(type == LB_CORNER) + (1 - Ny)*(type == RB_CORNER) + (-Ny - 1)*(type == RT_CORNER) + (Ny - 1)*(type == LT_CORNER);
			
			double boundaryForce = 0.0;
			
			if (F_bound)
			{
			    boundaryForce = 
				- (type == BOTTOM) * 0.5 * (ForceX[n] + ForceX[k]) * hx / gc
				+ (type == TOP) * 0.5 * (ForceX[n] + ForceX[k]) * hx / gc
				- (type == LEFT) * 0.5 * (ForceY[n] + ForceY[k]) * hy / gc
				+ (type == RIGHT) * 0.5 * (ForceY[n] + ForceY[k]) * hy / gc;
			}

			if (Ht[k]  > eps && !epsilon[n])
			{

				if (border[VELOCITY_X][type] != FROM_FILE)
					xUt[n] = border[VELOCITY_X][type] * xUt[k] + 2 * border_C[VELOCITY_X][type];
				if (border[VELOCITY_Y][type] != FROM_FILE)
					yUt[n] = border[VELOCITY_Y][type] * yUt[k] + 2 * border_C[VELOCITY_Y][type];
				if (border[HEIGHT][type] != FROM_FILE)
					Ht[n] = border[HEIGHT][type] * (Ht[k] + B[k]) - B[n] + boundaryForce + 2 * border_C[HEIGHT][type];

				if (TransportProblemFlag)
				{
					if (border[CONCENTRATION][type] != FROM_FILE)
						Ct[n] = border[CONCENTRATION][type] * Ct[k] + 2 * border_C[CONCENTRATION][type];
				}
				// COSTIL! COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL! 

				if (type == TOP)
				{
					yUt[n] = (2.0) * sin(pi * (Time_elapsed + dT) / (3600.0 * 6.0)) * (1.0 - exp(-(Time_elapsed + dT) / (6*3600.0) )) - yUt[k];
					Ct[n] = 2.0 * (1.0 - exp(-(Time_elapsed + dT) / 300.0)) - Ct[k];
				}
				
				// COSTIL! COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL!  COSTIL! 
			}
			else
			{
				xUt[n] = 0;
				yUt[n] = 0;

				Ht[n] = H[n];// Ht[k];
				if (TransportProblemFlag)
				{
					Ct[n] = C[n];
				}
			}
		}
	}
	else
	{
		// outer boundaries with internal boundaries
		#pragma omp parallel for
		for (int n = 0; n < Nx * Ny; n++)
		{
			int type = S[n];
			
			int k; 
			if (type >= INTERNALWALL)
			{
				type = type - INTERNALWALL;
				
				k = n + 1 * (type == BOTTOM) - Ny*(type == RIGHT) - 1 * (type == TOP) + Ny*(type == LEFT) + (Ny + 1)*(type == LB_CORNER) + (1 - Ny)*(type == RB_CORNER) + (-Ny - 1)*(type == RT_CORNER) + (Ny - 1)*(type == LT_CORNER);
				

				if (Ht[k]  > eps && !epsilon[n])
				{
					xUt[n] = border_WALL[VELOCITY_X][type] * xUt[k];
					yUt[n] = border_WALL[VELOCITY_Y][type] * yUt[k];
					Ht[n] =  (Ht[k] + B[k]) - B[n];
					if (TransportProblemFlag)
					{
						Ct[n] = Ct[k];
					}

				}
				else
				{
					Ht[n] = H[n];// Ht[k];
					if (TransportProblemFlag)
					{
						Ct[n] = C[n];
					}
				}
			}
			else if (type < INTERNAL && type != EXCLUDED)
			{
				
				k = n + 1 * (type == BOTTOM) - Ny*(type == RIGHT) - 1 * (type == TOP) + Ny*(type == LEFT) + (Ny + 1)*(type == LB_CORNER) + (1 - Ny)*(type == RB_CORNER) + (-Ny - 1)*(type == RT_CORNER) + (Ny - 1)*(type == LT_CORNER);

				if (Ht[k]  > eps/* && !epsilon[k]*/)
				{
					if (border[VELOCITY_X][type] != FROM_FILE)
						xUt[n] = border[VELOCITY_X][type] * xUt[k] + 2 * border_C[VELOCITY_X][type];
					if (border[VELOCITY_Y][type] != FROM_FILE)
						yUt[n] = border[VELOCITY_Y][type] * yUt[k] + 2 * border_C[VELOCITY_Y][type];
					if (border[HEIGHT][type] != FROM_FILE)
						Ht[n] = border[HEIGHT][type] * (Ht[k] + B[k]) - B[n] + 2 * border_C[HEIGHT][type];
					if (TransportProblemFlag)
					{
						if (border[CONCENTRATION][type] != FROM_FILE)
							Ct[n] = border[CONCENTRATION][type] * Ct[k] + 2 * border_C[CONCENTRATION][type];
					}
				}
				else
				{
					xUt[n] = 0;
					yUt[n] = 0;

					Ht[n] = H[n];// Ht[k];
					if (TransportProblemFlag)
					{
						Ct[n] = C[n];
					}
				}
			}
		}
	}

	if (tidesHarmonics)
	{
		addTidesHarmonicsBoundaryConditions();
	}
	
	// set array of epsilons to zero
	memset(epsilon, 0, Nx*Ny * sizeof(int)); 

	#pragma omp parallel for
	for (int m = 0; m < Nx*Ny; m++)
	{
		int  i = int(m / Ny);
		int j = m % Ny;

		H[m] = Ht[m];

		if (TransportProblemFlag)
		{
			C[m] = Ct[m];
		}

		// Dry-wet condition 
		// if H < eps -> point is dry
		// if (H+B)_dry > (H+B)_wet -> no mass flux, no momentum flux 
		// if (H+B)_dry < (H+B)_wet -> only mass flux

		if (H[m] <= eps ) {

			#pragma omp atomic
			epsilon[m] += 1;

			if (i != Nx - 1) {
				#pragma omp atomic
				epsilon[(i + 1)*Ny + j] += (int)(Ht[(i + 1)*Ny + j] + B[(i + 1)*Ny + j] < B[m] + eps);
				//xUt[(i + 1)*Ny + j] = 0.0;
				//yUt[(i + 1)*Ny + j] = 0.0;
			}
			if (i != 0) {
				#pragma omp atomic
				epsilon[(i - 1)*Ny + j] += (int)(Ht[(i - 1)*Ny + j] + B[(i - 1)*Ny + j] < B[m] + eps);
				//xUt[(i - 1)*Ny + j] = 0.0;
				//yUt[(i - 1)*Ny + j] = 0.0;
			}
			if (j != Ny - 1) {
				#pragma omp atomic
				epsilon[i*Ny + j + 1] += (int)(Ht[i*Ny + j + 1] + B[i*Ny + j + 1] < B[m] + eps);
				//xUt[i*Ny + j + 1] = 0.0;
				//yUt[i*Ny + j + 1] = 0.0;
			}
					
			if (j != 0) {
				#pragma omp atomic
				epsilon[i*Ny + j - 1] += (int)(Ht[i*Ny + j - 1] + B[i*Ny + j - 1] < B[m] + eps);
				//xUt[i*Ny + j - 1] = 0.0;
				//yUt[i*Ny + j - 1] = 0.0;
			}	
		}
		else
		{
			bool noMomentumFlux = false;

			if (i != Nx - 1) {
				noMomentumFlux = noMomentumFlux || (Ht[(i + 1)*Ny + j] < eps);
			}
			if (i != 0) {
				noMomentumFlux = noMomentumFlux || (Ht[(i - 1)*Ny + j] < eps);
			}
			if (j != Ny - 1) {
				noMomentumFlux = noMomentumFlux || (Ht[i*Ny + j + 1] < eps);
			}
			if (j != 0) {
				noMomentumFlux = noMomentumFlux || (Ht[i*Ny + j - 1] < eps);
			}

			if (!noMomentumFlux)
			{
				xU[m] = xUt[m];
				yU[m] = yUt[m];
			}
			else
			{
				xU[m] = 0.0;
				yU[m] = 0.0;
			}
		}
		
	}//end of for (m=0; m<(Ny*Nx); m++)
	
	// calculation of tau
	#pragma omp parallel for
	for (int m = 0; m < Nx*Ny; m++)
	{ 		
		if (!epsilon[m])
		{
			tau[m] = alpha*sqrt(hx*hy) / sqrt(gc*H[m]);//alpha*sqrt(hx*hy) / sqrt(gc*H[m]);
		}
		else
		{
			tau[m] = 0.0;
			xU[m] = 0.0;
			yU[m] = 0.0;
			//tau[m] = alpha*sqrt(hx*hy) / sqrt(gc*eps/2);
		}

		// checking if slgorithm become unstable
		double check = fabs(H[m]) + fabs(xU[m]) + fabs(yU[m]) + fabs(tau[m]);

		if (check > CriticalVal || check != check)
		{
			cout << "Critical Value: " << CriticalVal << " " << check << endl;

			Print_info_about_point("ERROR POINT", m);
			pause();
			/*
			Print_info_about_point("L:ERROR POINT", (i - 1)*Ny + j);
			Print_info_about_point("R:ERROR POINT", (i + 1)*Ny + j);
			Print_info_about_point("U:ERROR POINT", i*Ny + j + 1);
			Print_info_about_point("D:ERROR POINT", i*Ny + j - 1);
			Print_info_about_point("RU:ERROR POINT", (i + 1)*Ny + j + 1);
			Print_info_about_point("LU:ERROR POINT", (i - 1)*Ny + j + 1);
			Print_info_about_point("RD:ERROR POINT", (i + 1)*Ny + j - 1);
			Print_info_about_point("LD:ERROR POINT", (i - 1)*Ny + j - 1);
			*/
			//Visualization_to_techplot_fstream();
			//pause();
			//Raschet::Visualization_to_techplot();
			Stop_Raschet_Flag = 1;
		}
	}//end of for (m=0; m<(Ny*Nx); m++)
}
