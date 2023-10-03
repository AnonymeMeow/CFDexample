/*
 * tmarch.c
 *
 *  Created on: Jul 9, 2014
 *  Purpose: conduct the time march step
 *
 */

#include"comm.h"

/*------------------------------------------------------
 * update the conservative variables for perfect-gas
 * use Runge-kutta method
 * ------------------------------------------------*/
void updateq(int nc, double **q, double **rhs, double dt, int idt)
{
	int i, j, ii, jj, ic, ic1;
	double rho, rhoU, rhoV, rhoE, rrho, yas;
	double coef[3][3] = {{1., 3./4., 1./3.}, {0, 1./4., 2./3.},
							  {1., 1./4., 2./3.}};

	/* the (:,0), and (config1.ni-1, :) are assigned wall condition */

	for(i=0; i<config1.ni; i++)
	{
		ii = i + config1.Ng;
		for(j=0; j<config1.nj; j++)
		{
			jj  = j + config1.Ng;

			ic  = i*config1.nj + j;
			ic1 = ii*J0 + jj;
			yas = mesh.yaks[ic1];

			rho  = q[ic][0];
			rhoU = q[ic][0]*q[ic][1];
			rhoV = q[ic][0]*q[ic][2];
			rhoE = q[ic][0]*q[ic][3];

			q[ic][0] = coef[0][idt]*qo[ic][0]           + coef[1][idt]*rho  + coef[2][idt]*dt*rhs[ic][0]/yas;
			q[ic][1] = coef[0][idt]*qo[ic][0]*qo[ic][1] + coef[1][idt]*rhoU + coef[2][idt]*dt*rhs[ic][1]/yas;
			q[ic][2] = coef[0][idt]*qo[ic][0]*qo[ic][2] + coef[1][idt]*rhoV + coef[2][idt]*dt*rhs[ic][2]/yas;
			q[ic][3] = coef[0][idt]*qo[ic][0]*qo[ic][3] + coef[1][idt]*rhoE + coef[2][idt]*dt*rhs[ic][3]/yas;

			rrho = 1./q[ic][0];
			q[ic][1] = q[ic][1]*rrho;
			q[ic][2] = q[ic][2]*rrho;
			q[ic][3] = q[ic][3]*rrho;
		}
	}
}