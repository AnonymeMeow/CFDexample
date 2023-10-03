/*
 * thermochem.c
 *
 *  Created on: 2014.01.14
 *
 */

#include<math.h>

#include"comm.h"

/*---------------------------------------------------
 * Calculate the thermal properties
 * ------------------------------------------------*/
void gettherm(int nc, double **q, double *p,
		      double *t, double *gam1, double *rgas1, double *cv1)
{
	int ic;
	double temp, rgas0;

	if(config1.gasModel == 0)
	{
		/*-- for calorically perfect gas--*/
		rgas0 = (ru/config2.molWeight);
		for(ic=0; ic<nc; ic++)
		{
			temp = 0.5*(q[ic][1]*q[ic][1] + q[ic][2]*q[ic][2]);
			rgas1[ic] = rgas0;
			gam1[ic] = config2.gam0;
			cv1[ic]  = rgas0/(config2.gam0 - 1.);
			p[ic] = (gam1[ic] - 1)*(q[ic][3] - temp)*q[ic][0];
			t[ic] = p[ic]/q[ic][0]/rgas1[ic];
		}
	}
}
/*--------------------------------------------------------------
 * Calculate transport properties(viscosity, conductivity, ...
 * -------------------------------------------------------------*/
void gettrans(int nc, double *t, double *gam1,
		      double *cv1, double *mu, double *cond)
{
	int ic;
	double cp, ratio1, ratio2, tem;

	if(config1.gasModel == 0)
	{
		if(config1.visModel == 2)
		{
			/*---------- 1.2 Sutherland's law ----------*/
			for(ic=0; ic<nc; ic++)
			{
				tem = t[ic];
				cp  = gam1[ic]*cv1[ic];
				ratio1 = tem/config2.suthC1;
				ratio2 = (config2.suthC1 + config2.suthC2)/(tem + config2.suthC2);
				mu[ic] = config2.muRef*ratio2*pow(ratio1,1.5);
				cond[ic] = mu[ic]*cp/config2.Pr0;
			}
			for(ic=0; ic<nc; ic++)
			{
				mu[ic]   = mu[ic];
				cond[ic] = cond[ic];
			}
		}
	}
}