/*
 * vflux.c
 *
 *  Created on:    Apr 4, 2014
 *  last modified: Sep 10, 2014
 */

#include"comm.h"

/*--------------------------------------------------------------
 * Calculate viscous flux in x direction
 * -------------------------------------------------------------*/
void vfluxF(int start, int end)
{
	int i, j, il, ir, jr, ii, jj, k, ik, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, fv[12];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

    idr  = 0;
	idu  = 1;
    idv  = 2;
    idt  = 3;
    idm  = 4;
    idk  = 5;
    idux = 6;
    iduy = 7;
    idvx = 8;
    idvy = 9;
    idtx = 10;
    idty = 11;

    nvar = 12;
	/*--------------------------------------------------*/

   	il = start + config1.Ng - 1;
	ir = end + config1.Ng;
    jr = config1.nj + config1.Ng;

	Ec0 = 1.;
	coef_Re = 1.;
	coef_e  = 1.;

    for(j=config1.Ng; j<jr; j++)
    {
        for(i=start; i<end+2*config1.Ng; i++)
        {
        	/*---- convert to 1D-array ----*/
        	ic = i*J0 + j;

        	U1d.du[i] = Uv.u_xi[ic];
            U1d.dv[i] = Uv.v_xi[ic];
        	U1d.dt[i] = Uv.T_xi[ic];
        	U1d.rho[i]=  U.q[ic][0];
        	U1d.u[i]  =  U.q[ic][1];
        	U1d.v[i]  =  U.q[ic][2];
        	U1d.p[i]  =  U.pre[ic];
        	U1d.t[i]  =  U.tem[ic];
        	U1d.mu[i] =  U.mu[ic];
        	U1d.kt[i] =  U.kt[ic];
        }

    	for(i=il; i<ir; i++)
    	{
    		ic = i*J0 + j;

    		for(ivar=0; ivar<nvar; ivar++)
    			fv[ivar] = 0.;

	        /*---------- 1. get the face value ----------*/

	        for(k=0; k<6; k++)
	        {
	        	ik  = i - 2 + k;
    	/*---------- 1. interpolate i derivatives at the i+1/2 faces ----------*/
	    		fv[idr]  = fv[idr]  + interpo[k]*U1d.rho[ik];
	    		fv[idu]  = fv[idu]  + interpo[k]*U1d.u[ik];
	    		fv[idv]  = fv[idv]  + interpo[k]*U1d.v[ik];
	    		fv[idt]  = fv[idt]  + interpo[k]*U1d.t[ik];
	    		fv[idm]  = fv[idm]  + interpo[k]*U1d.mu[ik];
	    		fv[idk]  = fv[idk]  + interpo[k]*U1d.kt[ik];
	    		fv[iduy] = fv[iduy] + interpo[k]*U1d.du[ik];
	    		fv[idvy] = fv[idvy] + interpo[k]*U1d.dv[ik];
	    		fv[idty] = fv[idty] + interpo[k]*U1d.dt[ik];

        /*---------- 2. approximate j derivatives at i+1/2 faces ----------*/
	    		fv[idux] = fv[idux] + approxi[k]*U1d.u[ik]/dxc;
	    		fv[idvx] = fv[idvx] + approxi[k]*U1d.v[ik]/dxc;
	    		fv[idtx] = fv[idtx] + approxi[k]*U1d.t[ik]/dxc;
	    	}

        /*---------- 3.calculate the viscous Flux ----------*/
    		U1d.flux[i][0] = 0.;
    		U1d.flux[i][1] =
				coef_Re*fv[idm]*(
					Uv.fu1[ic]*fv[idux] +
					Uv.fu2[ic]*fv[iduy] +
					Uv.fuv[ic]*fv[idvx] +
					Uv.fu3[ic]*fv[idvy]
				);
    		U1d.flux[i][2] =
				coef_Re*fv[idm]*(
					Uv.fuv[ic]*fv[idux] +
					Uv.fv1[ic]*fv[iduy] +
					Uv.fv2[ic]*fv[idvx] +
					Uv.fv3[ic]*fv[idvy]
				);
    		U1d.flux[i][3] =
				fv[idu]*U1d.flux[i][1] +
				fv[idv]*U1d.flux[i][2] +
				coef_e*fv[idk]*(
					Uv.fe1[ic]*fv[idtx] +
					Uv.fe2[ic]*fv[idty]
				);
    	}
		/* For the same grids(without transformation), the metrics coefficients are:
		 * Uv.fu1 = 4./3, Uv.fu2 = 0, Uv.fu3 = -2./3, Uv.fuv = 0.;
		 * Uv.fv1 = 1,  Uv.fuv=0, Uv.fv2 = 1, Uv.fv3 = 0.;
		 * Uv.fe1 = 1., Uv.ge2 = 0.
		 * which reduced to the origin Navier-Storks Eq.
		 * Confirmed! */

		jj = j-config1.Ng;
		for(i=start+config1.Ng; i<ir; i++)
		{
			ii = i - config1.Ng;
			ic = ii*config1.nj + jj;
			for(k=0; k<neqv; k++)
				rhs[ic][k] += (U1d.flux[i][k] - U1d.flux[i-1][k])/dxc;
		}
	}
}


/*--------------------------------------------------------------
 * Calculate viscous Flux in j direction
 * -------------------------------------------------------------*/
void vfluxG(int start, int end)
{
	int i, j, jl, ir, jr, ii, jj, k, jk, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, fv[12];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

    idr  = 0;
	idu  = 1;
    idv  = 2;
    idt  = 3;
    idm  = 4;
    idk  = 5;
    idux = 6;
    iduy = 7;
    idvx = 8;
    idvy = 9;
    idtx = 10;
    idty = 11;

    nvar = 12;
	/*--------------------------------------------------*/

   	jl = config1.Ng - 1;
	ir = end + config1.Ng;
    jr = config1.nj + config1.Ng;

	Ec0 = 1.;
	coef_Re = 1.;
	coef_e  = 1.;

    for(i = start + config1.Ng; i<ir; i++)
    {
        for(j=0; j<J0; j++)
        {
        	/*---- convert to 1D-array ----*/
        	ic = i*J0 + j;

        	U1d.du[j] = Uv.u_xi[ic];
            U1d.dv[j] = Uv.v_xi[ic];
        	U1d.dt[j] = Uv.T_xi[ic];
        	U1d.rho[j]=  U.q[ic][0];
        	U1d.u[j]  =  U.q[ic][1];
        	U1d.v[j]  =  U.q[ic][2];
        	U1d.p[j]  =  U.pre[ic];
        	U1d.t[j]  =  U.tem[ic];
        	U1d.mu[j] =  U.mu[ic];
        	U1d.kt[j] =  U.kt[ic];
        }

    	for(j=jl; j<jr; j++)
    	{
    		ic  = i*J0 + j;

    		for(ivar=0; ivar<nvar; ivar++)
    			fv[ivar] = 0.;

    	    for(k = 0; k<6; k++)
    	    {
    	    	jk  = j - 2 + k;

    	/*---------- 1. interpolate at the j+1/2 faces ----------*/
    		    fv[idr]  = fv[idr]  + interpo[k]*U1d.rho[jk];
    		    fv[idu]  = fv[idu]  + interpo[k]*U1d.u[jk];
    		    fv[idv]  = fv[idv]  + interpo[k]*U1d.v[jk];
    		    fv[idt]  = fv[idt]  + interpo[k]*U1d.t[jk];
    			fv[idm]  = fv[idm]  + interpo[k]*U1d.mu[jk];
    			fv[idk]  = fv[idk]  + interpo[k]*U1d.kt[jk];
    		    fv[idux] = fv[idux] + interpo[k]*U1d.du[jk];
    			fv[idvx] = fv[idvx] + interpo[k]*U1d.dv[jk];
    			fv[idtx] = fv[idtx] + interpo[k]*U1d.dt[jk];

        /*---------- 2. approximate j derivatives at j+1/2 faces ----------*/
    			fv[iduy] = fv[iduy] + approxi[k]*U1d.u[jk]/dyc;
    			fv[idvy] = fv[idvy] + approxi[k]*U1d.v[jk]/dyc;
    			fv[idty] = fv[idty] + approxi[k]*U1d.t[jk]/dyc;
    		}

        /*---------- 3.calculate the viscous Flux ----------*/
			U1d.flux[j][0] = 0.;
			U1d.flux[j][1] =
				coef_Re*fv[idm]*(
					Uv.gu1[ic]*fv[idux] +
					Uv.gu2[ic]*fv[iduy] +
					Uv.gu3[ic]*fv[idvx] +
					Uv.guv[ic]*fv[idvy]
				);
			U1d.flux[j][2] =
				coef_Re*fv[idm]*(
					Uv.gv1[ic]*fv[idux] +
					Uv.guv[ic]*fv[iduy] +
					Uv.gv2[ic]*fv[idvx] +
					Uv.gv3[ic]*fv[idvy]
				);
			U1d.flux[j][3] =
				fv[idu]*U1d.flux[j][1] +
				fv[idv]*U1d.flux[j][2] +
				coef_e*fv[idk]*(
					Uv.ge1[ic]*fv[idtx] +
					Uv.ge2[ic]*fv[idty]
				);
    	}
		/* For the same grids(without transformation), the metrics coefficients are:
		 * Uv.gu1 = 0, Uv.gu2 = 1, Uv.gu3 = 1, Uv.guv = 0.;
		 * Uv.gv1 = -2./3, Uv.guv = 0, Uv.gv2 = 0, Uv.gv3 = 4./3;
		 * Uv.ge1 = 0., Uv.ge2 = 1.
		 * which reduced to the origin Navier-Storks Eq.
		 * Confirmed!! */

    	ii = i - config1.Ng;
		for(j=config1.Ng; j<jr; j++)
		{
			jj = j - config1.Ng;
			ic = ii*config1.nj + jj;
			for(k=0; k<neqv; k++)
				rhs[ic][k] += (U1d.flux[j][k] - U1d.flux[j-1][k])/dyc;
		}
	}
}

/*---------------------------------------------------
 * interpolation for dU/dy at (i+1/2, j)
 * ------------------------------------------------*/
void interpoDY(int start, int end)
{
    int i, ir, ii, j, jj, jr, ns, ic, ic1, icm, icp, k;
	double ujr, ujl, vjr, vjl, tjr, tjl;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};

	ir = end + config1.Ng;
    jr = config1.nj + config1.Ng;

    /* 1. interpolation for all the j derivatives. */
	for(i = start + config1.Ng; i<ir; i++)
    {
		for(j=config1.Ng; j<jr; j++)
    	{
    		ic = i*J0 + j;
    		ujr = 0.;
    		vjr = 0.;
    		tjr = 0.;
    		ujl = 0.;
    		vjl = 0.;
    		tjl = 0.;

			// 3. other region
			for(k=0; k<6; k++)
			{
				jj = j - config1.Ng + k;
				icp = i*J0 + jj+1;
				icm = i*J0 + jj;

				ujr = ujr + interpo[k]*U.q[icp][1];
				vjr = vjr + interpo[k]*U.q[icp][2];
				tjr = tjr + interpo[k]*U.tem[icp];
				ujl = ujl + interpo[k]*U.q[icm][1];
				vjl = vjl + interpo[k]*U.q[icm][2];
				tjl = tjl + interpo[k]*U.tem[icm];
			}
			/* 2. get the derivative on (i,j) cells */
			Uv.u_et[ic] = (ujr - ujl)/dyc;
			Uv.v_et[ic] = (vjr - vjl)/dyc;
			Uv.T_et[ic] = (tjr - tjl)/dyc;
		}
    }
}

/*---------------------------------------------------
 * interpolation for dU/dx at (i, j+1/2)
 * ------------------------------------------------*/
void interpoDX(int start, int end)
{
    int i, ir, ii, j, jr, ic, icm, icp;
	double uir, uil, vir, vil, tir, til;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};

	ir = end + config1.Ng;
    jr = config1.nj + config1.Ng;

    /* 1. interpolation for all the i derivatives. */

	for(j=config1.Ng; j<jr; j++)
	{
		for(i=start+config1.Ng; i<ir; i++)
    	{
    		ic = i*J0 + j;
    		uir = 0.;
    		vir = 0.;
    		tir = 0.;
    		uil = 0.;
    		vil = 0.;
    		til = 0.;

			// 3. other region
			for(int k=0; k<6; k++)
			{
				ii = i - 3 + k;
				icp = (ii+1)*J0 + j;
				icm = ii*J0 + j;

				uir = uir + interpo[k]*U.q[icp][1];
				vir = vir + interpo[k]*U.q[icp][2];
				tir = tir + interpo[k]*U.tem[icp];
				uil = uil + interpo[k]*U.q[icm][1];
				vil = vil + interpo[k]*U.q[icm][2];
				til = til + interpo[k]*U.tem[icm];
			}
			/* 2. get the derivative on (i,j) cells */
			Uv.u_xi[ic] = (uir - uil)/dxc;
			Uv.v_xi[ic] = (vir - vil)/dxc;
			Uv.T_xi[ic] = (tir - til)/dxc;
		}
	}
}