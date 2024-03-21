/*
 * vflux.c
 *
 *  Created on:    Apr 4, 2014
 *  last modified: Sep 10, 2014
 */

#include<string.h>
#include<math.h>
#include"comm.hpp"
#include"chemdata.hpp"

#define NSPEC_BUFFER 10

/*--------------------------------------------------------------
 * Calculate viscous flux in x direction
 * -------------------------------------------------------------*/
void vfluxF(double **rhs)
{
	int i, j, il, ir, jr, ii, jj, k, ik, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, fv[12];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void interpoDY();

	interpoDY();

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

   	il = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    if(config1.nonDi)
    {
    	Ec0 = uRef*uRef/(config2.gam0*cvRef*temRef);
    	coef_Re = 1./config2.Re0;
    	coef_e  = 1./(config2.Pr0*config2.Re0*Ec0);
    }
    else
    {
    	Ec0 = 1.;
        coef_Re = 1.;
        coef_e  = 1.;
    }

    for(j=config1.Ng; j<jr; j++)
    {
        for(i=0; i<I0; i++)
        {
        	/*---- convert to 1D-array ----*/
        	ic = i*J0 + j;

        	U1d.du[i] = Uv.u_xi[ic];
            U1d.dv[i] = Uv.v_xi[ic];
        	U1d.dt[i] = Uv.T_xi[ic];
        	U1d.rho[i]=  Ug.q[ic][0];
        	U1d.u[i]  =  Ug.q[ic][1];
        	U1d.v[i]  =  Ug.q[ic][2];
        	U1d.p[i]  =  Ug.pre[ic];
        	U1d.t[i]  =  Ug.tem[ic];
        	U1d.mu[i] =  Ug.mu[ic];
        	U1d.kt[i] =  Ug.kt[ic];
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
    		U1d.flux[i][1] = coef_Re*fv[idm]*(Uv.fu1[ic]*fv[idux] + Uv.fu2[ic]*fv[iduy]
    		               + Uv.fuv[ic]*fv[idvx]  + Uv.fu3[ic]*fv[idvy]);
    		U1d.flux[i][2] = coef_Re*fv[idm]*(Uv.fuv[ic]*fv[idux] + Uv.fv1[ic]*fv[iduy]
    		               + Uv.fv2[ic]*fv[idvx]  + Uv.fv3[ic]*fv[idvy]);
    		U1d.flux[i][3] = fv[idu]*U1d.flux[i][1] + fv[idv]*U1d.flux[i][2]
    		               + coef_e*fv[idk]*(Uv.fe1[ic]*fv[idtx] + Uv.fe2[ic]*fv[idty]);
    	}
		/* For the same grids(without transformation), the metrics coefficients are:
		 * Uv.fu1 = 4./3, Uv.fu2 = 0, Uv.fu3 = -2./3, Uv.fuv = 0.;
		 * Uv.fv1 = 1,  Uv.fuv=0, Uv.fv2 = 1, Uv.fv3 = 0.;
		 * Uv.fe1 = 1., Uv.ge2 = 0.
		 * which reduced to the origin Navier-Storks Eq.
		 * Confirmed! */

		jj = j-config1.Ng;
		for(i=config1.Ng; i<ir; i++)
		{
			ii = i - config1.Ng;
			ic = ii*config1.nj + jj;
			for(k=0; k<neqn; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[i][k] - U1d.flux[i-1][k])/dxc;
		}
	}
}


/*--------------------------------------------------------------
 * Calculate viscous Flux in j direction
 * -------------------------------------------------------------*/
void vfluxG(double **rhs)
{
	int i, j, jl, ir, jr, ii, jj, k, jk, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, fv[12];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void interpoDX();

	interpoDX();

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
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    if(config1.nonDi)
    {   // non-dimensional factor
    	Ec0 = uRef*uRef/(config2.gam0*cvRef*temRef);
    	coef_Re = 1./config2.Re0;
    	coef_e  = 1./(config2.Pr0*config2.Re0*Ec0);
    }
    else
    {
    	Ec0 = 1.;
        coef_Re = 1.;
        coef_e  = 1.;
    }

    for(i=config1.Ng; i<ir; i++)
    {
        for(j=0; j<J0; j++)
        {
        	/*---- convert to 1D-array ----*/
        	ic = i*J0 + j;

        	U1d.du[j] = Uv.u_xi[ic];
            U1d.dv[j] = Uv.v_xi[ic];
        	U1d.dt[j] = Uv.T_xi[ic];
        	U1d.rho[j]=  Ug.q[ic][0];
        	U1d.u[j]  =  Ug.q[ic][1];
        	U1d.v[j]  =  Ug.q[ic][2];
        	U1d.p[j]  =  Ug.pre[ic];
        	U1d.t[j]  =  Ug.tem[ic];
        	U1d.mu[j] =  Ug.mu[ic];
        	U1d.kt[j] =  Ug.kt[ic];
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
    			 U1d.flux[j][1] = coef_Re*fv[idm]*(Uv.gu1[ic]*fv[idux] + Uv.gu2[ic]*fv[iduy]
    		                           + Uv.gu3[ic]*fv[idvx]  + Uv.guv[ic]*fv[idvy]);
    			 U1d.flux[j][2] = coef_Re*fv[idm]*(Uv.gv1[ic]*fv[idux] + Uv.guv[ic]*fv[iduy]
    		                           + Uv.gv2[ic]*fv[idvx]  + Uv.gv3[ic]*fv[idvy]);
    			 U1d.flux[j][3] = fv[idu]*U1d.flux[j][1] + fv[idv]*U1d.flux[j][2]
    			                + coef_e*fv[idk]*(Uv.ge1[ic]*fv[idtx] + Uv.ge2[ic]*fv[idty]);
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
			for(k=0; k<neqn; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[j][k] - U1d.flux[j-1][k])/dyc;
		}
	}
}

/*--------------------------------------------------------------
 * Calculate viscous Flux in x direction of real-gas flow
 * -------------------------------------------------------------*/
void vfluxchemF(double **rhs)
{
	int i, j, il, ir, jr, ii, jj, k, ik, ic, ns, ivar, nvar, m,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double fvhs[NSPEC_BUFFER], fvdiff[NSPEC_BUFFER], gradx[NSPEC_BUFFER], grady[NSPEC_BUFFER], fv[NSPEC_BUFFER + 12];
	double hs, es, cpref, sumqs, Ec0, coef_rhos, coef_hs, coef_e, coef_Re;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void interpoDY();
	double getes(int ns, double t);

	interpoDY();

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

    m = config1.nspec;
    nvar = 12 + m;

	/*--------------------------------------------------*/

   	il = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    if(config1.nonDi)
    {
		cpref = config2.gam0*cvRef;
    	Ec0   = uRef*uRef/(cpref*temRef);
    	coef_Re = 1./config2.Re0;
		coef_hs = (rhoRef*diffRef*cpref/condRef)*Ec0;
    	coef_e = 1./(config2.Pr0*config2.Re0*Ec0);
    	coef_rhos = diffRef/(LRef*uRef);
    }
    else
    {
    	coef_Re = 1.;
        coef_hs = 1.;
        coef_e  = 1.;
        coef_rhos = 1.;
    }

    for(j=config1.Ng; j<jr; j++)
    {
		for(i=0; i<I0; i++)
		{
			/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.du[i] = Uv.u_et[ic];
			U1d.dv[i] = Uv.v_et[ic];
			U1d.dt[i] = Uv.T_et[ic];
			U1d.rho[i]=  Ug.q[ic][0];
			U1d.u[i]  =  Ug.q[ic][1];
			U1d.v[i]  =  Ug.q[ic][2];
			U1d.p[i]  =  Ug.pre[ic];
			U1d.t[i]  =  Ug.tem[ic];
			U1d.mu[i] =  Ug.mu[ic];
			U1d.kt[i] =  Ug.kt[ic];
			for(ns=0; ns<m; ns++)
			{
				U1d.qs[i][ns] = Ug.qs[ic][ns];
				U1d.Ds[i][ns] = Ug.di[ic][ns];
				U1d.dqs[i][ns] = Uv.qs_et[ic][ns];
			}
		}

		 /*---------- 1. get the face value ----------*/
    	for(i=il; i<ir; i++)
    	{
    		ic = i*J0 + j;

    		for(ivar=0; ivar<nvar; ivar++)
    			fv[ivar] = 0.0;

			for(ns=0; ns<m; ns++)
			{
				fvhs[ns]   = 0.;
				gradx[ns]  = 0.;
				grady[ns]  = 0.;
				fvdiff[ns] = 0.;
			}
	    		for(k=0; k<6; k++)
	    		{
	    			ik  = i - 2 + k;
    	/*---------- 1.3 high-order interpolate i derivatives at the i+1/2 faces ----------*/
	    			fv[idr]  = fv[idr]  + interpo[k]*U1d.rho[ik];
	    			fv[idu]  = fv[idu]  + interpo[k]*U1d.u[ik];
	    			fv[idv]  = fv[idv]  + interpo[k]*U1d.v[ik];
	    			fv[idt]  = fv[idt]  + interpo[k]*U1d.t[ik];
	    			fv[idm]  = fv[idm]  + interpo[k]*U1d.mu[ik];
	    			fv[idk]  = fv[idk]  + interpo[k]*U1d.kt[ik];
	    			fv[iduy] = fv[iduy] + interpo[k]*U1d.du[ik];
	    			fv[idvy] = fv[idvy] + interpo[k]*U1d.dv[ik];
	    			fv[idty] = fv[idty] + interpo[k]*U1d.dt[ik];

	    			for(ns=0; ns<m; ns++)
	    			{
	    				es = getes(ns, U1d.t[ik]);
	    				hs = es + (ru/specData[ns].wm/rgasRef)*U1d.t[ik]*Upsilon;
	    				fvhs[ns]   = fvhs[ns]   + interpo[k]*hs;
	    				grady[ns]  = grady[ns]  + interpo[k]*U1d.dqs[ik][ns];
	    				fvdiff[ns] = fvdiff[ns] + interpo[k]*U1d.Ds[ik][ns];
	    			}
        /*---------- 1.4 high-order approximate j derivatives at i+1/2 faces ----------*/
	    			fv[idux] = fv[idux] + approxi[k]*U1d.u[ik]/dxc;
	    			fv[idvx] = fv[idvx] + approxi[k]*U1d.v[ik]/dxc;
	    			fv[idtx] = fv[idtx] + approxi[k]*U1d.t[ik]/dxc;
	    			for(ns=0; ns<m; ns++)
	    				gradx[ns] = gradx[ns] + approxi[k]*U1d.qs[ik][ns]/dxc;
	    		}

        /*---------- 2.calculate the viscous Flux ----------*/
    		sumqs = 0.;
    		for(ns=0; ns<m; ns++)
    		{
    			U1d.flux[i][ns] = coef_rhos*fv[idr]*fvdiff[ns]
    			            *(Uv.fe1[ic]*gradx[ns] + Uv.fe2[ic]*grady[ns]);
				sumqs = sumqs + U1d.flux[i][ns];
    		}
    		for(ns=0; ns<m; ns++)
    		{
    			// correction of the diffusive flux
    			U1d.flux[i][ns] = U1d.flux[i][ns] - U1d.qs[i][ns]*sumqs;
    		}

    		U1d.flux[i][m]   = coef_Re*fv[idm]*(Uv.fu1[ic]*fv[idux] + Uv.fu2[ic]*fv[iduy]
    		                              + Uv.fuv[ic]*fv[idvx] + Uv.fu3[ic]*fv[idvy]);
    		U1d.flux[i][m+1] = coef_Re*fv[idm]*(Uv.fuv[ic]*fv[idux] + Uv.fv1[ic]*fv[iduy]
    		                              + Uv.fv2[ic]*fv[idvx] + Uv.fv3[ic]*fv[idvy]);
    		U1d.flux[i][m+2] = fv[idu]*U1d.flux[i][m] + fv[idv]*U1d.flux[i][m+1]
    		             + coef_e*fv[idk]*(Uv.fe1[ic]*fv[idtx] + Uv.fe2[ic]*fv[idty]);

    		sumqs = 0.;
    		for(ns=0; ns<m; ns++)
    		{
    			sumqs = sumqs + coef_hs*fvhs[ns]*fvdiff[ns]
    			              * (Uv.fe1[ic]*gradx[ns] + Uv.fe2[ic])*grady[ns];
    		}
    		U1d.flux[i][m+2] = U1d.flux[i][m+2] + coef_e*fv[idr]*sumqs;
    	}

		jj = j-config1.Ng;
		for(i=config1.Ng; i<ir; i++)
		{
			ii = i - config1.Ng;
			ic = ii*config1.nj + jj;
			for(k=0; k<neqn; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[i][k] - U1d.flux[i-1][k])/dxc;
		}
	}
}


/*--------------------------------------------------------------
 * Calculate viscous Flux in y direction of real-gas flow
 * -------------------------------------------------------------*/
void vfluxchemG(double **rhs)
{
	int i, j, jl, ir, jr, ii, jj, k, jk, ic, ns, ivar, m, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double fvhs[NSPEC_BUFFER], fvdiff[NSPEC_BUFFER], gradx[NSPEC_BUFFER], grady[NSPEC_BUFFER], fv[NSPEC_BUFFER + 12];
	double hs, es, sumqs, Ec0, cpref, coef_rhos, coef_hs, coef_e, coef_Re;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void interpoDX();
	double getes(int ns, double t);

	interpoDX();

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

    m = config1.nspec;
	nvar = 12 + m;
	/*--------------------------------------------------*/

   	jl = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    if(config1.nonDi)
    {
		cpref = config2.gam0*cvRef;
    	Ec0   = uRef*uRef/(cpref*temRef);
    	coef_Re = 1./config2.Re0;
		coef_hs = (rhoRef*diffRef*cpref/condRef)*Ec0;
    	coef_e  = 1./(config2.Pr0*config2.Re0*Ec0);
    	coef_rhos = diffRef/(LRef*uRef);
    }
    else
    {
    	coef_Re = 1.;
        coef_hs = 1.;
        coef_e  = 1.;
        coef_rhos = 1.;
    }

    for(i=config1.Ng; i<ir; i++)
    {
		for(j=0; j<J0; j++)
		{
			/*---- convert to 1D-array ----*/
			ic = i*J0 + j;

			U1d.du[j] = Uv.u_xi[ic];
			U1d.dv[j] = Uv.v_xi[ic];
			U1d.dt[j] = Uv.T_xi[ic];
			U1d.rho[j]=  Ug.q[ic][0];
			U1d.u[j]  =  Ug.q[ic][1];
			U1d.v[j]  =  Ug.q[ic][2];
			U1d.p[j]  =  Ug.pre[ic];
			U1d.t[j]  =  Ug.tem[ic];
			U1d.mu[j] =  Ug.mu[ic];
			U1d.kt[j] =  Ug.kt[ic];
			for(ns=0; ns<m; ns++)
			{
				U1d.qs[j][ns] = Ug.qs[ic][ns];
				U1d.Ds[j][ns] = Ug.di[ic][ns];
				U1d.dqs[j][ns] = Uv.qs_xi[ic][ns];
			}
		}

        /*---------- 1. get the face value ----------*/
    	for(j=jl; j<jr; j++)
    	{
    		ic = i*J0 + j;

    		for(ivar=0; ivar<nvar; ivar++)
    			fv[ivar] = 0.0;

    		for(ns=0; ns<m; ns++)
    		{
    			fvhs[ns]   = 0.;
    			gradx[ns]  = 0.;
    			grady[ns]  = 0.;
    			fvdiff[ns] = 0.;
    		}

    		for(k=0; k<6; k++)
    		{
    			jk  = j - 2 + k;
    	/*---------- 1.3 high-order interpolate i derivatives at the j+1/2 faces ----------*/
    			fv[idr]  = fv[idr]  + interpo[k]*U1d.rho[jk];
    			fv[idu]  = fv[idu]  + interpo[k]*U1d.u[jk];
    			fv[idv]  = fv[idv]  + interpo[k]*U1d.v[jk];
    			fv[idt]  = fv[idt]  + interpo[k]*U1d.t[jk];
    			fv[idm]  = fv[idm]  + interpo[k]*U1d.mu[jk];
    			fv[idk]  = fv[idk]  + interpo[k]*U1d.kt[jk];
    			fv[idux] = fv[idux] + interpo[k]*U1d.du[jk];
    			fv[idvx] = fv[idvx] + interpo[k]*U1d.dv[jk];
    			fv[idtx] = fv[idtx] + interpo[k]*U1d.dt[jk];

    			for(ns=0; ns<m; ns++)
    			{
    				es = getes(ns, U1d.t[jk]);
    				hs = es + (ru/specData[ns].wm/rgasRef)*U1d.t[jk]*Upsilon;
    				fvhs[ns]   = fvhs[ns]   + interpo[k]*hs;
    				gradx[ns]  = gradx[ns]  + interpo[k]*U1d.dqs[jk][ns];
    				fvdiff[ns] = fvdiff[ns] + interpo[k]*U1d.Ds[jk][ns];
    			}
        /*---------- 1.4 high-order approximate j derivatives at j+1/2 faces ----------*/
    			fv[iduy] = fv[iduy] + approxi[k]*U1d.u[jk]/dyc;
    			fv[idvy] = fv[idvy] + approxi[k]*U1d.v[jk]/dyc;
    			fv[idty] = fv[idty] + approxi[k]*U1d.t[jk]/dyc;
    			for(ns=0; ns<m; ns++)
    				grady[ns] = grady[ns] + approxi[k]*U1d.qs[jk][ns]/dyc;
    		}

    		/*---------- 2.calculate the viscous Flux ----------*/
    		sumqs = 0.;
    		for(ns=0; ns<m; ns++)
    		{
    			U1d.flux[j][ns] = coef_rhos*fv[idr]*fvdiff[ns]
    			            * (Uv.ge1[ic]*gradx[ns] + Uv.ge2[ic]*grady[ns]);
				sumqs = sumqs + U1d.flux[j][ns];
    		}

    		for(ns=0; ns<m; ns++)
    		{
    			// correction of the diffusive flux
    			U1d.flux[j][ns] = U1d.flux[j][ns] - U1d.qs[j][ns]*sumqs;
    		}

    		U1d.flux[j][m]   = coef_Re*fv[idm]*(Uv.gu1[ic]*fv[idux] + Uv.gu2[ic]*fv[iduy]
    		                              + Uv.gu3[ic]*fv[idvx] + Uv.guv[ic]*fv[idvy]);
    		U1d.flux[j][m+1] = coef_Re*fv[idm]*(Uv.gv1[ic]*fv[idux] + Uv.guv[ic]*fv[iduy]
    		                              + Uv.gv2[ic]*fv[idvx] + Uv.gv3[ic]*fv[idvy]);
    		U1d.flux[j][m+2] = fv[idu]*U1d.flux[j][m] + fv[idv]*U1d.flux[j][m+1]
    		                 + coef_e*fv[idk]*(Uv.ge1[ic]*fv[idtx] + Uv.ge2[ic]*fv[idty]);

    		sumqs = 0.;
    		for(ns=0; ns<config1.nspec; ns++)
    		{
    			sumqs = sumqs + coef_hs*fvhs[ns]*fvdiff[ns]
    			              * (Uv.ge1[ic]*gradx[ns] + Uv.ge2[ic]*grady[ns]);
    		}
    		U1d.flux[j][m+2] = U1d.flux[j][m+2] + coef_e*fv[idr]*sumqs;
    	}

    	ii = i - config1.Ng;
		for(j=config1.Ng; j<jr; j++)
		{
			jj = j - config1.Ng;
			ic = ii*config1.nj + jj;
			for(k=0; k<neqn; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[j][k] - U1d.flux[j-1][k])/dyc;
		}
	}
}

/*---------------------------------------------------
 * interpolation for dU/dy at (i+1/2, j)
 * ------------------------------------------------*/
void interpoDY()
{
    int i, ir, ii, j, jj, jr, ns, ic, ic1, icm, icp, k;
	double ujr, ujl, vjr, vjl, tjr, tjl, ysjr[NSPEC_BUFFER], ysjl[NSPEC_BUFFER];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};

	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    /* 1. interpolation for all the j derivatives. */
	for(i=config1.Ng; i<ir; i++)
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
    		for(ns=0; ns<config1.nspec; ns++)
    		{
    			ysjr[ns] = 0.;
    			ysjl[ns] = 0.;
    		}

			// 3. other region
			for(k=0; k<6; k++)
			{
				jj = j - config1.Ng + k;
				icp = i*J0 + jj+1;
				icm = i*J0 + jj;

				ujr = ujr + interpo[k]*Ug.q[icp][1];
				vjr = vjr + interpo[k]*Ug.q[icp][2];
				tjr = tjr + interpo[k]*Ug.tem[icp];
				ujl = ujl + interpo[k]*Ug.q[icm][1];
				vjl = vjl + interpo[k]*Ug.q[icm][2];
				tjl = tjl + interpo[k]*Ug.tem[icm];
				if(config1.gasModel != 0)
					for(ns=0; ns<config1.nspec; ns++)
					{
						ysjr[ns] = ysjr[ns] + interpo[k]*Ug.qs[icp][ns];
						ysjl[ns] = ysjl[ns] + interpo[k]*Ug.qs[icm][ns];
					}
			}
    		/* 2. get the derivative on (i,j) cells */
    		Uv.u_et[ic] = (ujr - ujl)/dyc;
    		Uv.v_et[ic] = (vjr - vjl)/dyc;
    		Uv.T_et[ic] = (tjr - tjl)/dyc;
			if(config1.gasModel != 0)
				for(ns=0; ns<config1.nspec; ns++)
					Uv.qs_et[ic][ns] = (ysjr[ns] - ysjl[ns])/dyc;
    	}
    }

	/* 3. assign boundary condition in direction i */
	for(j=config1.Ng; j<jr; j++)
	{
    	/*---- assign MPI boundary ----*/
    	/* left side, wall */
    	if(MyID == 0)
    	{
            ii = 2*config1.Ng - 1;
        	for(i=0; i<config1.Ng; i++)
		    {
    			ic = i*J0 + j;
    	    	ic1 = ii*J0 + j;
		    	Uv.u_et[ic] = -Uv.u_et[ic1];
		    	Uv.v_et[ic] = -Uv.v_et[ic1];
		    	Uv.T_et[ic] =  Uv.T_et[ic1];
				if(config1.gasModel != 0)
					for(ns=0; ns<config1.nspec; ns++)
						Uv.qs_et[ic][ns] = Uv.qs_et[ic1][ns];
			    ii -= 1;
		    }
    	}
    	else
    	{
    		/* MPI boundary */
    		jj = j - config1.Ng;
    		for(i=0; i<config1.Ng; i++)
		    {
    			ic  = i*J0 + j;
    			if(jj == config1.nj-1)
    			{
    				icm = i*config1.nj + jj-1;
    				icp = i*config1.nj + jj;
    			}
    			else
    			{
    				icm = i*config1.nj + jj;
    				icp = i*config1.nj + jj+1;
    			}
		    	Uv.u_et[ic] =  (mpiRecv_ql[icp].u - mpiRecv_ql[icm].u)/dyc;
		    	Uv.v_et[ic] =  (mpiRecv_ql[icp].v - mpiRecv_ql[icm].v)/dyc;
		    	Uv.T_et[ic] =  (mpiRecv_ql[icp].t - mpiRecv_ql[icm].t)/dyc;
				if(config1.gasModel != 0)
					for(ns=0; ns<config1.nspec; ns++)
						Uv.qs_et[ic][ns] = (mpiRecv_ql[icp].qs[ns] - mpiRecv_ql[icm].qs[ns])/dyc;
		    }
    	}
    	/* right side, solid Wall */
    	if(MyID == NMAXproc)
    	{
        	ii = ir -1;
        	for(i=ir; i<I0; i++)
        	{
    			ic  = i*J0  + j;
     		    ic1 = ii*J0 + j;
		    	Uv.u_et[ic] = -Uv.u_et[ic1];
		    	Uv.v_et[ic] = -Uv.v_et[ic1];
		    	Uv.T_et[ic] =  Uv.T_et[ic1];
				if(config1.gasModel != 0)
					for(ns=0; ns<config1.nspec; ns++)
						Uv.qs_et[ic][ns] = Uv.qs_et[ic1][ns];
			    ii -= 1;
		    }
    	}
    	else
    	{
    		/* MPI boundary */
    		jj = j - config1.Ng;
    		for(i=0; i<config1.Ng; i++)
    		{
    			ii = ir + i;

    			ic = ii*J0 + j;

    			if(jj == config1.nj-1)
    			{
    				icm = i*config1.nj + jj-1;
    				icp = i*config1.nj + jj;
    			}
    			else
    			{
    				icm = i*config1.nj + jj;
    				icp = i*config1.nj + jj+1;
    			}
    			Uv.u_et[ic] =  (mpiRecv_qr[icp].u - mpiRecv_qr[icm].u)/dyc;
    			Uv.v_et[ic] =  (mpiRecv_qr[icp].v - mpiRecv_qr[icm].v)/dyc;
    			Uv.T_et[ic] =  (mpiRecv_qr[icp].t - mpiRecv_qr[icm].t)/dyc;
    			if(config1.gasModel != 0)
    				for(ns=0; ns<config1.nspec; ns++)
    					Uv.qs_et[ic][ns] = (mpiRecv_qr[icp].qs[ns] - mpiRecv_qr[icm].qs[ns])/dyc;
    			ii -= 1;
		    }
		}
    }
}

/*---------------------------------------------------
 * interpolation for dU/dx at (i, j+1/2)
 * ------------------------------------------------*/
void interpoDX()
{
    int i, ir, ii, j, jr, ns, ic, icm, icp, k;
	double uir, uil, vir, vil, tir, til, ysir[NSPEC_BUFFER], ysil[NSPEC_BUFFER];
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};

	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

    /* 1. interpolation for all the i derivatives. */

	for(j=config1.Ng; j<jr; j++)
	{
		for(i=config1.Ng; i<ir; i++)
    	{
    		ic = i*J0 + j;
    		uir = 0.;
    		vir = 0.;
    		tir = 0.;
    		uil = 0.;
    		vil = 0.;
    		til = 0.;
    		for(ns=0; ns<config1.nspec; ns++)
    		{
    			ysir[ns] = 0.;
    			ysil[ns] = 0.;
    		}

			// 3. other region
			for(k=0; k<6; k++)
			{
				ii = i - 3 + k;
				icp = (ii+1)*J0 + j;
				icm = ii*J0 + j;

				uir = uir + interpo[k]*Ug.q[icp][1];
				vir = vir + interpo[k]*Ug.q[icp][2];
				tir = tir + interpo[k]*Ug.tem[icp];
				uil = uil + interpo[k]*Ug.q[icm][1];
				vil = vil + interpo[k]*Ug.q[icm][2];
				til = til + interpo[k]*Ug.tem[icm];
				if(config1.gasModel != 0)
					for(ns=0; ns<config1.nspec; ns++)
					{
						ysir[ns] = ysir[ns] + interpo[k]*Ug.qs[icp][ns];
						ysil[ns] = ysil[ns] + interpo[k]*Ug.qs[icm][ns];
					}
			}
			/* 2. get the derivative on (i,j) cells */
			Uv.u_xi[ic] = (uir - uil)/dxc;
			Uv.v_xi[ic] = (vir - vil)/dxc;
			Uv.T_xi[ic] = (tir - til)/dxc;
			if(config1.gasModel != 0)
				for(ns=0; ns<config1.nspec; ns++)
					Uv.qs_xi[ic][ns] = (ysir[ns] - ysil[ns])/dxc;
    	}
	}
}