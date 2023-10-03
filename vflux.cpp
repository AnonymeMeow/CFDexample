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
void vfluxF(double **rhs)
{
	int i, j, il, ir, jr, ii, jj, k, ik, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, *fv;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void allocatevFlux(int nlen, struct strct_flux *f);
	void freevFlux(int nlen, struct strct_flux *f);
	void interpoDY();

	allocatevFlux(I0,&U1d);
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

	fv = (double*)malloc(sizeof(double)*nvar);
	/*--------------------------------------------------*/

   	il = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

	Ec0 = 1.;
	coef_Re = 1.;
	coef_e  = 1.;

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
			for(k=0; k<neqv; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[i][k] - U1d.flux[i-1][k])/dxc;
		}
	}

    free(fv);
	freevFlux(I0, &U1d);
}


/*--------------------------------------------------------------
 * Calculate viscous Flux in j direction
 * -------------------------------------------------------------*/
void vfluxG(double **rhs)
{
	int i, j, jl, ir, jr, ii, jj, k, jk, ic, ivar, nvar,
		idr, idu, idv, idt, idm, idk, idux, iduy, idvx, idvy, idtx, idty;
	double  Ec0, coef_Re, coef_e, *fv;
	double interpo[6]={1./60., -2./15., 37./60., 37./60., -2./15., 1./60.};
	double approxi[6]={-1./90., 25./180., -245./180., 245./180., -25./180., 1./90.};

	void allocatevFlux(int nlen, struct strct_flux *f);
	void freevFlux(int nlen, struct strct_flux *f);
	void interpoDX();

	allocatevFlux(J0,&U1d);
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

	fv  = (double*)malloc(sizeof(double)*nvar);
	/*--------------------------------------------------*/

   	jl = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
    jr = config1.nj + config1.Ng;

	Ec0 = 1.;
	coef_Re = 1.;
	coef_e  = 1.;

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
			for(k=0; k<neqv; k++)
				rhs[ic][k] = rhs[ic][k] + (U1d.flux[j][k] - U1d.flux[j-1][k])/dyc;
		}
	}

    free(fv);
	freevFlux(J0, &U1d);
}

/*---------------------------------------------------
 * interpolation for dU/dy at (i+1/2, j)
 * ------------------------------------------------*/
void interpoDY()
{
    int i, ir, ii, j, jj, jr, ns, ic, ic1, icm, icp, k;
	double ujr, ujl, vjr, vjl, tjr, tjl;
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
			}
    	}

		/* 2. get the derivative on (i,j) cells */
		Uv.u_et[ic] = (ujr - ujl)/dyc;
		Uv.v_et[ic] = (vjr - vjl)/dyc;
		Uv.T_et[ic] = (tjr - tjl)/dyc;
    }
    	/* 3. assign boundary condition in direction i */
#ifdef MPI_RUN

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
    			ii -= 1;
		    }
		}
#else
		{
			// left side, wall
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

			// right side, solid wall
			ii = ir -1;
			for(i=ir; i<I0; i++)
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
#endif
    }
}

/*---------------------------------------------------
 * interpolation for dU/dx at (i, j+1/2)
 * ------------------------------------------------*/
void interpoDX()
{
    int i, ir, ii, j, jr, ns, ic, icm, icp, k;
	double uir, uil, vir, vil, tir, til;
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
			}
    	}
		/* 2. get the derivative on (i,j) cells */
		Uv.u_xi[ic] = (uir - uil)/dxc;
		Uv.v_xi[ic] = (vir - vil)/dxc;
		Uv.T_xi[ic] = (tir - til)/dxc;
	}
}

/*---------------------------------------------------
 * allocate memory for calculation inviscid flux
 * ------------------------------------------------*/
void allocatevFlux(int nlen, struct strct_flux *f)
{
	int i;

	f->du = (double*)malloc(sizeof(double)*nlen);
	f->dv = (double*)malloc(sizeof(double)*nlen);
	f->dt = (double*)malloc(sizeof(double)*nlen);
	f->rho= (double*)malloc(sizeof(double)*nlen);
	f->u  = (double*)malloc(sizeof(double)*nlen);
	f->v  = (double*)malloc(sizeof(double)*nlen);
	f->p  = (double*)malloc(sizeof(double)*nlen);
	f->t  = (double*)malloc(sizeof(double)*nlen);
	f->mu = (double*)malloc(sizeof(double)*nlen);
	f->kt = (double*)malloc(sizeof(double)*nlen);
	f->flux = (double**)malloc(sizeof(double*)*nlen);
	for(i=0; i<nlen; i++)
	{
		f->flux[i]  = (double*)malloc(sizeof(double)*neqv);
	}
}
/*---------------------------------------------------
 * free the memory
 * ------------------------------------------------*/
void freevFlux(int nlen, struct strct_flux *f)
{
	int i;

	free(f->du);
	free(f->dv);
	free(f->dt);
	free(f->rho);
	free(f->u);
	free(f->v);
	free(f->p);
	free(f->t);
	free(f->mu);
	free(f->kt);

	for(i=0; i<nlen; i++)
	{
		free(f->flux[i]);
	}
	free(f->flux);
}