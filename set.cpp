/*
 * set.c
 *
 *  Created on: 2014.01.19
 *
 */

#include"comm.h"

/*---------------------------------------------------
 * set and initialize the job
 * ------------------------------------------------*/
void setjob()
{
	void initjob();
	void setGeom();

	void allocateU(int nlen, struct strct_U *U);
	void allocateUv();
	void allocateU1d(int);
	void allocateOthers();

	int nc  = config1.ni*config1.nj;
	int nc1 = I0*J0;

	allocateU(nc, &U);
	allocateU(nc1, &Ug);
	allocateUv();
	allocateU1d(MAX(I0,J0));
	allocateOthers();

	initjob();

	setGeom();
}

/*---------------------------------------------------
 * Conduct the initialization for the simulation
 * ------------------------------------------------*/
void initjob()
{
	int ic, icp, x_sh;
	void assigncells(int i1, int in, int j1, int jn, double u,
			         double v, double t, double p);

	config2.t0 = 0;
	config1.iStep0 = 1;

	double dis = config2.x0*config2.Lx;
	for(x_sh=0; x_sh<config1.ni; x_sh++)
	{
		ic = x_sh*J0 + 0;
		icp = (x_sh+1)*J0 + 0;
		if((mesh.x[ic]<=dis) && (dis<mesh.x[icp]))
			break;
	}

	assigncells(0,x_sh, 0,config1.nj, config2.u1,config2.v1,config2.T1,config2.p1);
	assigncells(x_sh,config1.ni, 0,config1.nj, config2.u2,config2.v2,config2.T2,config2.p2);

	/*
	MPI_Barrier(MPI_COMM_WORLD); // Lock Here!
	*/
}
/*-----------------------------------------------------------
 * Assign initial value to each cells
 * ---------------------------------------------------------*/
void assigncells(int i1, int in, int j1, int jn, double u,
		         double v, double t, double p)
{
	double rgas1 = (ru/config2.molWeight);

	for(int i=i1; i<in; i++)
	{
		for(int j=j1; j<jn; j++)
		{
			int ic = i*config1.nj + j;

			U.q[ic][1] = u;
			U.q[ic][2] = v;
			U.tem[ic]  = t;
			U.pre[ic]  = p;

			double ek = 0.5*(u*u + v*v);

			double RT = rgas1*t;
			U.q[ic][0] = p/RT;
			U.q[ic][3] = ek + RT/(config2.gam0-1);
		}
	}
}

/*---------------------------------------------------
 * Set the geometry related variables
 * ------------------------------------------------*/
void setGeom()
{
	int i, j, il, jl, ir, jr, ic, icp;
	double  xxi, xet, yxi, yet, xix, xiy, etx, ety, sxsx, sxsy, sysy,
	        exex, exey, sxex, sxey, eyey, syex, syey, yas;

	il = config1.Ng - 1;
	jl = config1.Ng - 1;
	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

	for(j=config1.Ng; j<jr; j++)
	{
		for(i=il; i<ir; i++)
		{
			ic  = i*J0 + j;
			icp = (i+1)*J0 + j;

			/* the metrics for viscous term */
			xxi = 0.5*(mesh.x_xi[ic] + mesh.x_xi[icp]);
			yxi = 0.5*(mesh.y_xi[ic] + mesh.y_xi[icp]);
			xet = 0.5*(mesh.x_et[ic] + mesh.x_et[icp]);
			yet = 0.5*(mesh.y_et[ic] + mesh.y_et[icp]);
			yas = 0.5*(mesh.yaks[ic] + mesh.yaks[icp]);

			xix =  yet/yas;
			xiy = -xet/yas;
			etx = -yxi/yas;
			ety =  xxi/yas;

			sxsx = xix * xix;
			sysy = xiy * xiy;
			sxsy = xix * xiy;

			sxex = xix * etx;
			sxey = xix * ety;
			syex = xiy * etx;
			syey = xiy * ety;

			Uv.fu1[ic] = (sysy + 4./3.*sxsx) * yas;
			Uv.fu2[ic] = (syey + 4./3.*sxex) * yas;
			Uv.fu3[ic] = (syex - 2./3.*sxey) * yas;
			Uv.fv1[ic] = (sxey - 2./3.*syex) * yas;
			Uv.fv2[ic] = (sxsx + 4./3.*sysy) * yas;
			Uv.fv3[ic] = (sxex + 4./3.*syey) * yas;

			Uv.fuv[ic] = (sxsy /3.) * yas;

			Uv.fe1[ic] = (sxsx + sysy) * yas;
			Uv.fe2[ic] = (sxex + syey) * yas;
		}
	}

	/*---------- j direction ----------*/
	for(i=config1.Ng; i<ir; i++)
	{
		for(j=jl; j<jr; j++)
		{
			ic  = i*J0 + j;
			icp = i*J0 + j+1;

		/*---------- 3. average the values ----------*/
			xxi = 0.5*(mesh.x_xi[ic] + mesh.x_xi[icp]);
			yxi = 0.5*(mesh.y_xi[ic] + mesh.y_xi[icp]);
			xet = 0.5*(mesh.x_et[ic] + mesh.x_et[icp]);
			yet = 0.5*(mesh.y_et[ic] + mesh.y_et[icp]);
			yas = 0.5*(mesh.yaks[ic] + mesh.yaks[icp]);

		/*---------- 4. calculate geometry factor ----------*/
			xix =  yet/yas;
			xiy = -xet/yas;
			etx = -yxi/yas;
			ety =  xxi/yas;

			exex = etx * etx;
			eyey = ety * ety;
			exey = etx * ety;

			sxex = xix * etx;
			sxey = xix * ety;
			syex = xiy * etx;
			syey = xiy * ety;


			Uv.gu1[ic] = (syey + 4./3.*sxex) * yas;
			Uv.gu2[ic] = (eyey + 4./3.*exex) * yas;
			Uv.gu3[ic] = (sxey - 2./3.*syex) * yas;
			Uv.gv1[ic] = (syex - 2./3.*sxey) * yas;
			Uv.gv2[ic] = (sxex + 4./3.*syey) * yas;
			Uv.gv3[ic] = (exex + 4./3.*eyey) * yas;
			Uv.guv[ic] = (exey /3.) * yas;

			Uv.ge1[ic]  = (sxex + syey) * yas;
			Uv.ge2[ic]  = (exex + eyey) * yas;
		}
	}
}
/*---------------------------------------------------
 * allocate memory for U vector
 * ------------------------------------------------*/
void allocateU(int nlen, struct strct_U *U)
{
	U->mu   = (double*)malloc(sizeof(double)*nlen);
	U->kt   = (double*)malloc(sizeof(double)*nlen);
	U->cv   = (double*)malloc(sizeof(double)*nlen);
	U->rgas = (double*)malloc(sizeof(double)*nlen);
	U->pre  = (double*)malloc(sizeof(double)*nlen);
	U->tem  = (double*)malloc(sizeof(double)*nlen);
	U->gam  = (double*)malloc(sizeof(double)*nlen);
	U->q    = (double**)malloc(sizeof(double*)*nlen);

	for(int i=0; i<nlen; i++)
		U->q[i]  = (double*)malloc(sizeof(double)*neqv);
}

/*---------------------------------------------------
 * allocate memory for derivatives
 * ------------------------------------------------*/
void allocateUv()
{
	int nlen = I0*J0;

	Uv.fu1 = (double*)malloc(sizeof(double)*nlen);
	Uv.fu2 = (double*)malloc(sizeof(double)*nlen);
	Uv.fu3 = (double*)malloc(sizeof(double)*nlen);
	Uv.fv1 = (double*)malloc(sizeof(double)*nlen);
	Uv.fv2 = (double*)malloc(sizeof(double)*nlen);
	Uv.fv3 = (double*)malloc(sizeof(double)*nlen);
	Uv.fuv = (double*)malloc(sizeof(double)*nlen);
	Uv.fe1 = (double*)malloc(sizeof(double)*nlen);
	Uv.fe2 = (double*)malloc(sizeof(double)*nlen);

	Uv.gu1 = (double*)malloc(sizeof(double)*nlen);
	Uv.gu2 = (double*)malloc(sizeof(double)*nlen);
	Uv.gu3 = (double*)malloc(sizeof(double)*nlen);
	Uv.gv1 = (double*)malloc(sizeof(double)*nlen);
	Uv.gv2 = (double*)malloc(sizeof(double)*nlen);
	Uv.gv3 = (double*)malloc(sizeof(double)*nlen);
	Uv.guv = (double*)malloc(sizeof(double)*nlen);
	Uv.ge1 = (double*)malloc(sizeof(double)*nlen);
	Uv.ge2 = (double*)malloc(sizeof(double)*nlen);

	// 2. allocate U derivatives
	Uv.u_xi = (double*)malloc(sizeof(double)*nlen);
	Uv.u_et = (double*)malloc(sizeof(double)*nlen);
	Uv.v_xi = (double*)malloc(sizeof(double)*nlen);
	Uv.v_et = (double*)malloc(sizeof(double)*nlen);
	Uv.T_xi = (double*)malloc(sizeof(double)*nlen);
	Uv.T_et = (double*)malloc(sizeof(double)*nlen);
	for (int i = 0; i < nlen; i++)
	{
		Uv.u_xi[i] = 0;
		Uv.u_et[i] = 0;
		Uv.v_xi[i] = 0;
		Uv.v_et[i] = 0;
		Uv.T_xi[i] = 0;
		Uv.T_et[i] = 0;
	}
}

void allocateU1d(int nlen)
{
	U1d.xix = (double*)malloc(sizeof(double)*nlen);
	U1d.xiy = (double*)malloc(sizeof(double)*nlen);
	U1d.etx = (double*)malloc(sizeof(double)*nlen);
	U1d.ety = (double*)malloc(sizeof(double)*nlen);
	U1d.yas = (double*)malloc(sizeof(double)*nlen);
	U1d.rho = (double*)malloc(sizeof(double)*nlen);
	U1d.du = (double*)malloc(sizeof(double)*nlen);
	U1d.dv = (double*)malloc(sizeof(double)*nlen);
	U1d.dt = (double*)malloc(sizeof(double)*nlen);
	U1d.u = (double*)malloc(sizeof(double)*nlen);
	U1d.v = (double*)malloc(sizeof(double)*nlen);
	U1d.p = (double*)malloc(sizeof(double)*nlen);
	U1d.e = (double*)malloc(sizeof(double)*nlen);
	U1d.t = (double*)malloc(sizeof(double)*nlen);
	U1d.gam = (double*)malloc(sizeof(double)*nlen);
	U1d.mu = (double*)malloc(sizeof(double)*nlen);
	U1d.kt = (double*)malloc(sizeof(double)*nlen);
	U1d.flux = (double**)malloc(sizeof(double*)*nlen);
	for (int i=0; i<nlen; i++)
	{
		U1d.flux[i] = (double*)malloc(sizeof(double)*neqv);
	}
}

/*---------------------------------------------------
 * allocate memory for other variables
 * ------------------------------------------------*/
void allocateOthers()
{
	int nc = config1.ni*config1.nj;

	/*------------Other variables-------------*/
	rhs  = (double**)malloc(sizeof(double*)*nc);
	qo   = (double**)malloc(sizeof(double*)*nc);
	for(int ic=0; ic<nc; ic++)
	{
		qo[ic]   = (double*)malloc(sizeof(double)*neqv);
		rhs[ic]  = (double*)malloc(sizeof(double)*neqv);
	}
}