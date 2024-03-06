/*
 * set.c
 *
 *  Created on: 2014.01.19
 *
 */
#include<math.h>
#include<mpi.h>
#include"comm.h"
#include"chemdata.h"
/*---------------------------------------------------
 * set and initialize the job
 * ------------------------------------------------*/
void setjob()
{
	int nc, nc1;

	void initjob();
	void importjob();
	void nondimen();
	void setGeom();

	void allocateU(int nlen, struct strct_U *U);
	void allocateUv();
	void allocateOthers();

	if(config1.gasModel == 0)
		neqn = neqv;
	else
		neqn = neqv + config1.nspec -1;

	nc  = config1.ni*config1.nj;
	nc1 = I0*J0;

	allocateU(nc, &U);
	allocateU(nc1, &Ug);
	allocateUv();
	allocateOthers();

	nondimen();

	if(config1.newrun == 1)
		initjob();
	else
		importjob();

	setGeom();
}

/*---------------------------------------------------
 * Conduct the initialization for the simulation
 * ------------------------------------------------*/
void initjob()
{
	int i, mni, i0, ir, ic, icp;
	double dis;
	void assigncells(int i1, int in, int j1, int jn, double u,
			         double v, double t, double p, double *qs);

	config2.t0 = 0;
	config1.iStep0 = 1;

	if(MyID == 0)
	{
		/* Note: only MyID==0 carries the mesh data */
		mni = nproc*config1.ni;
		dis = config2.x0*config2.Lx;
		for(i=0; i<mni; i++)
		{
			ic = i*config1.nj + 0;
			icp = (i+1)*config1.nj + 0;
			if((mesh.x[ic]<=dis) && (dis<mesh.x[icp]))
				break;
		}
	    config1.x_sh = i; 
	}
#ifdef MPI_RUN
	MPI_Bcast(&config1.x_sh, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
	// get the position of the diaphragm
    i0 = config1.x_sh/config1.ni; // the integer part
    ir = config1.x_sh%config1.ni; // the remainder part

#ifdef MPI_RUN

	if(MyID < i0)
		assigncells(0,config1.ni, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
	else if((MyID == i0) && (ir != 0))
	{
		assigncells(0,ir, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
		assigncells(ir,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p, inc[1].ys);
	} 
	else
		assigncells(0,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p, inc[1].ys);

	MPI_Barrier(MPI_COMM_WORLD);
#else

	assigncells(i0,ir, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
	assigncells(ir,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p,inc[1].ys);

#endif
}
/*-----------------------------------------------------------
 * Assign initial value to each cells
 * ---------------------------------------------------------*/
void assigncells(int i1, int in, int j1, int jn, double u,
		         double v, double t, double p, double *qs)
{
	int i, j, ic, ns;
	double ein, es, ek, rgas1, RT;
	double getenergy(double qs[], double t);
	double getes(int ns, double t);
	double getrgas(double qs[]);

	rgas1 = (ru/config2.molWeight)/rgasRef;

	for(i=i1; i<in; i++)
		for(j=j1; j<jn; j++)
		{
			ic = i*config1.nj + j;

			U.q[ic][1] = u;
			U.q[ic][2] = v;
			U.tem[ic]  = t;
			U.pre[ic]  = p;

			ek = 0.5*(u*u + v*v);

			if(config1.gasModel == 0)
			{
				RT = rgas1*t;
				U.q[ic][0] = p/RT;
				U.q[ic][3] = ek + RT/(config2.gam0-1)*Upsilon;
			}
			else
			{
				ein = 0.;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					U.qs[ic][ns] = qs[ns];
					es = getes(ns, t);
					ein = ein + qs[ns]*es;
				}
				U.q[ic][3] = ek + ein;
				rgas1 = getrgas(qs);
				U.q[ic][0] = p/(rgas1*t);
			}
		}
}
/*-----------------------------------------------------------
 * Import the former solution.
 * ---------------------------------------------------------*/
void importjob()
{
	int i, j, ic, ic1, id, ns, mni, nc;
	double x, y, rho, u, v, p, T, e, dum;
	char linebuf[200], varname[200], filename[30];
	struct strct_field
	{
		double *rho, *u, *v, *p, *t, *e, *ga, **qs;
	} Uf;
	FILE *fp;
	void endjob();

	mni = nproc*config1.ni;
	nc  = mni*config1.nj;

	if(MyID == 0)
	{
		Uf.rho = (double*)malloc(sizeof(double)*nc);
		Uf.u   = (double*)malloc(sizeof(double)*nc);
		Uf.v   = (double*)malloc(sizeof(double)*nc);
		Uf.p   = (double*)malloc(sizeof(double)*nc);
		Uf.t   = (double*)malloc(sizeof(double)*nc);
		Uf.e   = (double*)malloc(sizeof(double)*nc);
		if(config1.gasModel != 0)
		{
			Uf.ga   = (double*)malloc(sizeof(double)*nc);
			Uf.qs = (double**)malloc(sizeof(double*)*nc);
			for(ic=0; ic<nc; ic++)
				Uf.qs[ic] = (double*)malloc(sizeof(double)*config1.nspec);
		}
		else
			Uf.qs = NULL;

		sprintf(filename, "nstep%d_field.dat", config1.iStep0 -1);
		fp = fopen(filename,"r");
		if(fp == NULL)
		{
			printf("%s not found! \n", filename);
#ifdef MPI_RUN
			MPI_Abort( MPI_COMM_WORLD, 21);
#else
			endjob();
#endif
		}
		printf("reading the field file... \n");
		fgets(linebuf,sizeof(linebuf),fp);
		fgets(varname,sizeof(linebuf),fp);
		fgets(linebuf,sizeof(linebuf),fp);

		/*----1. read the field file----*/
		for(j=0; j<config1.nj; j++)
			for(i=0; i<mni; i++)
			{
				ic = i*config1.nj + j;
				fscanf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",&dum, &dum,
						&Uf.rho[ic], &Uf.u[ic], &Uf.v[ic], &Uf.p[ic], &Uf.t[ic], &Uf.e[ic]);
				if(config1.gasModel != 0)
				{
					fscanf(fp," %lf",&Uf.ga[ic]);
					for(ns=0; ns<config1.nspec; ns++)
						fscanf(fp," %le",&Uf.qs[ic][ns]);
				}
			}
		fclose(fp);

		/*----2. write the new tcv.dat file for each processors----*/
		for(id=0; id<nproc; id++)
		{
			sprintf(filename, "tcv%d.dat", id);
			fp = fopen(filename,"w");
			fprintf(fp, "Title = \"Flow field of local processor\"\n");
			fprintf(fp, "%s", varname);
			fprintf(fp,"ZONE T='%d', I= %d, J= %d, f=point \n", MyID, config1.ni, config1.nj);

			for(j=0; j<config1.nj; j++) // without ghost cells
				for(i = 0; i<config1.ni; i++)
				{
					ic = (id*config1.ni + i)*config1.nj + j;
					ic1 = (i+config1.Ng)*J0 + j+config1.Ng;

					x   = mesh.xi[ic1];
					y   = mesh.et[ic1];
					fprintf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",x, y,
							Uf.rho[ic],Uf.u[ic],Uf.v[ic],Uf.p[ic],Uf.t[ic],Uf.e[ic]);
					if(config1.gasModel != 0)
					{
						fprintf(fp," %lf",Uf.ga[ic]);
						for(ns=0; ns<config1.nspec; ns++)
							fprintf(fp," %le",Uf.qs[ic][ns]);
					}
					fprintf(fp, "\n");
				}
			fclose(fp);
		}
	}

#ifdef MPI_RUN
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	/*----3. Each processor read the solution file----*/
    sprintf(filename, "tcv%d.dat", MyID);
    fp = fopen(filename,"r");
    if(fp == NULL)
    {
    	printf("plot file tcv%d.dat not found! \n", MyID);
        endjob();
    }
	fgets(linebuf, sizeof (linebuf), fp);
	fgets(linebuf, sizeof (linebuf), fp);
	fgets(linebuf, sizeof (linebuf), fp); // skip the title

	for(j=0; j<config1.nj; j++)
		for(i=0; i<config1.ni; i++)
		{
			ic = i*config1.nj + j;
			if(fscanf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",&dum, &dum, &rho, &u, &v, &p, &T, &e) != 8)
			{
				printf("format error in tcv.dat \n");
#ifdef MPI_RUN
				MPI_Abort( MPI_COMM_WORLD, 22);
#else
				endjob();
#endif
			}
			if(config1.gasModel != 0)
			{
				fscanf(fp," %lf",&U.gam[ic]);
				for(ns=0; ns<config1.nspec; ns++)
					if(fscanf(fp," %le",&U.qs[ic][ns]) != 1){printf("format error in tcv.dat \n");endjob();}
			}

			U.q[ic][0] = rho;
			U.q[ic][1] = u;
			U.q[ic][2] = v;
			U.q[ic][3] = e;
			U.pre[ic]  = p;
			U.tem[ic]  = T;
		}
    fclose(fp);

	if(MyID == 0)
	{
		free(Uf.rho);
		free(Uf.u);
		free(Uf.v);
		free(Uf.p);
		free(Uf.t);
		free(Uf.e);
		if(config1.gasModel != 0)
		{
			for(ic=0; ic<nc; ic++)
				free(Uf.qs[ic]);
			free(Uf.qs);
		}
	}
}

/*---------------------------------------------------
 * Conduct the non-dimensional process
 * ------------------------------------------------*/
void nondimen()
{
	int ns;
	double rgas0, temp, ratio1, ratio2;
	FILE *outId;

	if(config1.gasModel == 0)
		rgas0  = (ru/config2.molWeight);
	else
	{
		temp = 0.;
		for(ns = 0; ns<config1.nspec; ns++)
		{
			temp = temp + specData[ns].qsin[1]/specData[ns].wm;
		}
		rgas0 = ru*temp;
	}

	if(config1.nonDi)
	{
		rgasRef = rgas0;
		LRef    = config2.Lx;
		temRef  = config2.temRef;
		preRef  = config2.preRef;
		ratio1  = config2.temRef/config2.suthC1;
		ratio2  = (config2.suthC1 + config2.suthC2)/(config2.temRef+ config2.suthC2);
		muRef   = config2.muRef*ratio2*pow(ratio1,1.5);

		uRef    = config2.MaRef*sqrt(config2.gam0*rgasRef*temRef);
		tRef    = LRef/uRef;
		cvRef   = rgasRef/(config2.gam0 - 1.);
		condRef = muRef*(config2.gam0*cvRef)/config2.Pr0;

		if(config2.Re0 < 10.)
		{
			rhoRef  = preRef/(rgas0*temRef);
			config2.Re0 = rhoRef*uRef*LRef/muRef;
			if(MyID == 0)
			{
				outId = fopen("outInfo.dat", "a");
				fprintf(outId, "Reynolds number is recalculated by reference condition: \n");
				fprintf(outId,"Re = %lf: \n",config2.Re0);
				fclose(outId);
			}
		}
		else
		{
			/*To keep consistency, the reference pressure is
			 * calculated by reference density, gas constant and temperature,
			 * which means the input "config2.preRef" is unused */
			rhoRef  = muRef*config2.Re0/(uRef*LRef);
			preRef  = rhoRef*rgasRef*temRef;
			if(MyID == 0)
			{
				outId = fopen("outInfo.dat", "a");
				fprintf(outId, "Under the given Reynolds number, the reference pressure is: \n");
				fprintf(outId,"P = %lf: (N/m^2)\n", preRef);
				fclose(outId);
			}
		}
		diffRef = muRef/(rhoRef*config2.Sc0);
		Upsilon = 1./(config2.gam0*config2.MaRef*config2.MaRef);
	}
	else
	{
		LRef    = 1.;
		uRef    = 1.;
		tRef    = 1.;
		cvRef   = 1.;
		rhoRef  = 1.;
		temRef  = 1.;
		preRef  = 1.;
		rgasRef = 1.;
		muRef   = 1.;
		condRef = 1.;
		Upsilon = 1.;
		diffRef = 1.;
		if((MyID==0) && (config1.nonDi!=0))
		{
			outId = fopen("outInfo.dat", "a");
			fprintf(outId, "the velocity for the dt via CFL number is calculated by: \n");
			fprintf(outId, "config2.MaRef*sqrt(config2.gam0*(ru/config2.molWeight)*config2.temRef) \n");
			fclose(outId);
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

	if(config1.visModel != 0)
	{
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
}
/*---------------------------------------------------
 * allocate memory for U vector
 * ------------------------------------------------*/
void allocateU(int nlen, struct strct_U *U)
{
	int i;
	U->mu   = (double*)malloc(sizeof(double)*nlen);
	U->kt   = (double*)malloc(sizeof(double)*nlen);
	U->cv   = (double*)malloc(sizeof(double)*nlen);
	U->rgas = (double*)malloc(sizeof(double)*nlen);
	U->pre  = (double*)malloc(sizeof(double)*nlen);
	U->tem  = (double*)malloc(sizeof(double)*nlen);
	U->gam  = (double*)malloc(sizeof(double)*nlen);
	U->q    = (double**)malloc(sizeof(double*)*nlen);

	for(i=0; i<nlen; i++)
		U->q[i]  = (double*)malloc(sizeof(double)*neqv);

	if(config1.gasModel == 0)
	{
		U->di = NULL;
		U->qs = NULL;
	}
	else
	{
		U->qs  = (double**)malloc(sizeof(double*)*nlen);
		U->di  = (double**)malloc(sizeof(double*)*nlen);
		for(i=0; i<nlen; i++)
		{
			U->qs[i] = (double*)malloc(sizeof(double)*config1.nspec);
			U->di[i] = (double*)malloc(sizeof(double)*config1.nspec);
		}
	}
}

/*---------------------------------------------------
 * allocate memory for derivatives
 * ------------------------------------------------*/
void allocateUv()
{
	int i, nlen;

	nlen = I0*J0;

	if(config1.visModel != 0)
	{
		// allocate geometry derivatives coefficients
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
		if(config1.gasModel == 0)
		{
			Uv.qs_xi = NULL;
			Uv.qs_et = NULL;
		}
		else
		{
			Uv.qs_xi = (double**)malloc(sizeof(double*)*nlen);
			Uv.qs_et = (double**)malloc(sizeof(double*)*nlen);
			for(i=0; i<nlen; i++)
			{
				Uv.qs_xi[i] = (double*)malloc(sizeof(double)*config1.nspec);
				Uv.qs_et[i] = (double*)malloc(sizeof(double)*config1.nspec);
			}
		}
	}
	else
	{
		Uv.fu1 = NULL;
		Uv.fu2 = NULL;
		Uv.fu3 = NULL;
		Uv.fv1 = NULL;
		Uv.fv2 = NULL;
		Uv.fv3 = NULL;
		Uv.fuv = NULL;
		Uv.fe1 = NULL;
		Uv.fe2 = NULL;
		Uv.gu1 = NULL;
		Uv.gu2 = NULL;
		Uv.gu3 = NULL;
		Uv.gv1 = NULL;
		Uv.gv2 = NULL;
		Uv.gv3 = NULL;
		Uv.guv = NULL;
		Uv.ge1 = NULL;
		Uv.ge2 = NULL;

		Uv.u_xi  = NULL;
		Uv.u_et  = NULL;
		Uv.v_xi  = NULL;
		Uv.v_et  = NULL;
		Uv.T_xi  = NULL;
		Uv.T_et  = NULL;
		Uv.qs_xi = NULL;
		Uv.qs_et = NULL;
	}
}

/*---------------------------------------------------
 * allocate memory for other variables
 * ------------------------------------------------*/
void allocateOthers()
{
	int ic, ns, nc;

	nc = config1.ni*config1.nj;

	/*------------Other variables-------------*/
	rhs  = (double**)malloc(sizeof(double*)*nc);
	qo   = (double**)malloc(sizeof(double*)*nc);
	for(ic=0; ic<nc; ic++)
	{
		qo[ic]   = (double*)malloc(sizeof(double)*neqv);
		rhs[ic]  = (double*)malloc(sizeof(double)*neqn);
	}
	if(config1.gasModel == 0)
		qso  = NULL;
	else
	{
		qso  = (double**)malloc(sizeof(double*)*nc);
		for(ic=0; ic<nc; ic++)
			qso[ic]  = (double*)malloc(sizeof(double)*config1.nspec);

		dsdq = (double***)malloc(sizeof(double**)*nc);
		for(ic=0; ic<nc; ic++)
		{
			dsdq[ic] = (double**)malloc(sizeof(double*)*config1.nspec);
			for(ns=0; ns<config1.nspec; ns++)
				dsdq[ic][ns] = (double*)malloc(sizeof(double)*neqn);
		}
	}
}