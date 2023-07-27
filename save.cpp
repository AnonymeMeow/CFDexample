/*
 * save.c
 *
 *  Created on: 2014.1.19
 *
 */

#include<string.h>
#include<math.h>
#include"comm.h"
#include"chemdata.h"

/*---------------------------------------------------
 * output flow field variables
 * ------------------------------------------------*/
void saveData(int step)
{
	int ic, ic1, i, j, ns;
	char filename[20];
	double x, y, rho, u, v, p, T, e;
	FILE *fp;

#ifdef MPI_RUN
	sprintf(filename, "tcv%d.dat", MyID);
	fp = fopen(filename, "w");

	for(j=0; j<config1.nj; j++)
		for(i=0; i<config1.ni; i++)
		{
			ic  = i*config1.nj + j;
			ic1 = (i+config1.Ng)*J0 + j+config1.Ng;
			rho = U.q[ic][0];
			u  =  U.q[ic][1];
			v  =  U.q[ic][2];
			e  =  U.q[ic][3];
			p  =  U.pre[ic];
			T  =  U.tem[ic];
			x  =  mesh.xi[ic1];
			y  =  mesh.et[ic1];
			fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", x, y, rho, u, v, p, T, e);
			if(config1.gasModel != 0)
			{
				fprintf(fp, " %lf ", U.gam[ic]);
				for(ns=0; ns<config1.nspec; ns++)
					fprintf(fp, " %le ", U.qs[ic][ns]);
			}
			fprintf(fp, "\n");
		}
	fclose(fp);

#else

	sprintf(filename, "tcv%d.dat", MyID);

	fp = fopen(filename, "w");
	fprintf(fp, "Title = \"Flow field\"\n");
	fprintf(fp, "Variables = x, y, rho, u, v, p, T, e");
	if(config1.gasModel != 0)
	{
		fprintf(fp, ", gama");
		for(ns=0; ns<config1.nspec; ns++)
			fprintf(fp, ", %s", specData[ns].spname);
	}

	fprintf(fp, "\n");
	fprintf(fp,"ZONE T='1', I= %d, J= %d, f=point \n", config1.ni, config1.nj);

	/*-- tecplot IJ-ordered or finite-element surface data set --*/
	for(j = 0; j<config1.nj; j++)
		for(i = 0; i<config1.ni; i++)
		{
			ic = i*config1.nj + j;
			ic1 = (i+config1.Ng)*J0 + j+config1.Ng;
			rho = U.q[ic][0];
			u  =  U.q[ic][1];
			v  =  U.q[ic][2];
			e  =  U.q[ic][3];
			p  =  U.pre[ic];
			T  =  U.tem[ic];
			x  =  mesh.xi[ic1];
			y  =  mesh.et[ic1];
			fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", x, y, rho, u, v, p, T, e);
			if(config1.gasModel != 0)
			{
				fprintf(fp, " %lf ", U.gam[ic]);
				for(ns=0; ns<config1.nspec; ns++)
					fprintf(fp, " %le ", U.qs[ic][ns]);
			}
			fprintf(fp, "\n");
		}
    fclose(fp);
#endif

}

/*---------------------------------------------------
 *  postprocess subroutine
 *  Created on: Apr 7, 2014
 *  Author: chensong
 * ------------------------------------------------*/
void postprocess(int istep)

{
	int i, j, ic, nc, ns, id, mni, ni, nj, nproc, dum0;
	double x, y, rho, u, v, p, T, e, dum;
	char filename[50];
	FILE  *fp;
	struct strct_Up
	{
		double *rho, *u, *v, *p, *t, *e, *ga, **qs;
	} Up;
	void endjob();

	 /*--------1. read the grid information--------*/
	fp    = fopen("gridset.dat", "r");
	if(fp == NULL)
	{
		printf("gridset.dat file not found \n");
		endjob();
	}
	fscanf(fp, "%d %d %d %d", &ni, &nj, &dum0, &nproc);
	fclose(fp);

	 /*--------2. allocate memory--------*/
	mni = nproc*ni;
	nc  = mni*nj;
	Up.rho    = (double*)malloc(sizeof(double)*nc);
	Up.u      = (double*)malloc(sizeof(double)*nc);
	Up.v      = (double*)malloc(sizeof(double)*nc);
	Up.p      = (double*)malloc(sizeof(double)*nc);
	Up.t      = (double*)malloc(sizeof(double)*nc);
	Up.e      = (double*)malloc(sizeof(double)*nc);
	Up.ga     = (double*)malloc(sizeof(double)*nc);
	if(config1.gasModel != 0)
	{
		Up.qs = (double**)malloc(sizeof(double*)*nc);
		for(ic=0; ic<nc; ic++)
			Up.qs[ic] = (double*)malloc(sizeof(double)*config1.nspec);
	}
	else
		Up.qs = NULL;

    /*--------3. read the solution file--------*/
	 //printf("reading the solution file... \n");
    for(id=0; id<nproc; id++)
    {
    	sprintf(filename, "tcv%d.dat", id);
    	fp=fopen(filename,"r");
    	if(fp == NULL)
    	{
    		printf("plot file tcv%d.dat not found! \n", id);
    		printf("for single processor, the #define MPI_RUN in headfile 'comm.h' should be deactived \n");
    		exit(0);
    	}

		for(j=0; j<nj; j++) // without ghost cells
			for(i=0; i<ni; i++)
			{
				ic = (id*ni + i)*nj + j;
				fscanf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",&dum, &dum,
						&Up.rho[ic],&Up.u[ic],&Up.v[ic],&Up.p[ic],&Up.t[ic],&Up.e[ic]);
				if(config1.gasModel != 0)
				{
					fscanf(fp," %lf",&Up.ga[ic]);
					for(ns=0; ns<config1.nspec; ns++)
						fscanf(fp," %le",&Up.qs[ic][ns]);
				}
				fprintf(fp, "\n");
    	    }
        fclose(fp);
    }

	/*--------5. write the final field.dat--------*/
    sprintf(filename, "nstep%d_field.dat", istep);
    fp = fopen(filename, "w");
    fprintf(fp, "Title = \"Flow field\"\n");
	fprintf(fp, "Variables = x, y, rho, u, v, p, T, e");
	if(config1.gasModel != 0)
	{
		fprintf(fp, ", gama");
		for(ns=0; ns<config1.nspec; ns++)
			fprintf(fp, ", %s", spname[ns]);
	}
	fprintf(fp, "\n");
    fprintf(fp,"ZONE T='1', I= %d, J= %d, f=point \n", mni, nj);

	for(j=0; j<nj; j++)
		for(i=0; i<mni; i++)
		{
    		ic = i*nj + j;
    		x = mesh.x[ic];
    		y = mesh.y[ic];

    		rho = Up.rho[ic];
    		u   = Up.u[ic];
    		v   = Up.v[ic];
    		p   = Up.p[ic];
    		T   = Up.t[ic];
    		e   = Up.e[ic];
			fprintf(fp, "%8.4e %8.4e %6.4e %6.4e %6.4e %6.4e %10.4f %6.4e", x, y, rho, u, v, p, T, e);
			if(config1.gasModel != 0)
			{
				fprintf(fp, " %6.4f ", Up.ga[ic]);
				for(ns=0; ns<config1.nspec; ns++)
					fprintf(fp, " %8.4e ", Up.qs[ic][ns]);
			}
			fprintf(fp, "\n");
    	}
    fclose(fp);
    printf("nstep=%d postprocess complete!!! \n", istep);

    free(Up.rho);
    free(Up.u);
    free(Up.v);
    free(Up.p);
    free(Up.t);
    free(Up.e);
    if(config1.gasModel != 0)
    {
    	for(ic=0; ic<nc; ic++)
    	   free(Up.qs[ic]);
    	free(Up.qs);
    }
}

