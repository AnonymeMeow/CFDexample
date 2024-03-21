/*
 * save.c
 *
 *  Created on: 2014.1.19
 *
 */

#include<string.h>
#include<math.h>
#include"comm.hpp"
#include"chemdata.hpp"

/*---------------------------------------------------
 * output flow field variables
 * ------------------------------------------------*/
void saveData(int step)
{
	int ic, ic1, i, j, ns;
	char filename[20];
	double x, y, rho, u, v, p, T, e;
	FILE *fp;

	sprintf(filename, "tcv%d.dat", MyID);
	fp = fopen(filename, "w");

	for(j=0; j<config1.nj; j++)
		for(i=0; i<config1.ni; i++)
		{
			ic = (i + MyID * config1.ni) * config1.nj + j;
			ic1 = (i + MyID * config1.ni + config1.Ng) * J0 + j + config1.Ng;
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
}

/*---------------------------------------------------
 *  postprocess subroutine
 *  Created on: Apr 7, 2014
 * ------------------------------------------------*/
void postprocess(int istep)
{
	int i, j, ic, nc, ns, id, mni;
	double x, y, rho, u, v, p, T, e, dum;
	char filename[50];
	FILE  *fp;

	mni = nproc * config1.ni;
	nc = mni * config1.nj;

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
    fprintf(fp,"ZONE T='1', I= %d, J= %d, f=point \n", mni, config1.nj);

	for(j=0; j<config1.nj; j++)
		for(i=0; i<mni; i++)
		{
    		ic = i*config1.nj + j;
			int ic1 = (i + config1.Ng) * J0 + j + config1.Ng;
    		x = mesh.x[ic1];
    		y = mesh.y[ic1];

    		rho = U.q[ic][0];
    		u   = U.q[ic][1];
    		v   = U.q[ic][2];
    		p   = U.pre[ic];
    		T   = U.tem[ic];
    		e   = U.q[ic][3];
			fprintf(fp, "%8.4e %8.4e %6.4e %6.4e %6.4e %6.4e %10.4f %6.4e", x, y, rho, u, v, p, T, e);
			if(config1.gasModel != 0)
			{
				fprintf(fp, " %6.4f ", U.gam[ic]);
				for(ns=0; ns<config1.nspec; ns++)
					fprintf(fp, " %8.4e ", U.qs[ic][ns]);
			}
			fprintf(fp, "\n");
    	}
    fclose(fp);
    printf("nstep=%d postprocess complete!!! \n", istep);
}