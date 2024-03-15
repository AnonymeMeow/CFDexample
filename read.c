/*
 * read.c
 *
 *  Created on: 2014.1.14.
 *
 */

#define MSMPI_NO_DEPRECATE_20

#include <ctype.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include"comm.h"
#include"chemdata.h"

/*---------------------------------------------------
 * read the necessary information for simulation
 * ------------------------------------------------*/
void readjob()
{
	void cheminput();
	void readconfig();
	void readic();
	void readmesh();
	void BcastData();

	neqv = 4; // solution variables without chemical terms

	// if(MyID == 0)
	// {
		readconfig();
		if(config1.gasModel != 0)
			cheminput();
	// }

	// BcastData();
	readic();
	readmesh();
}

/*---------------------------------------------------
 * read configuration information
 * ------------------------------------------------*/
void readconfig()
{
	FILE *fp, *outId;

	fp    = fopen("config.dat", "r");
	outId = fopen("outInfo.dat", "a");
	if(fp == NULL){printf("config.dat file not found! \n");exit(0);}

	if(fscanf(fp, "%lf %lf", &config2.t0, &config2.x0) != 2)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%d %d %d %d %d %d %d %d", &config1.newrun, &config1.nonDi, &config1.useDt,
			  &config1.iStep0, &config1.nStep, &config1.nRamp, &config1.Samples, &config1.ifilm) != 8)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf",  &config2.dt0, &config2.dt1, &config2.CFL0, &config2.CFL1) != 4)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%d %d %d %d", &config1.gasModel, &config1.reacModel, &config1.visModel,&config1.transModel )!= 4)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf %lf", &config2.molWeight,
			&config2.gam0, &config2.Pr0, &config2.Sc0, &config2.Re0) != 5)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf", &config2.muRef,  &config2.suthC1,  &config2.suthC2) != 3)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf", &config2.MaRef, &config2.temRef,  &config2.preRef) != 3)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf", &config2.p1, &config2.T1, &config2.u1, &config2.v1) != 4)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf", &config2.p2, &config2.T2, &config2.u2, &config2.v2) != 4)
			{printf("format error in config.dat\n");exit(0);}
	fclose(fp);

	fp = fopen("gridset.dat", "r");
	if(fp == NULL){printf("gridset.dat file not found \n");exit(0);}
	if(fscanf(fp, "%d %d %d %d", &config1.ni, &config1.nj,  &config1.Ng, &config1.nblock) != 4)
			 {printf("format error in gridset.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf", &config2.Lx, &config2.Ly) != 2)
			 {printf("format error in gridset.dat\n");exit(0);}
	fclose(fp);

	nproc = config1.nblock;

	fprintf(outId, "\n/---------------------------configure data---------------------------/\n");
	fprintf(outId, "\nt0=%lf, x0=%lf, Lx=%lf, Ly=%lf\n", config2.t0, config2.x0, config2.Lx, config2.Ly);
	fprintf(outId, "\ngrids point: ni=%d, nj=%d, ghost cells Ng=%d, nblock=%d\n", config1.ni, config1.nj, config1.Ng, config1.nblock);
	fprintf(outId, "\nnewrun=%d, nonDi=%d, useDt=%d, \niStep0=%d, nStep=%d, nRamp=%d, Samples=%d, ifilm=%d \n",
			        config1.newrun, config1.nonDi, config1.useDt, config1.iStep0, config1.nStep,
			        config1.nRamp, config1.Samples, config1.ifilm);
	fprintf(outId, "dt0=%lf, dt1=%lf, CFL0=%lf, CFL1=%lf, \n", config2.dt0, config2.dt1, config2.CFL0, config2.CFL1);
	fprintf(outId, "\ngasModel=%d, reacModel=%d, visModel=%d, transModel=%d \n", config1.gasModel,
			       config1.reacModel, config1.visModel,config1.transModel);
	fprintf(outId, "molWeight=%lf, gam0=%lf, pr0=%lf, Sc0=%lf, Re0=%lf\n", config2.molWeight,
			        config2.gam0, config2.Pr0, config2.Sc0, config2.Re0);
	fprintf(outId, "mu0=%lf, suthC1=%lf, suthC2=%lfn",config2.muRef, config2.suthC1, config2.suthC2);
	fprintf(outId, "Ma0=%lf, temRef=%lf,  preRef=%lf\n", config2.MaRef, config2.temRef, config2.preRef);

	fprintf(outId, "\nInitial condition: \n");
	fprintf(outId, "p_1=%lf, T_1=%lf, u_1=%lf, v1=%lf \n",config2.p1, config2.T1, config2.u1, config2.v1);
	fprintf(outId, "p_2=%lf, T_2=%lf, u_2=%lf, v2=%lf \n",config2.p2, config2.T2, config2.u2, config2.v2);

	config1.thermo_base = 0;
	config1.timeOrder = 3;

	fclose(outId);
}
/*---------------------------------------------------
 * read initial condition
 * ------------------------------------------------*/
void readic()
{
	int ns, ib, nb;

	nb = 2; // No. of initial condition

	inc[0].p = config2.p1;
	inc[0].t = config2.T1;
	inc[0].u = config2.u1;
	inc[0].v = config2.v1;
	inc[1].p = config2.p2;
	inc[1].t = config2.T2;
	inc[1].u = config2.u2;
	inc[1].v = config2.v2;
	if(config1.gasModel !=0)
	{
		for(ib=0; ib<nb; ib++)
			for(ns=0; ns<config1.nspec; ns++)
				inc[ib].ys[ns] = specData[ns].qsin[ib];
	}
}

/*---------------------------------------------------
 * read mesh data
 * ------------------------------------------------*/
void readmesh()
{
	int i, j, ic,nc, mni;
	char linebuf[200], filename[20];
	FILE  *fp, *outId;

	/*---1. Grid information---*/
	I0 = config1.ni + 2*config1.Ng;
	J0 = config1.nj + 2*config1.Ng;
	nc = I0*J0;

	/*---2. Allocate memory---*/
	mesh.xi   = (double*)malloc(sizeof(double)*nc);
	mesh.et   = (double*)malloc(sizeof(double)*nc);
	mesh.x_xi = (double*)malloc(sizeof(double)*nc);
	mesh.x_et = (double*)malloc(sizeof(double)*nc);
	mesh.y_xi = (double*)malloc(sizeof(double)*nc);
	mesh.y_et = (double*)malloc(sizeof(double)*nc);
	mesh.yaks = (double*)malloc(sizeof(double)*nc);

	sprintf(filename, "set%d.dat", MyID);
	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		printf("grid file set%d.dat not found! \n", MyID);
		exit(11);
	}

	fgets(linebuf, sizeof(linebuf), fp);
	fgets(linebuf, sizeof(linebuf), fp);
	fgets(linebuf, sizeof(linebuf), fp); // skip the title

	/*---3. Read the grids---*/
	for(j=0; j<J0; j++)
		for(i=0; i<I0; i++)
		{
			ic = i*J0 + j;
			if(fscanf(fp," %lf %lf %lf %lf %lf %lf %lf",
				  &mesh.xi[ic], &mesh.et[ic], &mesh.x_xi[ic], &mesh.x_et[ic],
				  &mesh.y_xi[ic], &mesh.y_et[ic], &mesh.yaks[ic]) != 7)
			{
				printf("format error in grid file set%d.dat! \n", MyID);
				exit(12);
			}
		}
	fclose(fp);
	
	/* After coordinate transformation,
	 * the distance between cells are equal */
	dxc = mesh.xi[1*J0] - mesh.xi[0*J0];
	dyc = mesh.et[1]    - mesh.et[0];

	/*---4. Read the physical mesh---*/

	if(MyID == 0)
	{
		outId = fopen("outInfo.dat", "a");

	    mni = nproc*config1.ni;
		nc = config1.nj*mni;
		mesh.x = (double*)malloc(sizeof(double)*nc);
		mesh.y = (double*)malloc(sizeof(double)*nc);

		fp = fopen("mesh.dat","r");
		if(fp == NULL)
		{
			printf("mesh file mesh.dat not found! \n");
			exit(13);
		}
		fgets(linebuf, sizeof (linebuf), fp);
		fgets(linebuf, sizeof (linebuf), fp);
		fgets(linebuf, sizeof (linebuf), fp); // skip the title
		
		for(j=0; j<config1.nj; j++)
			for(i=0; i<mni; i++)
			{
		    	ic = i*config1.nj + j;
		    	if(fscanf(fp," %lf  %lf",&mesh.x[ic], &mesh.y[ic])!= 2)
		    	{
					printf("format error in mesh file! \n");
					exit(14);
		    	}
			}

	    fclose(fp);
		fprintf(outId,"\nRead grid data complete!\n");
		fprintf(outId,"\n/---------------------runtime information---------------------/\n");
		fclose(outId);
	}
}
/*---------------------------------------------------
 * read chemical and thermal data
 * ------------------------------------------------*/
void cheminput()
{
	void cheminit();
	void therminit();
	void transinit();
	FILE *outId;

	outId = fopen("outInfo.dat", "a");

	cheminit();
	therminit();
	if(config1.transModel != 0)
		transinit();
	fprintf(outId,"\nRead thermo-chemical data complete!\n");
	fclose(outId);
}

/*---------------------------------------------------
 * Read the reaction data from model.chem
 * ------------------------------------------------*/
void cheminit()
{
	int im, ns;
	double sumqs;
	char linebuf[200], tempbuf[15];
	FILE *fp, *outId;

	char *trim(char *str);
	void reaction(FILE *fp);

	fp    = fopen("scarf_chem.dat", "r");
	outId = fopen("outInfo.dat", "a");

	if(fp == NULL){printf("chemical file not found"); exit(0);}

	fprintf(outId,"\n/---------------------------mixtures data---------------------------/\n");

	fgets(linebuf,sizeof(linebuf),fp); trim(linebuf);
	while(linebuf[0] == '*')
	{
		fgets(linebuf,sizeof(linebuf),fp);
		trim(linebuf);
	}

	if(sscanf(linebuf, "%d %d %d",  &config1.nmix, &config1.nspec, &config1.nreac) != 3 )
	{
		printf("format error in nmix, nspec, nreac\n");
		exit(0);
	}
	fprintf(outId, "\nnmix=%d, nspec=%d, nreac=%d \n", config1.nmix, config1.nspec, config1.nreac);

	for(ns=0; ns<config1.nspec; ns++)
	{
		if(fscanf(fp, "%s", tempbuf) != 1){printf("format error in species names \n");exit(0);}
		strcpy(spname[ns], trim(tempbuf));
		fprintf(outId,"species %d: %s \n", ns+1, spname[ns]);
	}

	if(fgets(linebuf, sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
	for(im=0; im<config1.nmix; im++)
	{
		if(fgets(linebuf, sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		fprintf(outId,"the mixture %d: \n",im+1);

		sumqs = 0.;
		for(ns=0; ns<config1.nspec-1; ns++)
		{
			if(fscanf(fp, "%le", &specData[ns].qsin[im]) != 1){printf("format error in species mass fraction \n");exit(0);}
			fprintf(outId,"%le  ", specData[ns].qsin[im]);
			sumqs = sumqs + specData[ns].qsin[im];
		}
		if(sumqs <= 1)
		{
			specData[config1.nspec-1].qsin[im] = 1. - sumqs;
			fprintf(outId,"%le", specData[config1.nspec-1].qsin[im]);
		}
		else
		{
			printf("input mixture error! \n");
			exit(15);
		}
		fprintf(outId,"\n");
		if(fgets(linebuf, sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
	}

	fclose(outId);
	reaction(fp);
}

/*----------------------------------------------------------------------------
*     Read the thermodynamic data file. Sets molecular weights,
*     heat capacity coefficients, heats of formation etc.
*----------------------------------------------------------------------------*/
void therminit()
{
	int nsets, istate, nd, ns, icount, isspec;
	double rmw, hform;
	char linebuf[200], tempname[15];
	FILE *fp, *outId;

	char *trim(char *str);

	fp    = fopen("scarf_glenn.dat", "r");
	outId = fopen("outInfo.dat", "a");

	if(fp == NULL)
	{
		printf("thermini.dat file not found! \n");
		exit(0);
	}

	fprintf(outId,"\n/---------------------------thermal data---------------------------/\n");

	icount = 0;
	while(fgets(linebuf, sizeof (linebuf), fp) != NULL)
	{
		trim(linebuf);
		if(linebuf[0] == '*' )
			continue;

		sscanf(linebuf, "%s", tempname);

		if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
		sscanf(linebuf, "%2d", &nsets);
		if(sscanf(&linebuf[50], "%2d %lf %lf", &istate, &rmw, &hform) != 3 )
		{
			printf("[scarf_glenn] format error in istate, rmw, hfrom\n");
			exit(16);
		}

		isspec = 0;
		for (ns=0; ns<config1.nspec; ns++)
		{
			if(strcmp(trim(tempname),trim(spname[ns])) == 0)
			{
				icount++;
				specData[ns].tempmax =     0.;
				specData[ns].tempmin = 30000.;
				specData[ns].ntemrng = nsets;
				specData[ns].wm      = rmw*1.e-3; // 1.e-3 transform wm to standard unit
				specData[ns].hof     = hform/specData[ns].wm;
				fprintf(outId,"\n species %d: %s\n",ns+1, spname[ns]);
				fprintf(outId,"nsets= %d, rmw=%f  hform=%f\n",nsets, rmw, hform);

				for(nd = 0; nd<nsets; nd++)
				{
					if(fgets(linebuf,sizeof (linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
					if( sscanf(linebuf, "%lf %lf", &specData[ns].temrng[nd][0], &specData[ns].temrng[nd][1]) != 2 )
					{
						fprintf(outId,"[therminit] format error in temrng\n");
						exit(0);
					}
					fprintf(outId,"temrng[0]=%lf  temrng[1]=%lf\n",specData[ns].temrng[nd][0],specData[ns].temrng[nd][1]);

					if(fgets(linebuf,sizeof (linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
					if( sscanf(linebuf, "%le%le%le%le%le",  &specData[ns].acoef[nd][0],
						&specData[ns].acoef[nd][1], &specData[ns].acoef[nd][2], &specData[ns].acoef[nd][3], &specData[ns].acoef[nd][4]) != 5 )
					{
						printf("[therminit] format error in acoef\n");
						exit(0);
					}
					fprintf(outId,"acoef0 -> 4= %16.9e %16.9e %16.9e %16.9e %16.9e\n", specData[ns].acoef[nd][0],
							specData[ns].acoef[nd][1], specData[ns].acoef[nd][2], specData[ns].acoef[nd][3], specData[ns].acoef[nd][4]);

					if(fgets(linebuf,sizeof (linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
					if( sscanf(linebuf, "%le%le",  &specData[ns].acoef[nd][5], &specData[ns].acoef[nd][6]) != 2 )
					{
						fprintf(outId,"[therminit] format error in acoef\n");
						exit(0);
					} // according to the format of NASA_glenn, this line should split and skip one data...
					if( sscanf(&linebuf[48], "%le%le", &specData[ns].acoef[nd][7], &specData[ns].acoef[nd][8]) != 2 )
					{
						fprintf(outId,"[therminit] format error in acoef\n");
						exit(0);
					}
					fprintf(outId,"acoef5 ->8 =%16.9e %16.9e %16.9e %16.9e \n", specData[ns].acoef[nd][5], specData[ns].acoef[nd][6],
							specData[ns].acoef[nd][7], specData[ns].acoef[nd][8]);

					specData[ns].tempmax = MAX(specData[ns].tempmax,specData[ns].temrng[nd][1]);
					specData[ns].tempmin = MIN(specData[ns].tempmin,specData[ns].temrng[nd][0]);
				}
				isspec = 1;
			    break;
			}
		}
		if(isspec == 0)
			for(nd = 0; nd<nsets; nd++)
			{
				// skip the lines
				if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
				if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
				if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in thermal data \n");exit(0);}
			}
		if(icount == config1.nspec)
		{
			fprintf(outId,"\nFinish reading thermodynamic data!\n");
			break;
		}
	} // end while
	if(icount < config1.nspec)
	{
		printf("\n species missing in the file thermal file !!!\n");
		exit(17);
	}
	fclose(fp);
	fclose(outId);
}

/*----------------------------------------------------------------------------
*     Read transport related data.
*----------------------------------------------------------------------------*/
void transinit()
{
	void muBlottner();
	void muCollision();

	if(config1.transModel == 1)
	{
		muBlottner();
	}
	else if(config1.transModel == 2)
	{
		muCollision();
	}
}
/*---------------------------------------------------
 * Read the viscosity coefficients of Blottner Model
 * ------------------------------------------------*/
void muBlottner()
{
	int ns, icount;
	double As, Bs, Cs;
	char linebuf[200], tempname[15];
	FILE *fp, *outId;

	char *trim(char *str);

	fp    = fopen("scarf_blottner.dat", "r");
	outId = fopen("outInfo.dat", "a");

	if(fp == NULL)
	{
		printf("scarf_blottner.dat not found! \n");
		exit(18);
	}

	fprintf(outId,"\n/--------viscosity coefficients of Blottner Model--------/\n");

	icount = 0;
	while(fgets(linebuf, sizeof (linebuf), fp) != NULL)
	{
		trim(linebuf);
		if(linebuf[0] == '*' )
			continue;

		if(sscanf(linebuf, "%s %le %le %le", tempname, &As, &Bs, &Cs) != 4 )
		{
			printf("[scarf_blottner.dat] format error! \n");
			exit(0);
		}

		for (ns=0; ns<config1.nspec; ns++)
		{
			if(strcmp(trim(tempname),trim(spname[ns])) == 0)
			{
				icount++;
				specData[ns].vis[0] = As;
				specData[ns].vis[1] = Bs;
				specData[ns].vis[2] = Cs;
				fprintf(outId,"\n %s:  ",spname[ns]);
				fprintf(outId,"A= %f, B=%f  C=%f\n",As, Bs, Cs);

				break;
			}
		}

		if(icount == config1.nspec)
		{
			fprintf(outId,"\nFinish reading viscosity data!\n");
			break;
		}
	} // end while
	if(icount < config1.nspec)
	{
		printf("\n species missing in the file trans.dat !!!\n");
		exit(19);
	}
	fclose(fp);
	fclose(outId);
}
/*---------------------------------------------------
 * Initialize transport data based on Lennard-Jones values for
 * characteristic molecular diameter and energy potential
 * ------------------------------------------------*/
void muCollision()
{
	int ns, icount;
	double tempsigma, tempepok;
	char linebuf[200], tempname[20];
	FILE *fp, *outId;

	char *trim(char *str);

	fp    = fopen("scarf_collision.dat", "r");
	outId = fopen("outInfo.dat", "a");

	if(fp == NULL)
	{
		printf("scarf_collision.dat not found! \n");
		exit(0);
	}

	fprintf(outId,"\n/--------Collision integral data for viscosity --------/\n");

	icount = 0;
	while(fgets(linebuf, sizeof (linebuf), fp) != NULL)
	{
		trim(linebuf);
		if(linebuf[0] == '*' )
			continue;

		if(sscanf(linebuf, "%s %le %le", tempname, &tempsigma, &tempepok) != 3 )
		{
			printf("[scarf_collision.dat] format error! \n");
			exit(110);
		}

		for (ns=0; ns<config1.nspec; ns++)
		{
			if(strcmp(trim(tempname),trim(spname[ns])) == 0)
			{
				icount++;
				specData[ns].vis[0] = tempsigma;
				specData[ns].vis[1] = 1./tempepok;
				fprintf(outId,"\n %s:  ",spname[ns]);
				fprintf(outId,"length = %f, L-J potential=%f \n",tempsigma, tempepok);

				break;
			}
		}

		if(icount == config1.nspec)
		{
			fprintf(outId,"\nFinish reading collision data!\n");
			break;
		}
	} // end while
	if(icount < config1.nspec)
	{
		printf("\n species missing in the scarf_collision.dat !!!\n");
		exit(122);
	}

	fclose(fp);
	fclose(outId);
}

/*----------------------------------------------------------------------
*     Read the species coefficients directly from reaction expression
*----------------------------------------------------------------------*/
void reaction(FILE *fp)
{
	int nr, ns, i, j, k, ilen, reacsp, ichar, iskip, nchar, schar;
	double rcoeff, ttvf;
	char linebuf[80], tempbuf[80], rsp[15];
	FILE *outId;
	void DBconvert(char *str, int ilen, char *temp, double coef);
	char *trim(char *str);

	outId = fopen("outInfo.dat", "a");
	fprintf(outId,"\n/---------------------------reaction data---------------------------/\n");

	nchar = sizeof(linebuf)/sizeof(linebuf[0]);
	schar = sizeof(rsp)/sizeof(rsp[0]);

	if((fgets(linebuf,sizeof(linebuf),fp)) == NULL)
	{
		printf("reaction name missing in chemical file \n");
		exit(112);
	}

    for(nr=0; nr<config1.nreac; nr++)
    {
		reacData[nr].thirdbody = 0;
		reacData[nr].backrate  = 0;
    	for(ns=0; ns<config1.nspec; ns++)
    	{
    		reacData[nr].nrsr[ns] = 0;
    		reacData[nr].npsr[ns] = 0;
    	}
    }

	for(nr=0; nr<config1.nreac; nr++)
	{
		if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		trim(linebuf);
	    ilen = strlen(linebuf);
	    fprintf(outId,"%s \n",linebuf);

	    reacsp = 1;
	    rcoeff = 1.;
	    iskip  = 0;
	    for(i = 0; i<ilen; i++)
	    {
	    	if( i < iskip )
	    		continue;
	        if(linebuf[i] != ' ') // skip the blank
	        {
	        	for(j = i; j<ilen; j++)
	        	{
	        		if((linebuf[j] == '+' || linebuf[j] == '=') && reacsp )  // the left of '=', namely, reactant
	                {
	        			if(linebuf[j] == '=')
	        				reacsp = 0; // stop at right hand side
	        			ichar = j-i;

	        			for(k = 0; k<nchar; k++)
	        				tempbuf[k] = '\0';  // nchar blank space
	        			for(k = 0; k<schar; k++)
	        				rsp[k] =  '\0';  // schar blank space


	        			for(k = 0; k<ichar; k++)
	        				tempbuf[k] = linebuf[i+k];  // fill the char between two +

	        			iskip = j+1; // skip the reactanct before +, jump to another reactant

	        			DBconvert(tempbuf, ichar, rsp, rcoeff);
	        			for(ns=0; ns<config1.nspec; ns++)
	        				if (strcmp(trim(spname[ns]),rsp) == 0)
	        				{
	        					reacData[nr].nrsr[ns]  = reacData[nr].nrsr[ns] + rcoeff;
	        					break;
	        				}

	        			if(strcmp("M",rsp) == 0)
	        				reacData[nr].thirdbody = 1;
	        			break;
					}

					//next, handle the right of '=', namely, production

					if( linebuf[j] == '+' || linebuf[j] == ';')
					{
						ichar = j-i;

						for(k = 0; k<nchar; k++)
							tempbuf[k] = '\0'; // nchar blank space
						for(k = 0; k<schar; k++)
							rsp[k] =  '\0'; // schar blank space

						for (k = 0; k<ichar; k++)
							tempbuf[k] = linebuf[i+k];

						iskip = j+1;
						DBconvert(tempbuf, ichar, rsp, rcoeff);

						for(ns=0; ns<config1.nspec; ns++)
							if(strcmp(spname[ns],rsp) == 0)
							{
								reacData[nr].npsr[ns]  = reacData[nr].npsr[ns] + rcoeff;
								break;
							}

						if(linebuf[j] == ';')
							goto Label;
						break; // finish reading a reaction
	                }
	        	}
	        } //end if(linebuf[i] != ' ')
		}

Label:  if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}

		/*--- read the "forward" reaction coeffiecients ---*/
		if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		if( sscanf(linebuf,"%lf %lf %lf %lf", &reacData[nr].af, &reacData[nr].nf, &reacData[nr].thetaf, &ttvf) != 4)
		{
			printf("[reaction] format error in af, nr, thetaf \n");
			exit(113);
		}
		fprintf(outId,"forward coefficients: \n");
		fprintf(outId,"af=%lf nf=%lf thetaf=%lf ttvf=%lf \n",reacData[nr].af, reacData[nr].nf, reacData[nr].thetaf, ttvf);

		if(reacData[nr].thirdbody)
		{
			if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
			trim(linebuf);
			fprintf(outId,"%s: \n", linebuf);

			if(linebuf[0] != 't' && linebuf[0] != 'T')
			{
				printf("third body coefficients not found! \n");
				exit(114);
			}
			for(ns=0; ns<config1.nspec; ns++)
			{
				if(fscanf(fp, "%lf", &reacData[nr].thrdeff[ns]) != 1){printf("format error in third-body coefficients\n");exit(0);}
				fprintf(outId,"%lf ", reacData[nr].thrdeff[ns]);
			}
			fprintf(outId,"\n");
			if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		}

		/*--- read the "backward" reaction coeffiecients ---*/
		if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
		{
			trim(linebuf);
			if(linebuf[0] == 'b' || linebuf[0] == 'B')
			{
				reacData[nr].backrate = 1;
				fprintf(outId,"backward coefficients: \n");
				if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
				if( sscanf(linebuf,"%lf %lf %lf \n", &reacData[nr].ab, &reacData[nr].nb, &reacData[nr].thetab) != 3)
				{
					printf("[reaction] format error in ab,nb, thetab \n");
					exit(115);
				}
				fprintf(outId,"ab=%lf nb=%lf thetab=%lf \n",reacData[nr].ab, reacData[nr].nb, reacData[nr].thetab);
			}
			else if(linebuf[0] == 'e' || linebuf[0] == 'E')
			{
				// use equilibrium constant
				reacData[nr].backrate = 2;
				fprintf(outId,"Equilibrium curve fit data: \n");
				if(fgets(linebuf,sizeof(linebuf),fp) == NULL){printf("format error in chemical data \n");exit(0);}
				if( sscanf(linebuf,"%lf %lf %lf %lf %lf\n", &reacData[nr].Br[0], &reacData[nr].Br[1],
						&reacData[nr].Br[2], &reacData[nr].Br[3], &reacData[nr].Br[4]) != 5)
				{
					printf("[reaction] format error in Br[1] ~ Br[5]\n");
					exit(116);
				}
				fprintf(outId,"Br[1] ~ Br[5]=%lf, %lf, %lf, %lf, %lf\n",reacData[nr].Br[0], reacData[nr].Br[1],
						reacData[nr].Br[2], reacData[nr].Br[3], reacData[nr].Br[4]);
			}
			else
			{
				fprintf(outId, "Equilibrium constant calculated by Gibbs free energy. \n");
			}
		}
	}
	fclose(outId);
	fclose(fp);
}

/*----------------------------------------------------------------------
*     function for trimming the white spaces
*----------------------------------------------------------------------*/
char *trim(char *str)
{
    size_t len = 0;
    char *frontp = str - 1;
    char *endp = NULL;

    if( str == NULL )
		return NULL;

    if( str[0] == '\0' )
		return str;

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address
     * the first non-whitespace characters from
     * each end.
     */
    while( isspace(*(++frontp)) );
    while( isspace(*(--endp)) && endp != frontp );

    if( str + len - 1 != endp )
		*(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
		*str = '\0';

    endp = str;
    if( frontp != str )
    {
    	while( *frontp ) *endp++ = *frontp++;
        *endp = '\0';
    }

    return str;
}

/*----------------------------------------------------------------------
* converts any character into a character and real value
* e.g.  string: [][]2.5N2[][], where [] represent blank
* input: str->string; ilen->length of char; mcha->length of string
* output: temp-> N2 ; coef-> 2.5
*----------------------------------------------------------------------*/
void DBconvert(char *str, int ilen, char *temp, double coef)
{
	int i, idec, ireal , ichar, iplace;
	double fact;
	char isupper[27], isdigit[10];
	char *k;

	char *trim(char *str);

	strcpy(isupper, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	strcpy(isdigit, "123456789");

	coef  = 1.0;
	idec =  0;
	ireal = 0;
	ichar = 0;

	for(i = 0; i<ilen; i++)
	{
		if((str[i]>='A') &&(str[i]<='Z')) break;
		ichar ++; // the position of the species name

		if(str[i] == '.' )
		{
			idec  = i;
			continue;
		}

		ireal ++;
		if(str[i] == ' ' ) ireal = i - 1;
	}

	if( idec == 0 ) idec = ireal+1;

	fact   = 0.;
	for(i = 0; i<=idec-1; i++)
	{
		k = strchr(isupper, str[i]);
		iplace = idec - (i+1);
		fact   = fact + (double)(k[0])*pow(10,iplace);
	}

	iplace = 0;
	for(i = idec+1; i<=ireal; i++)
	{
		k = strchr(isdigit, str[i] );
		iplace = iplace + 1;
		fact   = fact + (double)(k[0])*pow(0.1,iplace);
	}

	if(fact > 0)
		coef = fact;
	else
		printf("[DBconvert] species coefficient error!");

	for( i = ichar; i<ilen; i++)
	{
		iplace = i-ichar;
		temp[iplace] = str[i];
	}

	temp = trim(temp);
}

/*---------------------------------------------------
 * tell other tasks the thermal-chemical data
 * ------------------------------------------------*/
void BcastData()
{
	int count1, count2, block[2];

	MPI_Aint offset[2], extent;
	MPI_Datatype reacDatatype, specDatatype, oldtype[2];

	/*-----------Broadcast configuration data-----------*/
	count1 = sizeof(config1)/sizeof(int);
	count2 = sizeof(config2)/sizeof(double);
	MPI_Bcast(&config1, count1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&config2, count2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*-----------Broadcast thermo-chemical data-----------*/

	offset[0] = 0;
	oldtype[0] = MPI_INT;
	block[0] = 2;
	MPI_Type_extent(MPI_INT, &extent);
	offset[1] = 2 * extent;
	oldtype[1] = MPI_DOUBLE;
	block[1] = (sizeof(specData[0])-block[0])/sizeof(double);
	MPI_Type_struct(2, block, offset, oldtype, &specDatatype);

	offset[0] = 0;
	oldtype[0] = MPI_INT;
	block[0] = 2;
	MPI_Type_extent(MPI_INT, &extent);
	offset[1] = 2 * extent;
	oldtype[1] = MPI_DOUBLE;
	block[1] = (sizeof(reacData[0])-block[0])/sizeof(double);
	MPI_Type_struct(2, block, offset, oldtype, &reacDatatype);

	MPI_Type_commit(&reacDatatype);
	MPI_Type_commit(&specDatatype);

	MPI_Bcast(reacData, maxreac, reacDatatype, 0, MPI_COMM_WORLD);
	MPI_Bcast(specData, maxspec, specDatatype, 0, MPI_COMM_WORLD);
}