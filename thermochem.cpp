/*
 * thermochem.c
 *
 *  Created on: 2014.01.14
 *
 */
#include<string.h>
#include<math.h>
#include"comm.hpp"
#include"chemdata.hpp"

/*---------------------------------------------------
 * Calculate the thermal properties
 * ------------------------------------------------*/
void gettherm(int nc, double **q, double **qs, double *p,
		      double *t, double *gam1, double *rgas1, double *cv1)
{
	int ic, ns, nr, nd, ntol, icount;
	double hs, tol, told, tnew, energy, dum, tdiff,
		   tem_min, tem_max, temp, temp2, temp3, temp4,
		   cp, cps, u, v, etot, tem, rgas0;
	double rwm[maxspec];
    double c2 = 0.5, c3 = 1./3., c4 = 0.25, c5 = 0.2;
    void endjob();

    /*-- Set a range for the temperature. Terminate the program
     * if un-physical conditions are encountered --*/
	tem_min = 20.0;
	tem_max = 20000.0;

	FILE *outId;

	if(config1.gasModel == 0)
	{
		/*-- for calorically perfect gas--*/
		rgas0 = (ru/config2.molWeight);
		for(ic=0; ic<nc; ic++)
		{
			temp = 0.5*(q[ic][1]*q[ic][1] + q[ic][2]*q[ic][2]);
			rgas1[ic] = rgas0/rgasRef;
			gam1[ic] = config2.gam0;
			cv1[ic]  = rgas0/(config2.gam0 - 1.)/cvRef;
			p[ic] = (gam1[ic] - 1)*(q[ic][3] - temp)*q[ic][0]/Upsilon;
			t[ic] = p[ic]/q[ic][0]/rgas1[ic];
			tem = t[ic]*temRef;
			if( !(tem_min<tem && tem<tem_max))
			{
				outId = fopen("outInfo.dat", "a");
				fprintf(outId, "T = %f (K) at at i=%d, j=%d\n", tem, ic/config1.nj, ic%config1.nj);
				fclose(outId);
				printf("T = %f (K) at at i=%d, j=%d\n", tem, ic/config1.nj, ic%config1.nj);
				printf("The temperature may be unphysical for ideal gas model \n");
				exit(41);
			}
		}
	}
	else
	{
		/*-- for none calorically perfect gas--*/
		ntol = 20;
		tol = 1.0e-10;

		for(ns = 0; ns<config1.nspec; ns++)
			rwm[ns] = 1.0/specData[ns].wm;

		/*-- calculate mixture heat capacity, ratio of
		 *-- specific heats and total energy per mass */
		switch(config1.thermo_base)
		{
		case 0 : // NASA_glenn database
			{
				for(ic = 0; ic<nc; ic++)
				{
					/*--recover dimension--*/
					tem = t[ic]*temRef;
					u = q[ic][1]*uRef;
					v = q[ic][2]*uRef;
					etot = q[ic][3]*uRef*uRef;

					temp = 0.;
					for(ns=0; ns<config1.nspec; ns++)
					{
						temp = temp + qs[ic][ns]*rwm[ns];
					}
					rgas1[ic] = ru*temp;

					icount = 0;
					dum = (etot - 0.5*(u*u + v*v));
					told = tem;
/*----------------------Start Newton's iteration Method----------------------*/
					do {
						nr = 0;
					    cp = 0. ;
						temp = 0.;
					    energy = 0. ;
					    for(ns = 0; ns<config1.nspec; ns++)
					    {
					    	cps = 0.;
					    	if( told > specData[ns].tempmax ) //exceeds maximum
					    	{
					    		nr   = specData[ns].ntemrng;
					            temp = specData[ns].tempmax;
					    	}
					        else if(told < specData[ns].tempmin) //exceeds minimum
					        {
					            nr   = 0;
					            temp = specData[ns].tempmin;
					    	}
					        else
					        {
					        	for(nd = 0; nd<=specData[ns].ntemrng; nd++)
					        		if((told >= specData[ns].temrng[nd][0])  && (told <= specData[ns].temrng[nd][1]) )
					        		{
					        			nr   = nd;
					                    temp = told;
					                    break;
					        		}
					        }
							temp2 = temp*temp;
							temp3 = temp*temp2;
							temp4 = temp*temp3;
							cps = ru*(specData[ns].acoef[nr][0]/temp2
									  + specData[ns].acoef[nr][1]/temp
									  + specData[ns].acoef[nr][2]
									  + specData[ns].acoef[nr][3]*temp
									  + specData[ns].acoef[nr][4]*temp2
									  + specData[ns].acoef[nr][5]*temp3
									  + specData[ns].acoef[nr][6]*temp4)*rwm[ns];

							hs  = ru*(-specData[ns].acoef[nr][0]/temp2
									  + specData[ns].acoef[nr][1]*log(temp)/temp
									  + specData[ns].acoef[nr][2]
									  + c2*specData[ns].acoef[nr][3]*temp
									  + c3*specData[ns].acoef[nr][4]*temp2
									  + c4*specData[ns].acoef[nr][5]*temp3
									  + c5*specData[ns].acoef[nr][6]*temp4
									  +    specData[ns].acoef[nr][7]/temp)*rwm[ns]*temp;

							hs = hs + cps*(told - temp);

							energy = energy + qs[ic][ns]*(hs - ru*rwm[ns]*told);
							cp = cp + qs[ic][ns]*cps;
					    }
					    cv1[ic] = cp - rgas1[ic];
					    tnew = told + (dum-energy)/(cv1[ic] + 1.e-20);

					    tdiff = tnew - told;
					    told = tnew; //
					    icount = icount + 1;

					//--Check for convergence or maximum trip count
					}while ((fabs(tdiff) > tol) && (icount <= ntol));
/*----------------------End Newton's iteration Method----------------------*/

					t[ic] = tnew;
					
					if(icount == ntol)
					{
						printf("----------\n");
						printf("myid%d Reached maximum iteration Nr. in determining temperature.\n",MyID);
						printf("T_new=%le, T_old=%le T_diff=%le at i=%d, j=%d. icount=%d\n", tnew, told, tdiff, ic/config1.nj, ic%config1.nj, icount);
						printf("----------\n");
					}
					if(!(tem_min<t[ic] && t[ic]<tem_max))
					{
						outId = fopen("outInfo.dat", "a");
						fprintf(outId, "----------\n");
						fprintf(outId, "myid%d Newton method T out of range. \n",MyID);
						fprintf(outId, "Tem = %e at i=%d, j=%d, ic=%d\n", U.tem[ic],  ic/config1.nj, ic%config1.nj,ic);
						fprintf(outId, "U.tem = %e \n", U.tem[ic]);
						fprintf(outId, "Maybe reduce time-step size...\n");
						fprintf(outId, "----------\n");
						fclose(outId);

						/* At this stage, one can set temperature to the last time step value,
							* i.e.  t[ic] = tem; But it is recommended to kill the program and
							* examine the running conditions*/
						exit(42);
					}
					/*--non-dimension--*/
					gam1[ic]  = cp/cv1[ic];
					t[ic]     = t[ic]/temRef;
					cv1[ic]   = cv1[ic]/cvRef;
					rgas1[ic] = rgas1[ic]/rgasRef;
					p[ic]     = q[ic][0]*rgas1[ic]*t[ic];
				}
			break;
			}
		default:
		{
			printf("the default thermo_base is NASA_glenn data\n");
			endjob();
			break;
    	}

    	}// end switch
    }
}

/*---------------------------------------------------
 * Calculate the mixture gas constant
 * ------------------------------------------------*/
double getrgas(double qs[])
{
	int ns;
	double rgas1;

	rgas1 = 0.0;
	if(config1.gasModel==0)
	{
		rgas1 =  (ru/config2.molWeight)/rgasRef;
	}
	else if(config1.gasModel==1)
	{
		for(ns = 0; ns<config1.nspec; ns++)
			rgas1 = rgas1 + qs[ns]/specData[ns].wm;

		rgas1 = ru*rgas1/rgasRef;
	}
	return(rgas1);
}

/*---------------------------------------------------
 * Calculate the internal energy of species 'ns'
 * ------------------------------------------------*/
double getes(int ns, double t)
{
	int nr, nd;
	double hs, es, temp, temp2, temp3, temp4, cps, rwm;
    double c2 = 0.5, c3 = 1./3., c4 = 0.25, c5 = 0.2;
    void endjob();

    t = t*temRef;

    if(config1.gasModel==1)
	{
    	nr = 0;
    	temp = 0.;
		rwm = 1.0/specData[ns].wm;
		if( t > specData[ns].tempmax)
		{
			nr   = specData[ns].ntemrng;
			temp = specData[ns].tempmax;
		}
		else if(t < specData[ns].tempmin)
		{
			nr   = 0;
			temp = specData[ns].tempmin;
		}
		else
		{
			for(nd = 0; nd<=specData[ns].ntemrng; nd++)
				if((t >= specData[ns].temrng[nd][0])  && (t <= specData[ns].temrng[nd][1]) )
				{
					nr   = nd;
					temp = t;
					break;
				}
		}
		temp2 = temp*temp;
		temp3 = temp*temp2;
		temp4 = temp*temp3;
		cps = ru*(specData[ns].acoef[nr][0]/temp2
				  + specData[ns].acoef[nr][1]/temp
				  + specData[ns].acoef[nr][2]
				  + specData[ns].acoef[nr][3]*temp
				  + specData[ns].acoef[nr][4]*temp2
				  + specData[ns].acoef[nr][5]*temp3
				  + specData[ns].acoef[nr][6]*temp4)*rwm;

		hs  = ru*(-specData[ns].acoef[nr][0]/temp2
				  + specData[ns].acoef[nr][1]*log(temp)/temp
				  + specData[ns].acoef[nr][2]
				  + c2*specData[ns].acoef[nr][3]*temp
				  + c3*specData[ns].acoef[nr][4]*temp2
				  + c4*specData[ns].acoef[nr][5]*temp3
				  + c5*specData[ns].acoef[nr][6]*temp4
				  +    specData[ns].acoef[nr][7]/temp)*rwm*temp;

		hs = hs + cps*(t-temp);
		es = hs - ru*rwm*t;
		es = es/(uRef*uRef); //non-dimension
	}
	else
	{
		printf("[getes] unsupported gasModel\n");
		exit(43);
	}

	return(es);
}
/*--------------------------------------------------------------
 * Calculate transport properties(viscosity, conductivity, ...
 * -------------------------------------------------------------*/
void gettrans(int nc, double **qs, double *t, double *gam1,
		      double *cv1, double *mu, double *cond, double **diff)
{
	int i, ic, ns, na, nb, nr;
	double cp, rho, ratio1, ratio2, tem, dum1, dum2, dum3, As, Bs, Cs,
	       mixwm, tmin, tmax, upper, lower, omegamkc, redtem, sumxs, patm,
	       avewmab, avesigab, avekepab, omegadab, tiny;
	double wm[maxspec], rwm[maxspec], mus[maxspec], phi[maxspec], xs[maxspec], ys[maxspec],
	       diffab[maxspec][maxspec];
	struct strct_coll
	{
		double tstar[82], omegamk[82], omegad[82];
	}coll;
    void endjob();
	void coll_integral(int len, double *f1, double *f2, double *f3);

	tiny = 1.e-20;

	if(config1.gasModel == 0)
	{
		/*---------- 1. ideal gas flow----------*/
		if(config1.visModel == 0)
			return;
		else if(config1.visModel == 1)
		{
			/*---------- 1.1 constant "mu" and "cond"----------*/
			cp  = config2.gam0*cvRef;
			ratio1  = config2.temRef/config2.suthC1;
			ratio2  = (config2.suthC1 + config2.suthC2)/(config2.temRef+ config2.suthC2);
			dum1  = config2.muRef*ratio2*pow(ratio1,1.5);
			dum2  = dum1*cp/config2.Pr0;
			for(ic=0; ic<nc; ic++)
			{
				mu[ic]   = dum1/muRef;
				cond[ic] = dum2/condRef;
			}
		}
		else if(config1.visModel == 2)
		{
			/*---------- 1.2 Sutherland's law ----------*/
			for(ic=0; ic<nc; ic++)
			{
				tem = t[ic]*temRef;
				cp  = gam1[ic]*(cv1[ic]*cvRef);
				ratio1 = tem/config2.suthC1;
				ratio2 = (config2.suthC1 + config2.suthC2)/(tem + config2.suthC2);
				mu[ic] = config2.muRef*ratio2*pow(ratio1,1.5);
				cond[ic] = mu[ic]*cp/config2.Pr0;
			}
			for(ic=0; ic<nc; ic++)
			{
				mu[ic]   = mu[ic]/muRef;
				cond[ic] = cond[ic]/condRef;
			}
		}
		else
		{
			printf("unsupported viscosity model... \n");
			exit(44);
		}
	}
	else
	{
		/*---------- 2. real gas flow----------*/
		if(config1.visModel == 0)
			return;

		for(ns = 0; ns<config1.nspec; ns++)
		{
			wm[ns]  = specData[ns].wm;
			rwm[ns] = 1./specData[ns].wm;
		}

		if(config1.transModel == 1)
		{
			/*---------- 2.1. Blonner's Curve fit data (for air)----------*/

			for(ic=0; ic<nc; ic++)
			{
				tem = t[ic]*temRef;
				rho = U.q[ic][0]*rhoRef;
				cp  = gam1[ic]*(cv1[ic]*cvRef);

				mixwm = 0.;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					ys[ns]  = qs[ic][ns];
					mixwm = mixwm + ys[ns]/wm[ns];
				}
				mixwm = 1./mixwm;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					As = specData[ns].vis[0];
					Bs = specData[ns].vis[1];
					Cs = specData[ns].vis[2];
					mus[ns] = 0.1*exp((As*log(tem) + Bs)*log(tem) + Cs);
					xs[ns] = ys[ns]*mixwm*rwm[ns];
				}
				// using Wilke's semi-empirical mixing rule
				for(ns = 0; ns<config1.nspec; ns++)
				{
					phi[ns] = 0.;
					for(nr = 0; nr<config1.nspec; nr++)
					{
						dum1 = sqrt(8.*(1. + wm[ns]*rwm[nr]));
						dum2 = sqrt(mus[ns]/mus[nr])*pow((wm[nr]*rwm[ns]),0.25);
						dum3 = pow((1. + dum2), 2);
						phi[ns] = phi[ns] + xs[nr]*dum3/dum1;
					}
					phi[ns] = 1./(phi[ns] + tiny);
				}

				// Mixture laminar viscosity
				mu[ic] = 0.;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					mu[ic] = mu[ic] + xs[ns]*mus[ns]*phi[ns];
				}
				cond[ic] = mu[ic]*cp/config2.Pr0;
				//assume all the species have the same diffusion coefficient
				diff[ic][0] = mu[ic]/(rho*config2.Sc0 + tiny);
			}
			for(ic=0; ic<nc; ic++)
			{
				mu[ic]   = mu[ic]/muRef;
				cond[ic] = cond[ic]/condRef;
				for(ns = 0; ns<config1.nspec; ns++)
			    {
					diff[ic][ns] = diff[ic][0]/diffRef;
			    }
			}
		}
		else if(config1.transModel == 2)
		{
			/* transport properties given by the kinetic theory
			 * ref: J.D. Anderson, 2th ed.,AIAA Inc. 2006
			 * page: 691-705 */

			coll_integral(82, coll.tstar, coll.omegamk, coll.omegad);

			/*--- 1. viscosity coefficient ---*/
			tmin = coll.tstar[0];
			tmax = coll.tstar[81];
			for(ic=0; ic<nc; ic++)
			{
				tem = t[ic]*temRef;
				rho = U.q[ic][0]*rhoRef;
				patm = U.pre[ic]*preRef/1.0132e5;
				cp  = gam1[ic]*(cv1[ic]*cvRef);
				mixwm = 0.;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					ys[ns]  = qs[ic][ns];
					mixwm = mixwm + ys[ns]/wm[ns];
				}
				mixwm = 1./mixwm;

				for(ns = 0; ns<config1.nspec; ns++)
				{
					redtem = tem*specData[ns].vis[1];
					// handle exceeds Min and Max limit
					if(redtem < tmin) redtem = coll.tstar[0];
					if(redtem > tmax) redtem = coll.tstar[81];

					// linear interpolation
					for(i=1; i<=81; i++)
					{
		                if(redtem <= coll.tstar[i])
		                {
		                	lower    = coll.tstar[i] - coll.tstar[i-1];
		                	upper    = coll.omegamk[i] - coll.omegamk[i-1];
		                	omegamkc = upper/lower*(redtem - coll.tstar[i-1]) + coll.omegamk[i-1];
		                	break;
		                }
					}
					mus[ns] = 2.6693e-06*sqrt((wm[ns]*1000)*tem);
					mus[ns] = mus[ns]/(specData[ns].vis[0]*specData[ns].vis[0]*omegamkc);
					xs[ns] = ys[ns]*mixwm*rwm[ns];
				}

				// using Wilke's semi-empirical mixing rule
				for(ns = 0; ns<config1.nspec; ns++)
				{
					phi[ns] = 0.;
					for(nr = 0; nr<config1.nspec; nr++)
					{
						dum1 = sqrt(8.*(1. + wm[ns]*rwm[nr]));
						dum2 = sqrt(mus[ns]/mus[nr])*pow((wm[nr]*rwm[ns]),0.25);
						dum3 = pow((1. + dum2), 2);
						phi[ns] = phi[ns] + xs[nr]*dum3/dum1;
					}
					phi[ns] = 1./(phi[ns] + tiny);
				}
				// Mixture laminar viscosity
				mu[ic] = 0.;
				for(ns = 0; ns<config1.nspec; ns++)
				{
					mu[ic] = mu[ic] + xs[ns]*mus[ns]*phi[ns];
				}
				//diff[ic][0] = mu[ic]/(rho*config2.Sc0 + tiny);

				/*--- 2. diffusive coefficient ---*/

				dum1  = 1.8583e-7*pow(tem,1.5)/patm;
				for(nb=0; nb<config1.nspec; nb++)
					for(na=0; na<config1.nspec; na++)
					{
						avewmab  = (rwm[na] + rwm[nb])/1000.; // convert the unit.
			            avesigab = 0.5*(specData[na].vis[0] + specData[nb].vis[0]);
			            avekepab = sqrt(specData[na].vis[1] * specData[nb].vis[1]);
			            redtem   = tem*avekepab;
					    if(redtem < tmin) redtem = coll.tstar[0];
						if(redtem > tmax) redtem = coll.tstar[81];
						for(i=1; i<=81; i++)
						{
							if(redtem <= coll.tstar[i])
				            {
								lower    = coll.tstar[i] - coll.tstar[i-1];
				                upper    = coll.omegad[i] - coll.omegad[i-1];
				                omegadab = upper/lower*(redtem - coll.omegad[i-1]) + coll.omegad[i-1];
				                break;
				            }
						}
						diffab[na][nb] = dum1*sqrt(avewmab)/(avesigab*avesigab*omegadab);
					}
				//  mixture diffusion coefficient
				for(ns=0; ns<config1.nspec; ns++)
				{
					sumxs = 0.;
					for(na=0; na<ns; na++)
						sumxs = sumxs + xs[na]/(diffab[na][ns] + tiny);
					for(na=ns+1; na<config1.nspec; na++)
						sumxs = sumxs + xs[na]/(diffab[na][ns] + tiny);

					diff[ic][ns] = (1.0 - xs[ns])/(sumxs + tiny);
				}

				/*--- 3. conductive coefficient ---*/
				cond[ic] = mu[ic]*cp/config2.Pr0;
			}
			for(ic=0; ic<nc; ic++)
			{
				mu[ic]   = mu[ic]/muRef;
				cond[ic] = cond[ic]/condRef;
				for(ns = 0; ns<config1.nspec; ns++)
			    {
					diff[ic][ns] = diff[ic][ns]/diffRef;
			    }
			}
		}
		else
		{
			printf("unsupported viscosity model... \n");
			exit(45);
		}
	}
}
/*---------------------------------------------------
 * collision integral, Hirschfelder, 1954
 * ------------------------------------------------*/
void coll_integral(int ik, double *ftstar, double *fomegamk, double *fomegad)
{
	int i;
	double tstar[] = {0.30,  0.35,   0.40,   0.45,   0.50,   0.55,  0.60,
		      0.65,  0.70 ,   0.75,   0.80,   0.85,   0.90,  0.95,
		      1.00,  1.05,   1.10,   1.15,   1.20,   1.25,  1.30,
		      1.35,  1.40,   1.45,   1.50,   1.55,   1.60,  1.65,
		      1.70,  1.75,   1.80,   1.85,   1.90,   1.95,  2.00,
		      2.10,  2.20,   2.30,   2.40,   2.50,   2.60,  2.70,
		      2.80,  2.90,   3.00,   3.10,   3.20,   3.30,  3.40,
		      3.50,  3.60,   3.70,   3.80,   3.90,   4.00,  4.10,
		      4.20,  4.30,   4.40,   4.50,   4.60,   4.70,  4.80,
		      4.90,  5.00,   6.00,   7.00,   8.00,   9.00, 10.00,
		     20.00, 30.00,  40.00,  50.00,  60.00,  70.00, 80.00,
		     90.00, 100.00, 200.00, 300.00, 400.00};
	double omegamk[] = {2.7850,2.6280,2.4920,2.3680,2.2570,2.1560, 2.0650,
			  1.9820,1.9080,1.8410,1.7800,1.7250,1.6750, 1.6290,
			  1.5870,1.5490,1.5140,1.4820,1.4520,1.4240, 1.3990,
			  1.3750,1.3530,1.3330,1.3140,1.2960,1.2790, 1.2640,
			  1.2480,1.2340,1.2210,1.2090,1.1970,1.1860, 1.1750,
			  1.1560,1.1380,1.1220,1.1070,1.0930,1.0810, 1.0690,
			  1.0580,1.0480,1.0390,1.0300,1.0220,1.0140, 1.0070,
		      0.9999,0.9932,0.9870,0.9811,0.9755,0.9700, 0.9649,
	          0.9600,0.9553,0.9507,0.9464,0.9422,0.9382, 0.9343,
		      0.9305,0.9269,0.8963,0.8727,0.8538,0.8379, 0.8242,
		      0.7432,0.7005,0.6718,0.6504,0.6335,0.6194, 0.6076,
		      0.5973,0.5882,0.5320,0.5016,0.4811};
	double omegad[] = {2.6620,2.4760,2.3180,2.1840,2.0660,1.9660,1.8770,
		      1.7980,1.7290,1.6670,1.6120,1.5620,1.5170,1.4760,
		      1.4390,1.4060,1.3750,1.3460,1.3200,1.2960,1.2730,
		      1.2530,1.2330,1.2150,1.1980,1.1820,1.1670,1.1530,
		      1.1400,1.1280,1.1160,1.1050,1.0940,1.0840,1.0750,
		      1.0570,1.0410,1.0260,1.0120,0.9996,0.9878,0.9770,
		      0.9672,0.9576,0.9490,0.9406,0.9328,0.9256,0.9186,
		      0.9120,0.9058,0.8998,0.8942,0.8888,0.8836,0.8788,
		      0.8740,0.8694,0.8652,0.8610,0.8568,0.8530,0.8492,
		      0.8456,0.8422,0.8124,0.7896,0.7712,0.7556,0.7424,
		      0.6640,0.6232,0.5960,0.5756,0.5596,0.5464,0.5352,
		      0.5256,0.5130,0.4644,0.4360,0.4170};

	for(i=0; i<ik; i++)
	{
		ftstar[i] = tstar[i];
		fomegamk[i] = omegamk[i];
		fomegad[i] = omegad[i];
	}
}
/*---------------------------------------------------
 * Calculate the chemical source terms and the
 * related jacobian matrix
 * ------------------------------------------------*/
void chemsource(int nc, double **q, double **qs, double *tem, double **rhs)
{
	int ic, i, j, ii, jj, ns, nr, m, ic1;
    double rdum, multf, multb, third, temp, sumrs, tiny, dkeqdt, tem1, rtem,
           Trigger_T, kfterm, kbterm, termU, termV, termE, zt, z2, alfa, yas,
           kf, kb, keq, gibbs[2];
    double rhos[maxspec], term[maxspec], rwm[maxspec], dkfdq[maxeqn], dkbdq[maxeqn];
    double dimen_factor[9] = {1.e6, 1.e6, 1.e6, 1.e6, 1.e6, 1.e6, 1.e6, 1., 1.};

	void gibbsenergy(int nr, int nc, double t, double f[]);
	void tem_derivative(double **dtdq);

	m  = config1.nspec;
	tiny = 1.e-20;

	// the lowest temperature to trigger the reaction.
	Trigger_T = 1800.;

    alfa = LRef/(rhoRef*uRef);

	for(ic = 0; ic<nc; ic++)
		for(ns = 0; ns<m; ns++)
		{
    		sour[ic][ns] = 0.;
			for(j = 0; j<neqn; j++)
				dsdq[ic][ns][j] = 0.;
		}

	/*--- 1. Calculated temperature derivatives ---*/

	tem_derivative(dtdq);

    for(ns = 0; ns<m; ns++)
         rwm[ns] = 1.0/(specData[ns].wm);

    switch(config1.gasModel)
    {
    case 1:
    {
		for(ic = 0; ic<nc; ic++)
		{
	    	tem1 = tem[ic]*temRef;
	    	for(ns = 0; ns<m; ns++)
	    		rhos[ns] = (q[ic][0]*rhoRef)*qs[ic][ns];

	    	rtem = 1./tem1;

			if(tem1 < Trigger_T)
				continue;

			for(nr = 0; nr<config1.nreac; nr++)
			{
				/*--- 2. Determine the forward reaction rate ---*/
    			kf = reacData[nr].af*pow(tem1,reacData[nr].nf)*
    				     exp(-reacData[nr].thetaf*rtem);

        		/*--- 3. Determine the backward reaction rate ---*/
        			if(reacData[nr].backrate == 1)
        			{
        				kb = reacData[nr].ab*pow(tem1,reacData[nr].nb)*
        				         exp(-reacData[nr].thetab*rtem);
        			}
        			else if(reacData[nr].backrate == 2)
        			{
        				//use park curve fit
            			zt = 10000*rtem;
            			z2 = zt*zt;
        				keq = exp(reacData[nr].Br[0] + reacData[nr].Br[1]*log(zt)
        						+ reacData[nr].Br[2]*zt + reacData[nr].Br[3]*z2
        						+ reacData[nr].Br[4]*zt*z2)*dimen_factor[nr];
        				kb  = kf/(keq + tiny);
        			}
        			else
        			{
        				//use Gibbs energy

        	    		gibbsenergy(nr, ic, tem1, gibbs);

        				sumrs = 0.;
        				for(ns = 0; ns<m; ns++)
        					sumrs = sumrs + reacData[nr].npsr[ns] - reacData[nr].nrsr[ns];
        				rdum = 101300.0/ru;
        				keq = pow((rdum*rtem),sumrs)*exp(-gibbs[0]);
        				kb  = kf/(keq + tiny);
        			}

    			if(reacData[nr].thirdbody)
    			{
                   /*---4.1 Determine chemical source term for
                    * reactions with third body efficiencies ---*/
    				multf = 1.;
    				multb = 1.;
                    third = 0.;
                    for(ns = 0; ns<m; ns++)
                    {
                    	rdum  = rhos[ns]*rwm[ns];
                    	third = third + reacData[nr].thrdeff[ns]*rdum;
                    	rdum  = MAX(rdum,tiny);
                    	multf = multf*pow(rdum,reacData[nr].nrsr[ns]);
                    	multb = multb*pow(rdum,reacData[nr].npsr[ns]);
                    }
                    for(ns = 0; ns<m; ns++)
                    	sour[ic][ns] = sour[ic][ns]
  									 + specData[nr].wm*(reacData[nr].npsr[ns]-reacData[nr].nrsr[ns])
  									 * third*(kf*multf-kb*multb)*alfa;
    			}
    			else
    			{
    				/*--- 4.2 Determine chemical source term for
    				 * reactions without third body efficiencies ---*/
    				multf = 1.;
                 	multb = 1.;
                 	for(ns = 0; ns<m; ns++)
                 	{
                 		rdum  = rhos[ns]*rwm[ns];
                 		rdum  = MAX(rdum,tiny);
                 		multf = multf*pow(rdum,reacData[nr].nrsr[ns]);
                 		multb = multb*pow(rdum,reacData[nr].npsr[ns]);
                 	}
                 	for(ns = 0; ns<m; ns++)
                 	{
                 		sour[ic][ns] = sour[ic][ns]
                 		             + specData[ns].wm*(reacData[nr].npsr[ns]-reacData[nr].nrsr[ns])
 									 *(kf*multf - kb*multb)*alfa;
                 	}
    			}

    			/*--- 5. Calculate the chemical source Jacobian ds/dq ---*/

    				//5.1 Determine d(k_f)/dU
    				rdum = kf*rtem*(reacData[nr].nf + reacData[nr].thetaf*rtem);

                    for(ns=0; ns<m; ns++)
                    	dkfdq[ns] = rdum*dtdq[ic][ns];

                	dkfdq[m]   = rdum*dtdq[ic][m];
                	dkfdq[m+1] = rdum*dtdq[ic][m+1];
                	dkfdq[m+2] = rdum*dtdq[ic][m+2];

                    //5.2 Determine d(k_b)/dU
                    if(reacData[nr].backrate == 1)
                        rdum = kb*rtem*(reacData[nr].nb + reacData[nr].thetab*rtem);
                    else if(reacData[nr].backrate == 2)
                    {
                    	 dkeqdt = keq*(reacData[nr].Br[1]/zt + reacData[nr].Br[2]
                    	                  + 2.*reacData[nr].Br[3]*zt + 3.*reacData[nr].Br[4]*z2);
                         rdum = kb*rtem*(reacData[nr].nf + reacData[nr].thetaf*rtem - dkeqdt/(keq + tiny));
                    }
                    else
                    {
                    	//use Gibbs energy
                        sumrs = 0.;
                        for(ns=0; ns<m; ns++)
                        	sumrs = sumrs + (reacData[nr].npsr[ns] - reacData[nr].nrsr[ns]);

                        dkeqdt = -keq*(sumrs*rtem  + gibbs[1]);
                        rdum = kb*rtem*(reacData[nr].nf + reacData[nr].thetaf*rtem - dkeqdt/(keq + tiny));
                    }

                	for(ns=0; ns<m; ns++)
                    	dkbdq[ns] = rdum*dtdq[ic][ns];
                	dkbdq[m]   = rdum*dtdq[ic][m];
                	dkbdq[m+1] = rdum*dtdq[ic][m+1];
                	dkbdq[m+2] = rdum*dtdq[ic][m+2];

                    //5.3 Calculate the concentration term
                	multf = 1.;
                	multb = 1.;
                	for(ns = 0; ns<m; ns++)
                	{
                		rdum  = rhos[ns]*rwm[ns];
                		rdum  = MAX(rdum,tiny);
                		multf = multf*pow(rdum,reacData[nr].nrsr[ns]);
                		multb = multb*pow(rdum,reacData[nr].npsr[ns]);
                	}
                    kfterm = kf*multf;
                    kbterm = kb*multb;

                    if(reacData[nr].thirdbody)
                    {
                    	third = 0.;
                        for(j=0; j<m; j++)
                        {
                        	rdum     = rhos[j]*rwm[j];
                        	third    = third + reacData[nr].thrdeff[j]*rdum;
                        	term[j] = reacData[nr].thrdeff[ns]*rwm[j]*(kfterm - kbterm);
                        }
                    }
                    else
                    {
                    	third = 1.;
                    	for(j=0; j<m; j++)
                    		term[j] = 0.;
                    }

                    //5.4 Determine the first term of dS/dU
                    for(j=0; j<m; j++)
                    {
                    	rdum  = MAX(rhos[j],tiny);
                    	term[j] = term[j] + third*(reacData[nr].nrsr[j]*kfterm
                    			 - reacData[nr].npsr[j]*kbterm)/rdum;
                    }

                    //5.5 Determine the second term of dS/dU
                    for(j=0; j<m; j++)
                    	term[j] = term[j] + third*(dkfdq[j]*multf - dkbdq[j]*multb);

                    termU = third*(dkfdq[m]*multf   - dkbdq[m]*multb);
                    termV = third*(dkfdq[m+1]*multf - dkbdq[m+1]*multb);
                    termE = third*(dkfdq[m+2]*multf - dkbdq[m+2]*multb);

                    /*5.6 Sum for all the reaction, and store the source Jacobian.
                         Mutiply non-dimensional factor for Jacobian term 'ds/dq' */
                    for(i=0; i<m; i++)
                    {
                    	temp = specData[i].wm*(reacData[nr].npsr[i] - reacData[nr].nrsr[i]);

                    	for(j=0; j<m; j++)
                    		dsdq[ic][i][j] = dsdq[ic][i][j] + temp*term[j]*(alfa*rhoRef);

                		dsdq[ic][i][m]   = dsdq[ic][i][m]   + temp*termU*(alfa*rhoRef*uRef);
                		dsdq[ic][i][m+1] = dsdq[ic][i][m+1] + temp*termV*(alfa*rhoRef*uRef);
                		dsdq[ic][i][m+2] = dsdq[ic][i][m+2] + temp*termE*(alfa*rhoRef*uRef*uRef);
                    }
			}// end nr -> nreac;
    	} // end ic -> nc;
    	break;
    } // end case1
    default:
    {
    	printf("unsupported chemical model [void chemsource]\n");
    	break;
    }
    }

    /*----- 6. Add source term to numerical right hand side-----*/

	for(i=0; i<config1.ni; i++)
	{
		ii = i + config1.Ng;
		for(j=0; j<config1.nj; j++) // without ghost cells
		{
			jj  = j + config1.Ng;

			ic  = i*config1.nj + j;
			ic1 = ii*J0 + jj;
			yas = mesh.yaks[ic1];
			for(ns=0; ns<m; ns++)
				rhs[ic][ns] = rhs[ic][ns] + sour[ic][ns]*yas;
		}
	}
}

/*--------------------------------------------------------------
 *  Determine the Gibbs free energy for the equilibrium constant
 * ------------------------------------------------------------*/
void gibbsenergy(int nr, int nc, double tem, double f[])
{
	int nd, ns, ntr , outofrange;
    double c2 = 0.5, c3 = 1./3., c4 = 0.25, c5 = 0.2;
    double temp, temp2, temp3, temp4, cps, gibbs, dgibbs;
    double hsbar, sbar, hsbar_im, dhsbardt, dsbardt, gibbsterm, dgibbsterm;
    void endjob();

    ntr = 0;
    temp = 0.;

    switch(config1.thermo_base)
    {
    case 0:
    {
    	gibbs  = 0.;
    	dgibbs = 0.;

        //-- check temperature range validity
		for(ns = 0; ns<config1.nspec; ns++)
		{
			outofrange = 0;
		    cps = 0.;
			if(tem > specData[ns].tempmax)
			{
				outofrange = 1;
				ntr   = specData[ns].ntemrng;
				temp = specData[ns].tempmax;
				temp2 = temp*temp;
				temp3 = temp*temp2;
				temp4 = temp*temp3;
				cps = (specData[ns].acoef[nr][0]/temp2
					+ specData[ns].acoef[nr][1]/temp
					+ specData[ns].acoef[nr][2]
					+ specData[ns].acoef[nr][3]*temp
					+ specData[ns].acoef[nr][4]*temp2
					+ specData[ns].acoef[nr][5]*temp3
					+ specData[ns].acoef[nr][6]*temp4);
			}
			else if(tem < specData[ns].tempmin)
			{
				outofrange = 1;
				ntr   = 0;
				temp = specData[ns].tempmin;
				temp2 = temp*temp;
				temp3 = temp*temp2;
				temp4 = temp*temp3;
				cps = (specData[ns].acoef[nr][0]/temp2
					+ specData[ns].acoef[nr][1]/temp
					+ specData[ns].acoef[nr][2]
					+ specData[ns].acoef[nr][3]*temp
					+ specData[ns].acoef[nr][4]*temp2
					+ specData[ns].acoef[nr][5]*temp3
					+ specData[ns].acoef[nr][6]*temp4);
			}
			else
			{
				for(nd = 0; nd<=specData[ns].ntemrng; nd++)
				if((tem >= specData[ns].temrng[nd][0])  && (tem <= specData[ns].temrng[nd][1]) )
				{
					ntr   = nd;
					temp = tem;
					temp2 = temp*temp;
					temp3 = temp*temp2;
					temp4 = temp*temp3;
					cps = 0.;
					break;
				}
			}
		    hsbar = (-specData[ns].acoef[ntr][0]/temp2
				+ specData[ns].acoef[ntr][1]*log(temp)/temp
				+ specData[ns].acoef[ntr][2]
				+ c2*specData[ns].acoef[ntr][3]*temp
				+ c3*specData[ns].acoef[ntr][4]*temp2
				+ c4*specData[ns].acoef[ntr][5]*temp3
				+ c5*specData[ns].acoef[ntr][6]*temp4
				+ specData[ns].acoef[ntr][7]/temp); // h/RT
		    hsbar_im = hsbar;

		    sbar  = (-c2*specData[ns].acoef[ntr][0]/temp2
		    		-     specData[ns].acoef[ntr][1]/temp
		    		+     specData[ns].acoef[ntr][2]*log(temp)
					+      specData[ns].acoef[ntr][3]*temp
					+  c2*specData[ns].acoef[ntr][4]*temp2
				    +  c3*specData[ns].acoef[ntr][5]*temp3
					+  c4*specData[ns].acoef[ntr][6]*temp4
					+     specData[ns].acoef[ntr][8]); // S/R

		    hsbar = (hsbar*temp + cps*(tem-temp))/tem;
		    sbar = sbar + cps*(log(tem/temp)); //
		    gibbsterm = hsbar - sbar;
		    gibbs = gibbs + (reacData[nr].npsr[ns]- reacData[nr].nrsr[ns])*gibbsterm;

	        /* Implicit terms for source Jacobian */
		    if(outofrange)
	        {
		    	dgibbsterm = (- temp*hsbar_im + cps*(tem + temp))/(tem*tem);
	        }
		    else
		    {
				dhsbardt = (2.*specData[ns].acoef[ntr][0]/temp3
						+ specData[ns].acoef[ntr][1]*(1.-log(temp))/temp2
						+ c2*specData[ns].acoef[ntr][3]
						+ 2.*c3*specData[ns].acoef[ntr][4]*temp
						+ 3.*c4*specData[ns].acoef[ntr][5]*temp2
						+ 4.*c5*specData[ns].acoef[ntr][6]*temp3
						-    specData[ns].acoef[ntr][7]/temp2);

		        dsbardt  = ( specData[ns].acoef[ntr][0]/temp3
						+     specData[ns].acoef[ntr][1]/temp2
						+     specData[ns].acoef[ntr][2]/temp
						+     specData[ns].acoef[ntr][3]
						+     specData[ns].acoef[ntr][4]*temp
						+     specData[ns].acoef[ntr][5]*temp2
						+     specData[ns].acoef[ntr][6]*temp3);
		        dgibbsterm = dhsbardt - dsbardt;
		    }
		    dgibbs = dgibbs + (reacData[nr].npsr[ns]- reacData[nr].nrsr[ns])*dgibbsterm;
    	}
    	break;
    }
	default:
	{
		printf("the default thermo_base is NASA_glenn data\n");
		endjob();
		break;
	}
    }// end switch

    f[0] = gibbs;
    f[1] = dgibbs;
}
/*---------------------------------------------------
 * calculate source term derivatives with respect to
 * species density & total density
 * ------------------------------------------------*/
void tem_derivative(double **dtdq)
{
	int ns, ic, m, nc;
	double rdum, rho, u, v, t, cv, es[maxspec];
	double getes(int ns, double t);

	m  = config1.nspec;
	nc = config1.ni*config1.nj;

	for(ic = 0; ic<nc; ic++)
	{
    	t   = U.tem[ic];
		for(ns = 0; ns<m; ns++)
			es[ns] = getes(ns,t)*uRef*uRef;

    	rho = U.q[ic][0]*rhoRef;
    	u   = U.q[ic][1]*uRef;
    	v   = U.q[ic][2]*uRef;
    	cv  = U.cv[ic]*cvRef;

		rdum = 1./(rho*cv);
		for(ns = 0; ns<m; ns++)
			dtdq[ic][ns] = (0.5*u*u - es[ns])*rdum;

		dtdq[ic][m]   = -rdum*u;
		dtdq[ic][m+1] = -rdum*v;
		dtdq[ic][m+2] =  rdum;
	}
}