/*
 * jobbody.c
 *
 *  Created on: 2014.01.14
 *  last modified: Sep 10, 2014
 *
 */

#include<math.h>

#include<fstream>
#include<iostream>
#include<chrono>

#include"comm.h"

/*---------------------------------------------------
 * The main loop for the simulation
 * ------------------------------------------------*/
void jobbody()
{
	std::ofstream log_file("outInfo.dat", std::ios_base::app);

	void RKtvd3(int ik, double dt);
	void gettherm(int nc, double **q, double *p,
			      double *t, double *gam1, double *rgas1, double *cv1);
	void gettrans(int nc, double *t, double *gam1,
			      double *cv1, double *mu, double *cond);
	void boundX();
	void boundY();
	void postprocess(int step);
	double getdt(int step);

	// mpiConfig();

	auto time_begin = std::chrono::steady_clock::now();

	int nc = config1.ni*config1.nj;
	gettherm(nc, U.q, U.pre, U.tem, U.gam, U.rgas, U.cv);
	gettrans(nc, U.tem, U.gam, U.cv, U.mu, U.kt);

	double sum_t = 0.;
	for(int iStep = config1.iStep0; iStep <= config1.nStep; iStep++)
	{
		for(int ic = 0; ic<nc; ic++)
		{
			for(int iv = 0; iv<neqv; iv++)
				qo[ic][iv] = U.q[ic][iv];
		}

		double dtc = getdt(iStep);

		for(int ik=0; ik<config1.timeOrder; ik++)
		{

			// mpiSendrecv();

			RKtvd3(ik, dtc);
			gettherm(nc, U.q, U.pre, U.tem, U.gam, U.rgas, U.cv);
			gettrans(nc, U.tem, U.gam, U.cv, U.mu, U.kt);
			boundX();
			boundY();
		}

		sum_t = sum_t + dtc;

		if(iStep%config1.Samples == 0)
		{
			// MPI_Barrier(MPI_COMM_WORLD); Lock Here!

			if(iStep%config1.ifilm==0)
			{
				auto time_cpu = std::chrono::steady_clock::now() - time_begin;
				double Ttot = config2.t0 + sum_t;

				std::cout << "output flow field result, flow_time = " << Ttot << " s." << std::endl;
				log_file << "istep = " << iStep << ", flow_time = " << Ttot << " s, cpu_time = " << time_cpu.count() / 1e9 << std::endl;
				postprocess(iStep);
			}
		}
	}

	std::cout << "simulation complete! " << std::endl;
	log_file << "\n Program exit normally! " << std::endl;
	log_file.close();
}

/*---------------------------------------------------
 * Time march use Runge-Kutta method
 * ------------------------------------------------*/
void RKtvd3(int ik, double dt)
{
	void flux(double **rhs);
	void updateq(double **q, double **RHS, double dt, int ik);

	flux(rhs);
	updateq(U.q, rhs, dt, ik);
}

/*---------------------------------------------------
 * Calculate the time step
 * ------------------------------------------------*/
double getdt(int step)
{
	// use CFL number
	double dx = MIN(dxc, dyc);
	double CFL;

	if(step > config1.nRamp)
		CFL = config2.CFL1;
	else
		CFL = config2.CFL0 + step*(config2.CFL1 - config2.CFL0)/config1.nRamp;
	
	double um = config2.MaRef*sqrt(config2.gam0*(ru/config2.molWeight)*config2.temRef);
	return CFL*dx/um;
}

/*------------------------------------------------------
 * update the conservative variables for perfect-gas
 * use Runge-kutta method
 * ------------------------------------------------*/
void updateq(double **q, double **rhs, double dt, int idt)
{
	int i, j, ii, jj, ic, ic1;
	double rho, rhoU, rhoV, rhoE, rrho, yas;
	double coef[3][3] = {
		{1., 3./4., 1./3.},
		{0, 1./4., 2./3.},
		{1., 1./4., 2./3.}
	};

	/* the (:,0), and (config1.ni-1, :) are assigned wall condition */

	for(i=0; i<config1.ni; i++)
	{
		ii = i + config1.Ng;
		for(j=0; j<config1.nj; j++)
		{
			jj  = j + config1.Ng;

			ic  = i*config1.nj + j;
			ic1 = ii*J0 + jj;
			yas = mesh.yaks[ic1];

			rho  = q[ic][0];
			rhoU = q[ic][0]*q[ic][1];
			rhoV = q[ic][0]*q[ic][2];
			rhoE = q[ic][0]*q[ic][3];

			q[ic][0] = coef[0][idt]*qo[ic][0]           + coef[1][idt]*rho  + coef[2][idt]*dt*rhs[ic][0]/yas;
			q[ic][1] = coef[0][idt]*qo[ic][0]*qo[ic][1] + coef[1][idt]*rhoU + coef[2][idt]*dt*rhs[ic][1]/yas;
			q[ic][2] = coef[0][idt]*qo[ic][0]*qo[ic][2] + coef[1][idt]*rhoV + coef[2][idt]*dt*rhs[ic][2]/yas;
			q[ic][3] = coef[0][idt]*qo[ic][0]*qo[ic][3] + coef[1][idt]*rhoE + coef[2][idt]*dt*rhs[ic][3]/yas;

			rrho = 1./q[ic][0];
			q[ic][1] = q[ic][1]*rrho;
			q[ic][2] = q[ic][2]*rrho;
			q[ic][3] = q[ic][3]*rrho;
		}
	}
}

/*---------------------------------------------------
 * Calculate the thermal properties
 * ------------------------------------------------*/
void gettherm(int nc, double **q, double *p,
		      double *t, double *gam1, double *rgas1, double *cv1)
{
	/*-- for calorically perfect gas--*/
	double rgas0 = (ru/config2.molWeight);
	for(int ic=0; ic<nc; ic++)
	{
		double temp = 0.5*(q[ic][1]*q[ic][1] + q[ic][2]*q[ic][2]);
		rgas1[ic] = rgas0;
		gam1[ic] = config2.gam0;
		cv1[ic]  = rgas0/(config2.gam0 - 1.);
		p[ic] = (gam1[ic] - 1)*(q[ic][3] - temp)*q[ic][0];
		t[ic] = p[ic]/q[ic][0]/rgas1[ic];
	}
}

/*--------------------------------------------------------------
 * Calculate transport properties(viscosity, conductivity, ...
 * -------------------------------------------------------------*/
void gettrans(int nc, double *t, double *gam1,
		      double *cv1, double *mu, double *cond)
{
	/*---------- 1.2 Sutherland's law ----------*/
	for(int ic=0; ic<nc; ic++)
	{
		double tem = t[ic];
		double cp  = gam1[ic]*cv1[ic];
		double ratio1 = tem/config2.suthC1;
		double ratio2 = (config2.suthC1 + config2.suthC2)/(tem + config2.suthC2);
		mu[ic] = config2.muRef*ratio2*pow(ratio1,1.5);
		cond[ic] = mu[ic]*cp/config2.Pr0;
	}
}

/*---------------------------------------------------
 *  postprocess subroutine
 *  Created on: Apr 7, 2014
 * ------------------------------------------------*/
void postprocess(int istep)
{
	char filename[50];
	FILE  *fp;

    sprintf(filename, "nstep%d_field.dat", istep);
    fp = fopen(filename, "w");
    fprintf(fp, "Title = \"Flow field\"\n");
	fprintf(fp, "Variables = x, y, rho, u, v, p, T, e\n");
    fprintf(fp,"ZONE T='1', I= %d, J= %d, f=point \n", config1.ni, config1.nj);

	for(int j=0; j<config1.nj; j++)
		for(int i=0; i<config1.ni; i++)
		{
    		int ic = i*config1.nj + j;
			int ic1 = (i + config1.Ng) * J0 + j + config1.Ng;
    		double x = mesh.x[ic1];
    		double y = mesh.y[ic1];

    		double rho = U.q[ic][0];
			double u  =  U.q[ic][1];
			double v  =  U.q[ic][2];
			double p  =  U.pre[ic];
			double T  =  U.tem[ic];
			double e  =  U.q[ic][3];
			fprintf(fp, "%8.4e %8.4e %6.4e %6.4e %6.4e %6.4e %10.4f %6.4e", x, y, rho, u, v, p, T, e);
			fprintf(fp, "\n");
    	}
    fclose(fp);
    printf("nstep=%d postprocess complete!!! \n", istep);
}