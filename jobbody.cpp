/*
 * jobbody.c
 *
 *  Created on: 2014.01.14
 *  last modified: Sep 10, 2014
 *
 */

#include<string.h>
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
	void mpiConfig();
	void gettherm(int nc, double **q, double *p,
			      double *t, double *gam1, double *rgas1, double *cv1);
	void mpiSendrecv();
	void gettrans(int nc, double *t, double *gam1,
			      double *cv1, double *mu, double *cond);
	void boundX();
	void boundY();
	void saveData(int step);
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

		if(MyID == 0)
		{
			sum_t = sum_t + dtc;
			if(iStep%1000 == 0)
			{
				// run on the cluster platform, output all the information to a file.
				log_file << "iStep: " << iStep << " / total: " << config1.nStep << std::endl;
			}
		}

		if(iStep%config1.Samples == 0)
		{
			saveData(iStep);

			// MPI_Barrier(MPI_COMM_WORLD); Lock Here!

			if((MyID==0) && (iStep%config1.ifilm==0))
			{
				auto time_cpu = std::chrono::steady_clock::now() - time_begin;
				double Ttot = config2.t0 + sum_t;

				std::cout << "output flow field result, flow_time = " << Ttot << " s." << std::endl;
				log_file << "istep = " << iStep << ", flow_time = " << Ttot << " s, cpu_time = " << time_cpu.count() / 1e9 << std::endl;
				postprocess(iStep);
			}
		}
	}

	if(MyID == 0)
	{
		std::cout << "simulation complete! " << std::endl;
		log_file << "\n Program exit normally! " << std::endl;
	}
	log_file.close();
}

/*---------------------------------------------------
 * Time march use Runge-Kutta method
 * ------------------------------------------------*/
void RKtvd3(int ik, double dt)
{
	int nc;
	void flux(double **rhs);
	void updateq(int nc, double **q, double **RHS, double dt, int ik);

	nc = config1.ni*config1.nj;

	if(config1.gasModel == 0) /* Perfect Gas Solver */
	{
		flux(rhs);
		updateq(nc, U.q, rhs, dt, ik);
	}
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