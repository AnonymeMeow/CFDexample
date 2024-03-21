/*
 * jobbody.c
 *
 *  Created on: 2014.01.14
 *  last modified: Sep 10, 2014
 *
 */

#include<string.h>
#include<math.h>

#include<chrono>
#include<thread>

#include"comm.hpp"

void sync()
{
	static int count = 0;
	std::unique_lock<std::mutex> unique(mutex);
	count++;
	if (count == nproc)
	{
		count = 0;
		conditional.notify_all();
	}
	else
	{
		conditional.wait(unique);
	}
}

/*---------------------------------------------------
 * The main loop for the simulation
 * ------------------------------------------------*/
void jobbody()
{
	FILE *outId;

	void RKtvd3(int ik, double dt);
	void gettherm(int nc, double **q, double **qs, double *p,
			      double *t, double *gam1, double *rgas1, double *cv1);
	void gettrans(int nc, double **qs, double *t, double *gam1,
			      double *cv1, double *mu, double *cond, double **diff);
	void boundX();
	void boundY();
	void saveData(int step);
	void postprocess(int step);
	double getdt(int step);
	void allocateU1d(int);
	void allocateOthers();

	auto time_begin = std::chrono::steady_clock::now();
	outId = fopen("outInfo.dat", "a");

	int nc = config1.ni*config1.nj;
	double sum_t = 0.;

	std::thread *threads = new std::thread[nproc];

	for (int i = 0; i < nproc; i++)
	{
		threads[i] = std::thread(
			[&](int thread_id) {
				MyID = thread_id;
				allocateU1d(MAX(config1.ni + 2 * config1.Ng, J0));
				allocateOthers();

				int ic, ik, iv, ns, iStep;
				double dtc, Ttot;

				gettherm(nc, U.q, U.qs, U.pre, U.tem, U.gam, U.rgas, U.cv);
				gettrans(nc, U.qs, U.tem, U.gam, U.cv, U.mu, U.kt, U.di);
				for(iStep = config1.iStep0; iStep <= config1.nStep; iStep++)
				{
					for(ic = 0; ic<nc; ic++)
					{
						int ic1 = ic + MyID * nc;
						for(iv = 0; iv<neqv; iv++)
							qo[ic][iv] = U.q[ic1][iv];

						if(config1.gasModel != 0)
							for(ns = 0; ns<config1.nspec; ns++)
								qso[ic][ns] = U.qs[ic1][ns];
					}

					dtc = getdt(iStep);

					for(ik=0; ik<config1.timeOrder; ik++)
					{
						RKtvd3(ik, dtc);
						gettherm(nc, U.q, U.qs, U.pre, U.tem, U.gam, U.rgas, U.cv);
						gettrans(nc, U.qs, U.tem, U.gam, U.cv, U.mu, U.kt, U.di);
						boundX();
						boundY();
					}

					if(MyID == 0)
					{
						sum_t = sum_t + dtc;
						if(iStep%100 == 0)
						{
							// run on the cluster platform, output all the information to a file.
							fprintf(outId, "iStep: %d / total: %d\n",iStep, config1.nStep);
							fflush(outId);

							// run on the local computer, monitor output on the screen.
							//printf("iStep: %d / total: %d\n", iStep, config1.nStep);
						}
					}

					if(iStep%config1.Samples == 0)
					{
						saveData(iStep);

						if((MyID==0) && (iStep%config1.ifilm==0))
						{
							auto time_cpu = std::chrono::steady_clock::now() - time_begin;
							Ttot = config2.t0 + sum_t*tRef;

							printf("output flow field result, flow_time = %e s.\n", Ttot);
							fprintf(outId, "istep = %d, flow_time = %e s, cpu_time = %e s\n", iStep, Ttot, time_cpu.count() / 1e9);
							postprocess(iStep);
						}
					}
				}
			}, i
		);
	}

	for (int i = 0; i < nproc; i++)
	{
		threads[i].join();
	}

	printf("simulation complete! \n");
	fprintf(outId, "\n Program exit normally! \n");

	fclose(outId);
}

/*---------------------------------------------------
 * Time march use Runge-Kutta method
 * ------------------------------------------------*/
void RKtvd3(int ik, double dt)
{
	int nc;
	void flux(double **rhs);
	void fluxchem(double **rhs);
	void updateq(int nc, double **q, double **RHS, double dt, int ik);
	void updateqs(int nc, double **q, double **qs, double **RHS, double dt, int ik);
	void chemsource(int nc, double **q, double **qs, double *tem, double **RHS);

	nc = config1.ni*config1.nj;

	if(config1.gasModel == 0)
	{
		/* Perfect Gas Solver */
		flux(rhs);
		updateq(nc, U.q, rhs, dt, ik);
	}
	else
	{
		/* Real Gas Solver */
		fluxchem(rhs);
		if(config1.reacModel != 0)
			chemsource(nc, U.q, U.qs, U.tem, rhs);

		updateqs(nc, U.q, U.qs, rhs, dt, ik);
	}
}

/*---------------------------------------------------
 * Calculate the time step
 * ------------------------------------------------*/
double getdt(int step)
{
	double um, dti, CFL, dt, dx;

	if(config1.useDt)
	{
		// use constant dt
		if(step > config1.nRamp)
			dti = config2.dt1;
	    else
	    	dti = config2.dt0 + step*(config2.dt1 - config2.dt0)/config1.nRamp;

		dt = dti;
	}
	else
	{
		// use CFL number

		if(dxc < dyc)
			dx = dxc;
		else
			dx = dyc;

		if(step > config1.nRamp)
			CFL = config2.CFL1;
	    else
	    	CFL = config2.CFL0 + step*(config2.CFL1 - config2.CFL0)/config1.nRamp;

		if(config1.nonDi)
			dt = CFL*dx;
		else
		{
			um = config2.MaRef*sqrt(config2.gam0*(ru/config2.molWeight)*config2.temRef);
			dt = CFL*dx/um;
		}
	}
	return(dt);
}