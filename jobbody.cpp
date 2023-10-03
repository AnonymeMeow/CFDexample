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

#undef MPI_RUN

/*---------------------------------------------------
 * The main loop for the simulation
 * ------------------------------------------------*/
void jobbody()
{
	std::ofstream log_file("outInfo.dat", std::ios_base::app);

	void RKtvd3(int ik, double dt);
	void mpiConfig();
	void gettherm(int nc, double **q, double **qs, double *p,
			      double *t, double *gam1, double *rgas1, double *cv1);
	void mpiSendrecv();
	void gettrans(int nc, double **qs, double *t, double *gam1,
			      double *cv1, double *mu, double *cond, double **diff);
	void boundX();
	void boundY();
	void saveData(int step);
	void postprocess(int step);
	double getdt(int step);

#ifdef MPI_RUN
	mpiConfig();
#endif
	auto time_begin = std::chrono::steady_clock::now();

	int nc = config1.ni*config1.nj;
	gettherm(nc, U.q, U.qs, U.pre, U.tem, U.gam, U.rgas, U.cv);
	gettrans(nc, U.qs, U.tem, U.gam, U.cv, U.mu, U.kt, U.di);

	double sum_t = 0.;
	for(int iStep = config1.iStep0; iStep <= config1.nStep; iStep++)
	{
		for(int ic = 0; ic<nc; ic++)
		{
			for(int iv = 0; iv<neqv; iv++)
				qo[ic][iv] = U.q[ic][iv];

			if(config1.gasModel != 0)
				for(int ns = 0; ns<config1.nspec; ns++)
					qso[ic][ns] = U.qs[ic][ns];
		}

		double dtc = getdt(iStep);

		for(int ik=0; ik<config1.timeOrder; ik++)
		{

#ifdef MPI_RUN
			mpiSendrecv();
#endif
			RKtvd3(ik, dtc);
			gettherm(nc, U.q, U.qs, U.pre, U.tem, U.gam, U.rgas, U.cv);
			gettrans(nc, U.qs, U.tem, U.gam, U.cv, U.mu, U.kt, U.di);
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

			/*
			#ifdef MPI_RUN
			MPI_Barrier(MPI_COMM_WORLD); // Lock Here!
			#endif
			*/

			if((MyID==0) && (iStep%config1.ifilm==0))
			{
				auto time_cpu = std::chrono::steady_clock::now() - time_begin;
				double Ttot = config2.t0 + sum_t*tRef;

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
	void fluxchem(double **rhs);
	void updateq(int nc, double **q, double **RHS, double dt, int ik);
	void updateqs(int nc, double **q, double **qs, double **RHS, double dt, int ik);
	void chemsource(int nc, double **q, double **qs, double *tem, double **RHS);

	nc = config1.ni*config1.nj;

	if(config1.gasModel == 0) /* Perfect Gas Solver */
	{
		flux(rhs);
		updateq(nc, U.q, rhs, dt, ik);
	}
	else /* Real Gas Solver */
	{
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
	if(config1.useDt) // use constant dt
	{
		if(step > config1.nRamp)
			return config2.dt1;
		return config2.dt0 + step*(config2.dt1 - config2.dt0)/config1.nRamp;
	}
	// use CFL number
	double dx = MIN(dxc, dyc);
	double CFL;

	if(step > config1.nRamp)
		CFL = config2.CFL1;
	else
		CFL = config2.CFL0 + step*(config2.CFL1 - config2.CFL0)/config1.nRamp;

	if(config1.nonDi)
		return CFL*dx;
	double um = config2.MaRef*sqrt(config2.gam0*(ru/config2.molWeight)*config2.temRef);
	return CFL*dx/um;
}

#ifdef MPI_RUN
/*---------------------------------------------------
 * MPI configuration
 * ------------------------------------------------*/
void mpiConfig()
{
	NMAXproc = nthread - 1;
	if(MyID == NMAXproc)
		nachbar_rechts = MPI_PROC_NULL;
	else
		nachbar_rechts = MyID + 1;

	if(MyID == 0)
		nachbar_links = MPI_PROC_NULL;
	else
		nachbar_links = MyID - 1;
}

/*---------------------------------------------------
 * Conduct MPI send and receive
 * ------------------------------------------------*/
void mpiSendrecv()
{
    MPI_Status status;
	int counts = sizeof(mpiSend_ql)/(sizeof(double));

	if(counts < config1.Ng*config1.nj)
	{
		/* Since is not convenient to use Dynamic Array for MPI transfer,
		 * the size of mpiSend_ql, mpiSend_qr etc. is set as constant
		 * in the head file comm.h. As for very fine grid it may have to
		 * modified.   */
		printf("The size of MPI array is smaller than required! \n");
		printf("Modified the value of 'mpicell' in comm.h... \n");
		exit(0);
	}

	for(int i=0; i<config1.Ng; i++)
		for(int j=0; j<config1.nj; j++)
		{
			/* assign U.q[i][j], (i=0,1,2) to transfer block */
			int ic = i*config1.nj + j;

			mpiSend_ql[ic].rho = U.q[ic][0];
			mpiSend_ql[ic].u   = U.q[ic][1];
			mpiSend_ql[ic].v   = U.q[ic][2];
			mpiSend_ql[ic].e   = U.q[ic][3];
			mpiSend_ql[ic].p   = U.pre[ic];
			mpiSend_ql[ic].t   = U.tem[ic];

			mpiSend_ql[ic].ga  = U.gam[ic];
			mpiSend_ql[ic].mu  = U.mu[ic];
			mpiSend_ql[ic].kt  = U.kt[ic];
			if(config1.gasModel != 0)
				for(int ns=0; ns<config1.nspec; ns++)
				{
					mpiSend_ql[ic].qs[ns] = U.qs[ic][ns];
					mpiSend_ql[ic].di[ns] = U.di[ic][ns];
				}
		}

	for(int i=0; i<config1.Ng; i++)
	{
		int ii = config1.ni - config1.Ng + i;
		for(int j=0; j<config1.nj; j++)
		{
			int ic  =  i*config1.nj + j;
			int ic1 = ii*config1.nj + j; // N-3, N-2, N-1 -> 0, 1, 2
			mpiSend_qr[ic].rho = U.q[ic1][0];
			mpiSend_qr[ic].u   = U.q[ic1][1];
			mpiSend_qr[ic].v   = U.q[ic1][2];
			mpiSend_qr[ic].e   = U.q[ic1][3];
			mpiSend_qr[ic].p   = U.pre[ic1];
			mpiSend_qr[ic].t   = U.tem[ic1];

			mpiSend_qr[ic].ga  = U.gam[ic1];
			mpiSend_qr[ic].mu  = U.mu[ic1];
			mpiSend_qr[ic].kt  = U.kt[ic1];
			if(config1.gasModel != 0)
				for(int ns=0; ns<config1.nspec; ns++)
				{
					mpiSend_qr[ic].qs[ns] = U.qs[ic1][ns];
					mpiSend_qr[ic].di[ns] = U.di[ic1][ns];
				}
		}
	}

	// mpiSend_qr -> mpiRecv_ql
	MPI_Sendrecv(&mpiSend_qr, counts, MPI_DOUBLE, nachbar_rechts, 1, &mpiRecv_ql, counts, MPI_DOUBLE, nachbar_links,  1, MPI_COMM_WORLD, &status);
	MPI_Sendrecv(&mpiSend_ql, counts, MPI_DOUBLE, nachbar_links,  2, &mpiRecv_qr, counts, MPI_DOUBLE, nachbar_rechts, 2, MPI_COMM_WORLD, &status);

	MPI_Barrier(MPI_COMM_WORLD);
}
#endif