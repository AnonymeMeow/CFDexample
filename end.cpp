/*
 * endjob.c
 *
 *  Created on: 2014.01.19
 *
 */

#include"comm.h"

/*---------------------------------------------------
 * End the simulation, clear all the space
 * ------------------------------------------------*/
void endjob()
{
	void freeU(int nlen, struct strct_U *U);
	void freeUv();
	void freeU1d(int);
	void freeOthers();

	freeU(I0*J0, &U);

	freeUv();
	freeU1d(MAX(I0,J0));
	freeOthers();

	free(mesh.x);
	free(mesh.y);
	printf("program exits! \n");
}

/*---------------------------------------------------
 * free memory of U vector
 * ------------------------------------------------*/
void freeU(int nlen, struct strct_U *U)
{
	free(U->mu);
	free(U->kt);
	free(U->cv);
    free(U->rgas);
	free(U->pre);
	free(U->tem);
	free(U->gam);

	for(int i=0; i<nlen; i++)
		free(U->q[i]);

	free(U->q);
}

/*---------------------------------------------------
 * free memory of dU variables
 * ------------------------------------------------*/
void freeUv()
{
	free(Uv.fu1);
	free(Uv.fu2);
	free(Uv.fu3);
	free(Uv.fv1);
	free(Uv.fv2);
	free(Uv.fv3);
	free(Uv.fuv);
	free(Uv.fe1);
	free(Uv.fe2);
	free(Uv.gu1);
	free(Uv.gu2);
	free(Uv.gu3);
	free(Uv.gv1);
	free(Uv.gv2);
	free(Uv.gv3);
	free(Uv.guv);
	free(Uv.ge1);
	free(Uv.ge2);

	free(Uv.u_xi);
	free(Uv.u_et);
	free(Uv.v_xi);
	free(Uv.v_et);
	free(Uv.T_xi);
	free(Uv.T_et);
}

void freeU1d(int nlen)
{
	free(U1d.xix);
	free(U1d.xiy);
	free(U1d.etx);
	free(U1d.ety);
	free(U1d.yas);
	free(U1d.rho);
	free(U1d.du);
	free(U1d.dv);
	free(U1d.dt);
	free(U1d.u);
	free(U1d.v);
	free(U1d.p);
	free(U1d.e);
	free(U1d.t);
	free(U1d.gam);
	free(U1d.mu);
	free(U1d.kt);
	for (int i=0; i<nlen; i++)
	{
		free(U1d.flux[i]);
	}
	free(U1d.flux);
}

/*---------------------------------------------------
 * free memory of other variables
 * ------------------------------------------------*/
void freeOthers()
{
	int nc = config1.ni*config1.nj;

	free(mesh.xi);
	free(mesh.et);
	free(mesh.x_xi);
	free(mesh.x_et);
	free(mesh.y_xi);
	free(mesh.y_et);
	free(mesh.yaks);

	for(int ic=0; ic<nc; ic++)
	{
		free(qo[ic]);
		free(rhs[ic]);
	}
	free(qo);
	free(rhs);
}