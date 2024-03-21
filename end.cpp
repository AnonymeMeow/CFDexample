/*
 * endjob.c
 *
 *  Created on: 2014.01.19
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include"comm.hpp"

/*---------------------------------------------------
 * End the simulation, clear all the space
 * ------------------------------------------------*/
void endjob()
{
	int nc, nc1;

	void freeU(int nlen, struct strct_U *U);
	void freeUv();
	void freeOthers();

	nc = config1.ni*config1.nj;
	nc1 = I0*J0;

	freeU(nc,  &U);
	freeU(nc1, &Ug);

	freeUv();
	freeOthers();

	if(MyID==0)
	{
		free(mesh.x);
		free(mesh.y);
		printf("program exits! \n");
	}
}

/*---------------------------------------------------
 * free memory of U vector
 * ------------------------------------------------*/
void freeU(int nlen, struct strct_U *U)
{
	int i;

	free(U->mu);
	free(U->kt);
	free(U->cv);
    free(U->rgas);
	free(U->pre);
	free(U->tem);
	free(U->gam);

	for(i=0; i<nlen; i++)
		free(U->q[i]);

	free(U->q);

	if(config1.gasModel != 0)
	{
		for(i=0; i<nlen; i++)
		{
			free(U->qs[i]);
			free(U->di[i]);
		}
	}
	free(U->qs);
	free(U->di);
}
/*---------------------------------------------------
 * free memory of dU variables
 * ------------------------------------------------*/
void freeUv()
{
	int i, nlen;

	nlen = I0*J0;

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

	if((config1.visModel != 0) && (config1.gasModel != 0))
	{
		for(i=0; i<nlen; i++)
		{
			free(Uv.qs_xi[i]);
			free(Uv.qs_et[i]);
		}
		free(Uv.qs_xi);
		free(Uv.qs_et);
	}
	else
	{
		free(Uv.qs_xi);
		free(Uv.qs_et);
	}
}
/*---------------------------------------------------
 * free memory of other variables
 * ------------------------------------------------*/
void freeOthers()
{
	int ic, ns, nc;

	nc = config1.ni*config1.nj;

	free(mesh.xi);
	free(mesh.et);
	free(mesh.x_xi);
	free(mesh.x_et);
	free(mesh.y_xi);
	free(mesh.y_et);
	free(mesh.yaks);

	for(ic=0; ic<nc; ic++)
	{
		free(qo[ic]);
		free(rhs[ic]);
	}
	free(qo);
	free(rhs);
	if(config1.gasModel != 0)
	{
		for(ic=0; ic<nc; ic++)
			free(qso[ic]);

		for(ic=0; ic<nc; ic++)
		{
			for(ns=0; ns<config1.nspec; ns++)
			{
				free(dsdq[ic][ns]);
			}
			free(dsdq[ic]);
		}
		free(dsdq);
	}
	free(qso);
}