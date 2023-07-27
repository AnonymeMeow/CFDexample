/*
 * bound.c
 *
 *  Created on: 2014.01.14
 *
 */

#include"comm.h"

/*---------------------------------------------------
 * Apply boundary condition in i direction
 * ------------------------------------------------*/
void boundX()
{
	int i, ir, j, jr, ii, ns, ic, ic1;

	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

#ifdef MPI_RUN
	/*---- assign MPI boundary ----*/

	 /* left side, solid wall */
	if(MyID == 0)
	{
        ii = 2*config1.Ng - 1; // 5 ->  0, 1, 2
    	for(i=0; i<config1.Ng; i++)
    	{
    		for(j=config1.Ng; j<jr; j++)
		    {
    			ic = i*J0 + j;
		    	ic1 = ii*J0 + j;
		    	Ug.q[ic][0] =  Ug.q[ic1][0];
		    	Ug.q[ic][1] =  -Ug.q[ic1][1];
		    	Ug.q[ic][2] =  -Ug.q[ic1][2];
		    	Ug.q[ic][3] =  Ug.q[ic1][3];
		    	Ug.tem[ic]  =  Ug.tem[ic1];
		    	Ug.pre[ic]  =  Ug.pre[ic1];

				Ug.gam[ic]  =  Ug.gam[ic1];
			    Ug.mu[ic]   =  Ug.mu[ic1];
			    Ug.kt[ic]   =  Ug.kt[ic1];
	    		if(config1.gasModel != 0)
	    			for(ns = 0; ns<config1.nspec; ns++)
	    			{
	    				Ug.qs[ic][ns] =  Ug.qs[ic1][ns];
	    				Ug.di[ic][ns] =  Ug.di[ic1][ns];
	    			}
		    }
		    ii -= 1;
    	}
	}
	else
	{
		for(i=0; i<config1.Ng; i++)
			for(j=0; j<config1.nj; j++)
		    {
				ic = i*J0 + (j+config1.Ng); // the left block's N-3->0, N-2->1, N-1->2
				ic1 = i*config1.nj + j;

				Ug.q[ic][0] = mpiRecv_ql[ic1].rho;
		    	Ug.q[ic][1] = mpiRecv_ql[ic1].u;
		    	Ug.q[ic][2] = mpiRecv_ql[ic1].v;
		    	Ug.q[ic][3] = mpiRecv_ql[ic1].e;
		    	Ug.pre[ic]  = mpiRecv_ql[ic1].p;
		    	Ug.tem[ic]  = mpiRecv_ql[ic1].t;

		    	Ug.gam[ic]  = mpiRecv_ql[ic1].ga;
		    	Ug.mu[ic]   = mpiRecv_ql[ic1].mu;
		    	Ug.kt[ic]   = mpiRecv_ql[ic1].kt;
	    		if(config1.gasModel != 0)
	    			for(ns = 0; ns<config1.nspec; ns++)
	    			{
	    				Ug.qs[ic][ns] =  mpiRecv_ql[ic1].qs[ns];
	    				Ug.di[ic][ns] =  mpiRecv_ql[ic1].di[ns];
	    			}
		    }
	}

	/* right side, solid wall */
	if(MyID == NMAXproc)
	{
		ii = ir -1;  //N+2 ->  N+3,N+4,N+5
				for(i=ir; i<I0; i++)
				{
					for(j=config1.Ng; j<jr; j++)
				    {
						ic = i*J0 + j;
				    	ic1 = ii*J0 + j;

				    	Ug.q[ic][0] =  Ug.q[ic1][0];
				    	Ug.q[ic][1] =  -Ug.q[ic1][1];
				    	Ug.q[ic][2] =  -Ug.q[ic1][2];
				    	Ug.q[ic][3] =  Ug.q[ic1][3];
				    	Ug.tem[ic]  =  Ug.tem[ic1];
				    	Ug.pre[ic]  =  Ug.pre[ic1];

						Ug.gam[ic]  =  Ug.gam[ic1];
				    	Ug.mu[ic]   =  Ug.mu[ic1];
				    	Ug.kt[ic]   =  Ug.kt[ic1];
			    		if(config1.gasModel != 0)
			    			for(ns = 0; ns<config1.nspec; ns++)
			    			{
			    				Ug.qs[ic][ns] =  Ug.qs[ic1][ns];
			    				Ug.di[ic][ns] =  Ug.di[ic1][ns];
			    			}

				    }
					ii -= 1;
				}

	}
	else
	{
		// the right block's  0->N+3, 1->N+4, 2->N+5
		for(i=0; i<config1.Ng; i++)
		{
			ii = ir + i;
			for(j=0; j<config1.nj; j++)
			{
				ic = ii*J0 + j+config1.Ng;
				ic1 = i*config1.nj + j;

				Ug.q[ic][0] = mpiRecv_qr[ic1].rho;
				Ug.q[ic][1] = mpiRecv_qr[ic1].u;
				Ug.q[ic][2] = mpiRecv_qr[ic1].v;
				Ug.q[ic][3] = mpiRecv_qr[ic1].e;
				Ug.pre[ic]  = mpiRecv_qr[ic1].p;
				Ug.tem[ic]  = mpiRecv_qr[ic1].t;

				Ug.gam[ic]  = mpiRecv_qr[ic1].ga;
				Ug.mu[ic]   = mpiRecv_qr[ic1].mu;
				Ug.kt[ic]   = mpiRecv_qr[ic1].kt;
	    		if(config1.gasModel != 0)
	    			for(ns = 0; ns<config1.nspec; ns++)
	    			{
	    				Ug.qs[ic][ns] = mpiRecv_qr[ic1].qs[ns];
	    				Ug.di[ic][ns] = mpiRecv_qr[ic1].di[ns];
	    			}
			}
		}
	}
#else
	for(j=config1.Ng; j<jr; j++) // for actual j cells, 3 -> NJ + 2
	{
		/* left side, wall*/

		ii = 2*config1.Ng - 1; // 5, 4, 3 ->  0, 1, 2
		for(i=0; i<config1.Ng; i++)
		{
			ic  = i*J0 + j;
			ic1 = ii*J0 + j;

			Ug.q[ic][0] =  Ug.q[ic1][0];
			Ug.q[ic][1] = -Ug.q[ic1][1];
			Ug.q[ic][2] = -Ug.q[ic1][2];
			Ug.q[ic][3] =  Ug.q[ic1][3];
			Ug.tem[ic]  =  Ug.tem[ic1];
			Ug.pre[ic]  =  Ug.pre[ic1];

			Ug.gam[ic]  =  Ug.gam[ic1];
			Ug.mu[ic]   =  Ug.mu[ic1];
			Ug.kt[ic]   =  Ug.kt[ic1];
			Ug.c[ic]    = -1.;
			if(config1.gasModel != 0)
				for(ns=0; ns<config1.nspec; ns++)
				{
					Ug.qs[ic][ns] = Ug.qs[ic1][ns];
					Ug.di[ic][ns] = Ug.di[ic1][ns];
				}
			ii -= 1;
		}
			/* right side, wall*/
		ii = ir - 1;
		for(i=ir; i<I0; i++) //N+2, N+1, N ->  N+3,N+4,N+5
		{
			ic = i*J0 + j;
			ic1 = ii*J0 + j;

			Ug.q[ic][0] =  Ug.q[ic1][0];
			Ug.q[ic][1] = -Ug.q[ic1][1];
			Ug.q[ic][2] = -Ug.q[ic1][2];
			Ug.q[ic][3] =  Ug.q[ic1][3];
			Ug.tem[ic]  =  Ug.tem[ic1];
			Ug.pre[ic]  =  Ug.pre[ic1];

			Ug.gam[ic]  =  Ug.gam[ic1];
			Ug.mu[ic]   =  Ug.mu[ic1];
			Ug.kt[ic]   =  Ug.kt[ic1];
			Ug.c[ic]    = -1.;
			if(config1.gasModel != 0)
				for(ns=0; ns<config1.nspec; ns++)
				{
					Ug.qs[ic][ns] = Ug.qs[ic1][ns];
					Ug.di[ic][ns] = Ug.di[ic1][ns];
				}

			ii -= 1;
		}
	}
#endif
}

/*---------------------------------------------------
 * Apply boundary condition in j direction
 * ------------------------------------------------*/
void boundY()
{
	int i, ir, j, jr, jj, ns, ic, ic1;

	 ir = config1.ni + config1.Ng;
	 jr = config1.nj + config1.Ng;

	for(i=config1.Ng; i<ir; i++)
	{
		/* lower side, solid wall */
		jj = 2*config1.Ng - 1;
		for(j=0; j<config1.Ng; j++)
		{
					ic  = i*J0 + j;
					ic1 = i*J0 + jj;

		    		Ug.q[ic][0] =  Ug.q[ic1][0];
		    		Ug.q[ic][1] = -Ug.q[ic1][1];
		    		Ug.q[ic][2] = -Ug.q[ic1][2];
		    		Ug.q[ic][3] =  Ug.q[ic1][3];
					Ug.tem[ic]  =  Ug.tem[ic1];
					Ug.pre[ic]  =  Ug.pre[ic1];

		    		Ug.gam[ic]  =  Ug.gam[ic1];
		    		Ug.mu[ic]   =  Ug.mu[ic1];
		    		Ug.kt[ic]   =  Ug.kt[ic1];
		    		if(config1.gasModel != 0)
		    			for(ns = 0; ns<config1.nspec; ns++)
		    			{
		    				Ug.qs[ic][ns] =  Ug.qs[ic1][ns];
		    				Ug.di[ic][ns] =  Ug.di[ic1][ns];
		    			}
		    		jj -= 1;
		    }

    	/* upper side, plate */
		jj = jr -1;  //N+2, N+1, N ->  N+3,N+4,N+5
		for(j=jr; j<J0; j++)
		{
			ic = i*J0 + j;
			ic1 = i*J0 + jj;

		    Ug.q[ic][0] =  Ug.q[ic1][0];
		    Ug.q[ic][1] =  1.0;
		    Ug.q[ic][2] =  -Ug.q[ic1][2];
		    Ug.q[ic][3] =  Ug.q[ic1][3];
		    Ug.tem[ic]  =  Ug.tem[ic1];
		    Ug.pre[ic]  =  Ug.pre[ic1];

		    Ug.gam[ic]  =  Ug.gam[ic1];
		    Ug.mu[ic]   =  Ug.mu[ic1];
		    Ug.kt[ic]   =  Ug.kt[ic1];
	    	if(config1.gasModel != 0)
	    		for(ns = 0; ns<config1.nspec; ns++)
	    		{
	    			Ug.qs[ic][ns] =  Ug.qs[ic1][ns];
	    			Ug.di[ic][ns] =  Ug.di[ic1][ns];
	    		}
	    	jj -= 1;
		}
    }

}