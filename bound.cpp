/*
 * bound.c
 *
 *  Created on: 2014.01.14
 *
 */

#include"comm.hpp"

/*---------------------------------------------------
 * Apply boundary condition in i direction
 * ------------------------------------------------*/
void boundX()
{
	int i, j, jr, ii, ns, ic, ic1;

	jr = config1.nj + config1.Ng;

	/*---- assign MPI boundary ----*/

	/* left side */
	if(MyID == 0)
	{
        ii = 2*config1.Ng - 1; // 5, 4, 3 ->  0, 1, 2
    	for(i=0; i<config1.Ng; i++)
    	{
    		for(j=config1.Ng; j<jr; j++)
		    {
    			ic = i*J0 + j;
		    	ic1 = ii*J0 + j;
#if defined Cavity || defined Tube
		    	Ug.q[ic][0] =  Ug.q[ic1][0];
		    	Ug.q[ic][1] =  -Ug.q[ic1][1]; // solid wall
		    	Ug.q[ic][2] =  -Ug.q[ic1][2]; // solid wall
		    	Ug.q[ic][3] =  Ug.q[ic1][3];
		    	Ug.tem[ic]  =  Ug.tem[ic1]; // solid wall
		    	Ug.pre[ic]  =  Ug.pre[ic1]; // solid wall
#elif defined Plate
		    	Ug.q[ic][0] =  Ug.q[ic1][0];
		    	Ug.q[ic][1] =  1361.12; // inlet free
		    	Ug.q[ic][2] =  0.; // inlet free
		    	Ug.q[ic][3] =  Ug.q[ic1][3];
		    	Ug.tem[ic]  =  288.16; // inlet free
		    	Ug.pre[ic]  =  1.01325e5; // inlet free
#endif

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

	/* right side */
	if(MyID == nproc - 1)
	{
		int Ir = config1.ni * nproc + config1.Ng;
		ii = Ir - 1; //N+2, N+1, N -> N+3, N+4, N+5
		for (i = Ir; i < I0; i++)
		{
			for (j = config1.Ng; j < jr; j++)
			{
				ic = i * J0 + j;
				ic1 = ii * J0 + j;

				Ug.q[ic][0] =  Ug.q[ic1][0];
#if defined Cavity || defined Tube
				Ug.q[ic][1] =  -Ug.q[ic1][1]; // solid wall
				Ug.q[ic][2] =  -Ug.q[ic1][2]; // solid wall
#elif defined Plate
				Ug.q[ic][1] =  Ug.q[ic1][1]; // inlet free
				Ug.q[ic][2] =  Ug.q[ic1][2]; // inlet free
#endif
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
}

/*---------------------------------------------------
 * Apply boundary condition in j direction
 * ------------------------------------------------*/
void boundY()
{
	int i, ir, j, jr, jj, ns, ic, ic1;

	ir = config1.ni + config1.Ng;
	jr = config1.nj + config1.Ng;

#if defined Cavity || defined Tube
	for(i=config1.Ng; i<ir; i++)
	{
		/* lower side, solid wall */
		jj = 2*config1.Ng - 1;
		for(j=0; j<config1.Ng; j++)
		{
			ic = (i + MyID * config1.ni) * J0 + j;
			ic1 = (i + MyID * config1.ni) * J0 + jj;

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

    	/* upper side */
		jj = jr -1;  //N+2, N+1, N ->  N+3,N+4,N+5
		for(j=jr; j<J0; j++)
		{
			ic = (i + MyID * config1.ni) * J0 + j;
			ic1 = (i + MyID * config1.ni) * J0 + jj;

		    Ug.q[ic][0] =  Ug.q[ic1][0];
#ifdef Cavity
		    Ug.q[ic][1] =  1.0; // plate
		    Ug.q[ic][2] =  -Ug.q[ic1][2]; // plate
#elif defined Tube
		    Ug.q[ic][1] =  Ug.q[ic1][1];
		    Ug.q[ic][2] =  Ug.q[ic1][2];
#endif
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
#elif defined Plate
	if(MyID <= 2)
	{
		for(i=config1.Ng; i<ir; i++)
	    {
		    jj = 2*config1.Ng - 1;// 5 -> 0, 1, 2
		    for(j=0; j<config1.Ng; j++)
		    {
				ic = (i + MyID * config1.ni) * J0 + j;
				ic1 = (i + MyID * config1.ni) * J0 + jj;

                //left down
		    	Ug.q[ic][0] =  Ug.q[ic1][0];
		    	Ug.q[ic][1] =  Ug.q[ic1][1];
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

    	    /* upper side */
		    jj = jr -1;  //N+2 ->  N+3,N+4,N+5
		    for(j=jr; j<J0; j++)
		    {
				ic = (i + MyID * config1.ni) * J0 + j;
				ic1 = (i + MyID * config1.ni) * J0 + jj;

		        Ug.q[ic][0] =  Ug.q[ic1][0];
		        Ug.q[ic][1] =  Ug.q[ic1][1];
		        Ug.q[ic][2] =  Ug.q[ic1][2];
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
	else
	{
		for(i=config1.Ng; i<ir; i++)
	    {
		    jj = 2*config1.Ng-1;
		    for(j=0; j<config1.Ng; j++)
		    {
				ic = (i + MyID * config1.ni) * J0 + j;
				ic1 = (i + MyID * config1.ni) * J0 + jj;

                //right down
		    	Ug.q[ic][0] =  Ug.q[ic1][0];
		    	Ug.q[ic][1] = -Ug.q[ic1][1];
		    	Ug.q[ic][2] = -Ug.q[ic1][2];
		    	Ug.q[ic][3] =  Ug.q[ic1][3];
				// Ug.tem[ic]  =  Ug.tem[ic1];//adiabatic
				Ug.tem[ic]  =  576.32-Ug.tem[ic1];//constant 288.16
				// Ug.tem[ic]  =  288.16;
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

    	    /* upper side */
		    jj = jr -1;  //N+2  ->  N+3,N+4,N+5
		    for(j=jr; j<J0; j++)
		    {
				ic = (i + MyID * config1.ni) * J0 + j;
				ic1 = (i + MyID * config1.ni) * J0 + jj;

		        Ug.q[ic][0] =  Ug.q[ic1][0];
		        Ug.q[ic][1] =  Ug.q[ic1][1];
		        Ug.q[ic][2] =  Ug.q[ic1][2];
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
#endif
}