/*
 * retangle_grids.c
 *  Splitted in j direction
 *  Rectangle grids Created on: May 8, 2014
 *  Update for Cylinder grids On Nov.03, 2014
 */

#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "comm.h"

void grid()
{
	int i, j, ic, ic1, ic2, nc,id, ir1, ii, jj, gtype,
		jr, jr1, i1, ik, ni, ir;
	double xm, ym, dx, dy, d_xi, d_et, R_c, rrr, pi, jm1, sumj1;
	double M_1, M_2, Mn_1, Mn_2, alpha1, alpha2;
	char filename[20];
	FILE  *fp;

	/*--------0. read grid configure data --------*/
	fp    = fopen("gridset.dat", "r");
	if(fp == NULL)
	{
		printf("gridset.dat file not found \n");
		exit(0);
	}
	fscanf(fp, "%d %d %d %d", &config1.ni, &config1.nj, &config1.Ng, &nproc);
	fscanf(fp, "%lf %lf %lf %lf", &config2.Lx, &config2.Ly, &alpha1, &alpha2);
	fscanf(fp, "%lf %lf %d", &xm, &ym, &gtype);
	if(gtype != 0)
	{
		//cylinder grid
		fscanf(fp, "%lf %lf %lf", &M_1, &M_2, &R_c);
		fscanf(fp, "%lf %lf", &alpha1, &alpha2);
	}

	fclose(fp);

	ni  = config1.ni * nproc;
	I0  = ni + 2*config1.Ng;
	J0  = config1.nj + 2*config1.Ng;
	ir  = ni + config1.Ng;
	jr  = config1.nj + config1.Ng;
	nc  = I0*J0;

	mesh.x     = (double*)malloc(sizeof(double)*nc);
	mesh.y     = (double*)malloc(sizeof(double)*nc);
	mesh.xi    = (double*)malloc(sizeof(double)*nc);
	mesh.et    = (double*)malloc(sizeof(double)*nc);
	mesh.x_xi  = (double*)malloc(sizeof(double)*nc);
	mesh.x_et  = (double*)malloc(sizeof(double)*nc);
	mesh.y_xi  = (double*)malloc(sizeof(double)*nc);
	mesh.y_et  = (double*)malloc(sizeof(double)*nc);
	mesh.yaks  = (double*)malloc(sizeof(double)*nc);


    if(gtype == 0)
    {
    	/*----------------------1.1 Create a Rectangle mesh (x, y)----------------------*/

        if( fabs(alpha1 - 1.0) < 1.e-8 )
        {
            printf("uniform grid in x direction! \n");
        }
        else
        {
        	M_1 = config2.Lx*(1. - alpha1)/(1. - pow(alpha1, ni-1));
        	Mn_1 = M_1*pow(alpha1, ni-1);
            printf("Geometric progression grid of x: X_1 = %le, X_N=%le\n", M_1, Mn_1);
        }
        if( fabs(alpha2 - 1.0) < 1.e-10 )
        {
            printf("uniform grid in y direction! \n");
        }
        else
        {
        	M_2 = config2.Ly*(1. - alpha2)/(1. - pow(alpha2, config1.nj-1));
        	Mn_2 = M_2*pow(alpha2, config1.nj-1);
        	printf("Geometric progression grid of y: Y_1 = %le, Y_n=%le\n", M_2, Mn_2);
        }
    	for(j=config1.Ng; j<jr; j++)
    	{
    		if( fabs(alpha2 - 1.0) < 1.e-10 )
    		{
    			// Uniform grid
				dy = (double)(j - config1.Ng);
				dy = config2.Ly*dy/config1.nj;
    		}
    		else
    		{
    			// Geometric sequence
#if defined Cavity || defined Plate
    			if(j == config1.Ng)
    				dy = 0.;
    			else
    				dy = M_2*(1. - pow(alpha2, j-config1.Ng))/(1. - alpha2);
#elif defined Tube
    			sumj1 = (j - config1.Ng)*1.e-4;
    			if(sumj1 <= 0.01)
    			{
    				dy = (double)(j - config1.Ng);
    				dy = 1.e-4*dy;
    				sumj1 = dy;
    				jm1 = j;
    			}
    			else
    			{
    				jj = j - jm1;
    				dy = sumj1 + 1.e-4*jj+ (jj-1)*jj/2.*alpha2;
    			}
#endif
    		}

    		for(i=config1.Ng; i<ir; i++)
    		{
    			ic = i*J0 + j;
    			if( fabs(alpha1 - 1.0) < 1.e-8 )
    			{
    				// Uniform grid
    				dx = (double)(i - config1.Ng);
    				dx = config2.Lx*dx/ni;
    			}
    			else
    			{
    				if(i == config1.Ng)
    					dx = 0.;
    				else
    					dx = M_1*(1. - pow(alpha1, i-config1.Ng))/(1. - alpha1);
    			}
    			mesh.x[ic] = dx;
    			mesh.y[ic] = dy;
            }
    	}
    }
	else
	{
		pi = 3.1415926;
    	alpha1 = alpha1*(pi/2);
    	alpha2 = alpha2*(pi/2);

    	// Note: one can only create grid points, not cells...
    	// i=config1.Ng line: an ellipse with long-radius config2.Lx, short-radius config2.Ly
    	for(j=config1.Ng; j<jr; j++)
    	{
    		ic = config1.Ng*J0 + j;
    		dy = j - config1.Ng;
    		dy = dy/(config1.nj-1);
    		rrr = 2.*dy - 1.; // range [-1, 1]
			
			mesh.x[ic] =  config2.Lx * cos(pi - alpha2*rrr) + M_2;
			mesh.y[ic] =  config2.Ly * sin(pi - alpha2*rrr);
    	}
    	// i=ni+config1.Ng line: a circle with radius R_c, center M_1
    	// for(j=config1.Ng; j<jr; j++)
    	// {
    	// 	ic = (ir-1)*J0 + j;
    	// 	dy = j - config1.Ng;
    	// 	dy = dy/(config1.nj-1);
    	// 	rrr = 2.*dy - 1.; // range [-1, 1]
    	// 	mesh.x[ic] = R_c * cos(pi - alpha1*rrr) + M_1;
	    // 	mesh.y[ic] = R_c * sin(pi - alpha1*rrr);
    	// }
	}
    // else
    // {

    // 	/*----------------------1.2 Create a Cylinder mesh (x, y)----------------------*/
    // 	pi = 3.1415926;
    // 	alpha1 = alpha1*(pi/2);
    // 	alpha2 = alpha2*(pi/2);

    // 	// Note: one can only create grid points, not cells...
    // 	// i=config1.Ng line: an ellipse with long-radius config2.Lx, short-radius config2.Ly
    // 	for(j=config1.Ng; j<jr; j++)
    // 	{
    // 		ic = config1.Ng*J0 + j;
    // 		dy = j - config1.Ng;
    // 		dy = dy/(config1.nj-1);
    // 		rrr = 2.*dy - 1.; // range [-1, 1]
	// 		mesh.x[ic] =  config2.Lx * cos(pi - alpha2*rrr) + M_2;
	// 		mesh.y[ic] =  config2.Ly * sin(pi - alpha2*rrr);

    // 	}
    // 	// i=ni+config1.Ng line: a circle with radius R_c, center M_1
    // 	for(j=config1.Ng; j<jr; j++)
    // 	{
    // 		ic = (ir-1)*J0 + j;
    // 		dy = j - config1.Ng;
    // 		dy = dy/(config1.nj-1);
    // 		rrr = 2.*dy - 1.; // range [-1, 1]
    // 		mesh.x[ic] = R_c * cos(pi - alpha1*rrr) + M_1;
	//     	mesh.y[ic] = R_c * sin(pi - alpha1*rrr);
    // 	}
	// }

	// Interpolation
    /*
	for(j=config1.Ng; j<jr; j++)
	{
		ic1 = config1.Ng*J0     + j;
		ic2 = (ir-1)*J0 + j;
		for(i=config1.Ng+1; i<ir-1; i++)
		{
			ic = i*J0 + j;
	    		dx = (i - config1.Ng);
	    		dx = dx/(ni-1);
	    		mesh.x[ic] = mesh.x[ic1] + dx*(mesh.x[ic2] - mesh.x[ic1]);
	    		mesh.y[ic] = mesh.y[ic1] + dx*(mesh.y[ic2] - mesh.y[ic1]);
		}
	} */

	/*----------2. set the new coordination (xi, et)----------*/

	d_xi = xm/ni;
	d_et = ym/config1.nj; // uniform meshmesh.x[ic];mesh.x[ic];
	for(j=config1.Ng; j<jr; j++)
	{
		for(i=config1.Ng; i<ir; i++)
		{
			ic = i*J0 + j;

			mesh.xi[ic] = (i-config1.Ng)*d_xi;
			mesh.et[ic] = (j-config1.Ng)*d_et;
		}
	}

	/*----------3. Interpolation for ghost cells----------*/
	/* keep the second derivative equals zero
	 * e.g. [x(i+1) - 2x(i) + x(i-1)] = 0 */

	for(j=config1.Ng; j<jr; j++)
	{
		// left side
	    ii = config1.Ng - 1; // set ghost cells: 2, 1, 0
	    for(i=0; i<config1.Ng; i++)
	    {
	    	ic = ii*J0 + j;
	    	ic1 = (ii+1)*J0 + j;
	    	ic2 = (ii+2)*J0 + j;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];

	    	ii -= 1;
	    }

	    // right side
	    for(i=ir; i<I0; i++) // set ghost cells: N+3, N+4, N+5
	    {
	    	ic  = i*J0 + j;
	    	ic1 = (i-1)*J0 + j;
	    	ic2 = (i-2)*J0 + j;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];
	    }
	}

	/* after interpolate ghost cells in left and right side,
	   interpolate ghost cells for bottom and upper from 0 to I0-1 */
	for(i=0; i<I0; i++)
	{
		// bottom side
		jj = config1.Ng - 1;
	    for(j=0; j<config1.Ng; j++)
	    {
	    	ic  = i*J0 + jj;
	    	ic1 = i*J0 + jj+1;
	    	ic2 = i*J0 + jj+2;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];

	    	jj -= 1;
	    }

	    // upper side
	    for(j=jr; j<J0; j++)
	    {
	    	ic  = i*J0 + j;
	    	ic1 = i*J0 + j-1;
	    	ic2 = i*J0 + j-2;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];
	    }
	}

	/*----------4. calculate the derivatives and metric Jacobian----------*/

	// x direction
	ir1 = I0 -1;
	for(j=0; j<J0; j++)
	{
		ic  =  0*J0 + j;
    	ic1 =  1*J0 + j;
    	mesh.x_xi[ic]  = (mesh.x[ic1] - mesh.x[ic])/d_xi;
    	mesh.y_xi[ic]  = (mesh.y[ic1] - mesh.y[ic])/d_xi;

	    for(i=1; i<ir1; i++)
	    {
	    	ic  = i*J0 + j;
	    	ic1 = (i-1)*J0 + j;
	    	ic2 = (i+1)*J0 + j;

	    	mesh.x_xi[ic]  = (mesh.x[ic2] - mesh.x[ic1])/d_xi/2;
	    	mesh.y_xi[ic]  = (mesh.y[ic2] - mesh.y[ic1])/d_xi/2;
	    }

		ic  =  ir1*J0 + j;
    	ic1 =  (ir1-1)*J0 + j;
    	mesh.x_xi[ic]  = (mesh.x[ic] - mesh.x[ic1])/d_xi;
    	mesh.y_xi[ic]  = (mesh.y[ic] - mesh.y[ic1])/d_xi;
	}

	// y direction
	jr1 = J0 -1;
	for(i=0; i<I0; i++)
	{
		ic  =  i*J0 + 0;
    	ic1 =  i*J0 + 1;
    	mesh.x_et[ic]  = (mesh.x[ic1] - mesh.x[ic])/d_et;
    	mesh.y_et[ic]  = (mesh.y[ic1] - mesh.y[ic])/d_et;

	    for(j=1; j<jr1; j++)
	    {
	    	ic  = i*J0 + j;
	    	ic1 = i*J0 + j-1;
	    	ic2 = i*J0 + j+1;

	    	mesh.x_et[ic]  = (mesh.x[ic2] - mesh.x[ic1])/d_et/2;
	    	mesh.y_et[ic]  = (mesh.y[ic2] - mesh.y[ic1])/d_et/2;
	    }

		ic  =  i*J0 + jr1;
    	ic1 =  i*J0 + jr1-1;
    	mesh.x_et[ic]  = (mesh.x[ic] - mesh.x[ic1])/d_et;
    	mesh.y_et[ic]  = (mesh.y[ic] - mesh.y[ic1])/d_et;
	}

	for(i=0; i<I0; i++)
	{
		for(j=0; j<J0; j++)
		{
			ic = i*J0 + j;
			mesh.yaks[ic] = mesh.x_xi[ic]*mesh.y_et[ic] -
					        mesh.y_xi[ic]*mesh.x_et[ic];
		}
	}
    printf("complete!!! \n");
}