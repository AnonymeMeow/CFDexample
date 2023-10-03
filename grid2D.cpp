/*
 * retangle_grids.c
 *  Splitted in j direction
 *  Rectangle grids Created on: May 8, 2014
 *  Update for Cylinder grids On Nov.03, 2014
 */

#include<math.h>

#include<iostream>
#include<fstream>

#include"comm.h"

void grid()
{
	/*--------0. read grid configure data --------*/
	std::ifstream file("gridset.dat"); // Gridset configuration file
	if(file.fail())
	{
		std::cout << "gridset.dat file not found" << std::endl;
		exit(0);
	}

	int Ng; // Number of ghost cells
	int nthread; // Number of threads
	file >> config1.ni >> config1.nj >> Ng >> nthread;

	double L_a; // Length of the rectangle in x direction
	double L_b; // Length of the rectangle in y direction
	double alpha1; // The common ratio of the geometric progression grid in x direction
	double alpha2; // The common ratio of the geometric progression grid in y direction
	file >> L_a >> L_b >> alpha1 >> alpha2;

	double xm; // The non-dimensionalized value of L_a
	double ym; // The non-dimensionalized value of L_b
	int gtype; // Type of the grid
	file >> xm >> ym >> gtype;

	double M_1;
	double M_2;
	double R_c;
	if(gtype != 0) //cylinder grid
	{
		file >> M_1 >> M_2 >> R_c;
		file >> alpha1 >> alpha2;
	}

	file.close();

	int ni  = config1.ni * nthread; // Number of cells in x direction
	int I0  = ni + 2*Ng; // Number of cells in x direction, including ghost cells
	int J0  = config1.nj + 2*Ng; // Number of cells in y direction, including ghost cells
	int ir  = ni + Ng; // Number of cells in x direction, including ghost cells at one side
	int jr  = config1.nj + Ng; // Number of cells in y direction, including ghost cells at one side
	int nc  = I0*J0; // Number of all cells, including ghost cells

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
            std::cout << "uniform grid in x direction!" << std::endl; 
        }
        else
        {
        	M_1 = L_a*(1. - alpha1)/(1. - pow(alpha1, ni-1));
        	double Mn_1 = M_1*pow(alpha1, ni-1);
            printf("Geometric progression grid of x: X_1 = %le, X_N=%le\n", M_1, Mn_1);
        }
        if( fabs(alpha2 - 1.0) < 1.e-10 )
        {
            printf("uniform grid in y direction! \n");
        }
        else
        {
        	M_2 = L_b*(1. - alpha2)/(1. - pow(alpha2, config1.nj-1));
        	double Mn_2 = M_2*pow(alpha2, config1.nj-1);
        	printf("Geometric progression grid of y: Y_1 = %le, Y_n=%le\n", M_2, Mn_2);
        }
    	for(int j=Ng; j<jr; j++)
    	{
			double dy;
    		if( fabs(alpha2 - 1.0) < 1.e-10 ) // Uniform grid
    		{
				dy = (double)(j - Ng);
				dy = L_b*dy/config1.nj;
    		}
    		else // Geometric sequence
    		{
    			if(j == Ng)
    				dy = 0.;
    			else
    				dy = M_2*(1. - pow(alpha2, j-Ng))/(1. - alpha2); 				
    		}

    		for(int i=Ng; i<ir; i++)
    		{
    			int ic = i*J0 + j;
				double dx;
    			if( fabs(alpha1 - 1.0) < 1.e-8 ) // Uniform grid
    			{
    				dx = (double)(i - Ng);
    				dx = L_a*dx/ni;
    			}
    			else // Geometric sequence
    			{
    				if(i == Ng)
    					dx = 0.;
    				else
    					dx = M_1*(1. - pow(alpha1, i-Ng))/(1. - alpha1);
    			}
    			mesh.x[ic] = dx;
    			mesh.y[ic] = dy;
            }
    	}
    }
	else
	{
		const double pi = 3.14159265358979;
    	alpha1 = alpha1*(pi/2);
    	alpha2 = alpha2*(pi/2);

    	// Note: one can only create grid points, not cells...
    	// i=Ng line: an ellipse with long-radius L_a, short-radius L_b
    	for(int j=Ng; j<jr; j++)
    	{
    		int ic = Ng*J0 + j;
    		double dy = j - Ng;
    		dy = dy/(config1.nj-1);
    		double rrr = 2.*dy - 1.; // range [-1, 1]
			
			mesh.x[ic] =  L_a * cos(pi - alpha2*rrr) + M_2;
			mesh.y[ic] =  L_b * sin(pi - alpha2*rrr);

    	}
	}

	/*----------2. set the new coordination (xi, et)----------*/

	// uniform mesh
	double d_xi = xm/ni, d_et = ym/config1.nj; // Non-dimensional length of each cell
	for(int j=Ng; j<jr; j++)
	{
		for(int i=Ng; i<ir; i++)
		{
			int ic = i*J0 + j;

			mesh.xi[ic] = (i-Ng)*d_xi;
			mesh.et[ic] = (j-Ng)*d_et;
		}
	}

	/*----------3. Interpolation for ghost cells----------*/
	/* keep the second derivative equals zero
	 * e.g. [x(i+1) - 2x(i) + x(i-1)] = 0 */

	for(int j=Ng; j<jr; j++)
	{
		// left side
	    int ii = Ng - 1; // set ghost cells: 2, 1, 0
	    for(int i=0; i<Ng; i++)
	    {
	    	int ic = ii*J0 + j;
	    	int ic1 = (ii+1)*J0 + j;
	    	int ic2 = (ii+2)*J0 + j;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];

	    	ii -= 1;
	    }

	    // right side
	    for(int i=ir; i<I0; i++) // set ghost cells: N+3, N+4, N+5
	    {
	    	int ic  = i*J0 + j;
	    	int ic1 = (i-1)*J0 + j;
	    	int ic2 = (i-2)*J0 + j;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];
	    }
	}

	/* after interpolate ghost cells in left and right side,
	   interpolate ghost cells for bottom and upper from 0 to I0-1 */
	for(int i=0; i<I0; i++)
	{
	     // bottom side
		int jj = Ng - 1;
	    for(int j=0; j<Ng; j++)
	    {
	    	int ic  = i*J0 + jj;
	    	int ic1 = i*J0 + jj+1;
	    	int ic2 = i*J0 + jj+2;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];

	    	jj -= 1;
	    }

	    // upper side
	    for(int j=jr; j<J0; j++)
	    {
	    	int ic  = i*J0 + j;
	    	int ic1 = i*J0 + j-1;
	    	int ic2 = i*J0 + j-2;

	    	mesh.x[ic]  = 2.*mesh.x[ic1] - mesh.x[ic2];
	    	mesh.y[ic]  = 2.*mesh.y[ic1] - mesh.y[ic2];

	    	mesh.xi[ic] = 2.*mesh.xi[ic1] - mesh.xi[ic2];
	    	mesh.et[ic] = 2.*mesh.et[ic1] - mesh.et[ic2];
	    }
	}

	/*----------4. calculate the derivatives and metric Jacobian----------*/

	// x direction
	int ir1 = I0 -1;
	for(int j=0; j<J0; j++)
	{
		int ic  =  0*J0 + j;
    	int ic1 =  1*J0 + j;
    	mesh.x_xi[ic]  = (mesh.x[ic1] - mesh.x[ic])/d_xi;
    	mesh.y_xi[ic]  = (mesh.y[ic1] - mesh.y[ic])/d_xi;

	    for(int i=1; i<ir1; i++)
	    {
	    	ic  = i*J0 + j;
	    	ic1 = (i-1)*J0 + j;
	    	int ic2 = (i+1)*J0 + j;

	    	mesh.x_xi[ic]  = (mesh.x[ic2] - mesh.x[ic1])/d_xi/2;
	    	mesh.y_xi[ic]  = (mesh.y[ic2] - mesh.y[ic1])/d_xi/2;
	    }

		ic  =  ir1*J0 + j;
    	ic1 =  (ir1-1)*J0 + j;
    	mesh.x_xi[ic]  = (mesh.x[ic] - mesh.x[ic1])/d_xi;
    	mesh.y_xi[ic]  = (mesh.y[ic] - mesh.y[ic1])/d_xi;

	}

	// y direction
	int jr1 = J0 -1;
	for(int i=0; i<I0; i++)
	{
		int ic  =  i*J0 + 0;
    	int ic1 =  i*J0 + 1;
    	mesh.x_et[ic]  = (mesh.x[ic1] - mesh.x[ic])/d_et;
    	mesh.y_et[ic]  = (mesh.y[ic1] - mesh.y[ic])/d_et;

	    for(int j=1; j<jr1; j++)
	    {
	    	ic  = i*J0 + j;
	    	ic1 = i*J0 + j-1;
	    	int ic2 = i*J0 + j+1;

	    	mesh.x_et[ic]  = (mesh.x[ic2] - mesh.x[ic1])/d_et/2;
	    	mesh.y_et[ic]  = (mesh.y[ic2] - mesh.y[ic1])/d_et/2;
	    }

		ic  =  i*J0 + jr1;
    	ic1 =  i*J0 + jr1-1;
    	mesh.x_et[ic]  = (mesh.x[ic] - mesh.x[ic1])/d_et;
    	mesh.y_et[ic]  = (mesh.y[ic] - mesh.y[ic1])/d_et;
	}

	for(int i=0; i<I0; i++)
	{
		for(int j=0; j<J0; j++)
		{
			int ic = i*J0 + j;
			mesh.yaks[ic] = mesh.x_xi[ic]*mesh.y_et[ic] - mesh.y_xi[ic]*mesh.x_et[ic];
		}
	}
    std::cout << "complete!!!" << std::endl;
}