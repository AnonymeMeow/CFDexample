#include<iostream>
#include<fstream>

#include"comm.h"

void readjob()
{
	std::ifstream gridset_file("gridset.dat");
	if(gridset_file.fail())
	{
		std::cout << "gridset.dat file not found" << std::endl;
		exit(0);
	}

	double xm; // The non-dimensionalized value of L_a
	double ym; // The non-dimensionalized value of L_b
	gridset_file
		>> config1.ni >> config1.nj >> config1.Ng >> config1.nblock
		>> config2.Lx >> config2.Ly
		>> xm >> ym;

	gridset_file.close();


	std::ifstream config_file("config.dat");
	if(config_file.fail())
	{
		std::cout << "config.dat file not found!" << std::endl;
		exit(0);
	}

	config_file
		>> config2.t0 >> config2.x0
		>> config1.newrun >> config1.nonDi >> config1.useDt
		>> config1.iStep0 >> config1.nStep >> config1.nRamp >> config1.Samples >> config1.ifilm
		>> config2.dt0 >> config2.dt1 >> config2.CFL0 >> config2.CFL1
		>> config1.gasModel >> config1.reacModel >> config1.visModel
		>> config2.molWeight >> config2.gam0 >> config2.Pr0
		>> config2.muRef >> config2.suthC1 >> config2.suthC2
		>> config2.MaRef >> config2.temRef >> config2.preRef
		>> config2.p1 >> config2.T1 >> config2.u1 >> config2.v1
		>> config2.p2 >> config2.T2 >> config2.u2 >> config2.v2;

	config_file.close();

	std::ofstream log_file("outInfo.dat");

	log_file
		<< "\n/---------------------------configure data---------------------------/\n"

		<< "\nt0=" << config2.t0 << ", x0=" << config2.x0 << ", Lx=" << config2.Lx << ", Ly=" << config2.Ly << "\n"

		<< "\ngrids point: ni=" << config1.ni << ", nj=" << config1.nj << ", ghost cells Ng=" << config1.Ng << ", nblock=" << config1.nblock << "\n"

		<< "\nnewrun=" << config1.newrun << ", nonDi=" << config1.nonDi << ", useDt=" << config1.useDt << " \n"
		<< "iStep0=" << config1.iStep0 << ", nStep=" << config1.nStep << ", nRamp=" << config1.nRamp << ", Samples=" << config1.Samples << ", ifilm=" << config1.ifilm << " \n"
		<< "dt0=" << config2.dt0 << ", dt1=" << config2.dt1 << ", CFL0=" << config2.CFL0 << ", CFL1=" << config2.CFL1 << " \n"

		<< "\ngasModel=" << config1.gasModel << ", reacModel=" << config1.reacModel << ", visModel=" << config1.visModel << " \n"
		<< "molWeight=" << config2.molWeight << ", gam0=" << config2.gam0 << ", pr0=" << config2.Pr0 << "\n"
		<< "mu0=" << config2.muRef << ", suthC1=" << config2.suthC1 << ", suthC2=" << config2.suthC2 << "\n"
		<< "Ma0=" << config2.MaRef << ", temRef=" << config2.temRef << ", preRef=" << config2.preRef << "\n"

		<< "\nInitial condition: \n"
		<< "p_1=" << config2.p1 << ", T_1=" << config2.T1 << ", u_1=" << config2.u1 << ", v1=" << config2.v1 << " \n"
		<< "p_2=" << config2.p2 << ", T_2=" << config2.T2 << ", u_2=" << config2.u2 << ", v2=" << config2.v2 << " \n";

	log_file.close();
		
	neqv = 4; // solution variables without chemical terms
	config1.timeOrder = 3;

	/*---1. Grid information---*/
	I0 = config1.ni + 2*config1.Ng;
	J0 = config1.nj + 2*config1.Ng;

	int ir  = config1.ni + config1.Ng; // Number of cells in x direction, including ghost cells at one side
	int jr  = config1.nj + config1.Ng; // Number of cells in y direction, including ghost cells at one side
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

	/*----------------------1.1 Create a Rectangle mesh (x, y)----------------------*/
	for(int j=config1.Ng; j<jr; j++)
	{
		double dy;

		dy = (double)(j - config1.Ng);
		dy = config2.Ly*dy/config1.nj;

		for(int i=config1.Ng; i<ir; i++)
		{
			int ic = i*J0 + j;
			double dx;

			dx = (double)(i - config1.Ng);
			dx = config2.Lx*dx/config1.ni;

			mesh.x[ic] = dx;
			mesh.y[ic] = dy;
		}
	}

	/*----------2. set the new coordination (xi, et)----------*/

	// uniform mesh
	double d_xi = xm/config1.ni, d_et = ym/config1.nj; // Non-dimensional length of each cell
	for(int j=config1.Ng; j<jr; j++)
	{
		for(int i=config1.Ng; i<ir; i++)
		{
			int ic = i*J0 + j;

			mesh.xi[ic] = (i-config1.Ng)*d_xi;
			mesh.et[ic] = (j-config1.Ng)*d_et;
		}
	}

	/*----------3. Interpolation for ghost cells----------*/
	/* keep the second derivative equals zero
	 * e.g. [x(i+1) - 2x(i) + x(i-1)] = 0 */

	for(int j=config1.Ng; j<jr; j++)
	{
		// left side
	    int ii = config1.Ng - 1; // set ghost cells: 2, 1, 0
	    for(int i=0; i<config1.Ng; i++)
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
		int jj = config1.Ng - 1;
	    for(int j=0; j<config1.Ng; j++)
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
		
	/* After coordinate transformation,
	 * the distance between cells are equal */
	dxc = mesh.xi[1*J0] - mesh.xi[0*J0];
	dyc = mesh.et[1]    - mesh.et[0];
}