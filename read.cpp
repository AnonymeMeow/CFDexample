/*
 * read.c
 *
 *  Created on: 2014.1.14.
 *
 */

#include <ctype.h>
#include <string.h>
#include <math.h>

#include <fstream>
#include <iostream>

#include "comm.h"

/*---------------------------------------------------
 * read the necessary information for simulation
 * ------------------------------------------------*/
void readjob()
{
	void readconfig();

	neqv = 4; // solution variables without chemical terms

	if(MyID == 0)
	{
		readconfig();
	}
}

/*---------------------------------------------------
 * read configuration information
 * ------------------------------------------------*/
void readconfig()
{
	std::ifstream config_file("config.dat");
	std::ofstream log_file("outInfo.dat");
	if(config_file.fail())
	{
		std::cout << "config.dat file not found!" << std::endl;
		exit(0);
	}

	config_file >> config2.t0 >> config2.x0
				>> config1.newrun >> config1.nonDi >> config1.useDt >> config1.iStep0 >> config1.nStep >> config1.nRamp >> config1.Samples >> config1.ifilm
				>> config2.dt0 >> config2.dt1 >> config2.CFL0 >> config2.CFL1
				>> config1.gasModel >> config1.reacModel >> config1.visModel
				>> config2.molWeight >> config2.gam0 >> config2.Pr0
				>> config2.muRef >> config2.suthC1 >> config2.suthC2
				>> config2.MaRef >> config2.temRef >> config2.preRef
				>> config2.p1 >> config2.T1 >> config2.u1 >> config2.v1
				>> config2.p2 >> config2.T2 >> config2.u2 >> config2.v2;
	config_file.close();

	std::ifstream gridset_file("gridset.dat");
	if(gridset_file.fail())
	{
		std::cout << "gridset.dat file not found" << std::endl;
		exit(0);
	}

	gridset_file >> config1.ni >> config1.nj >> config1.Ng >> config1.nblock
				 >> config2.Lx >> config2.Ly;
	gridset_file.close();

	log_file << "\n/---------------------------configure data---------------------------/\n"
			 
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

	config1.timeOrder = 3;

	log_file.close();

	/*---1. Grid information---*/
	I0 = config1.ni + 2*config1.Ng;
	J0 = config1.nj + 2*config1.Ng;
	
	/* After coordinate transformation,
	 * the distance between cells are equal */
	dxc = mesh.xi[1*J0] - mesh.xi[0*J0];
	dyc = mesh.et[1]    - mesh.et[0];

	int ns, ib, nb;

	nb = 2; // No. of initial condition

	inc[0].p = config2.p1;
	inc[0].t = config2.T1;
	inc[0].u = config2.u1;
	inc[0].v = config2.v1;
	inc[1].p = config2.p2;
	inc[1].t = config2.T2;
	inc[1].u = config2.u2;
	inc[1].v = config2.v2;
}