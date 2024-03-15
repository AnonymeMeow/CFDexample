/*
 * startjob.c
 *
 *  Created on: 2014.01.09
 *
 */

#include<string.h>
#include<mpi.h>
#include"comm.h"


/*---------------------------------------------------
 * Start the simulation. Write output information
 * ------------------------------------------------*/
void startjob(int argc, char *argv[])
{

	FILE *outId;

    // MPI_Init(&argc,&argv);
	// MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	// MPI_Comm_rank(MPI_COMM_WORLD, &MyID);

	/*------write the ouput file-------*/
	// if(MyID == 0)
	// {
		printf("simulation start...\n");
		outId = fopen("outInfo.dat", "w");
		fprintf(outId, "/***************************** simulation output *****************************/ \n");
		fprintf(outId, " \n");
		fclose(outId);
	// }
}