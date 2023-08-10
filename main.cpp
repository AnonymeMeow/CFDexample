
/*
 *  main.c
 *
 *  A Structure-based Cfd Application for Reactive Flow (SCARF).
 *  Constructed in SWL lab. on Jan 07, 2014
 *  Last updated on Oct. 27, 2014
 */

int scarf_main(int argc, char *argv[])
{
	void startjob(int argc, char *argv[]);
	void readjob();
	void setjob();
	void jobbody();
	void endjob();

	startjob(argc, argv);

	readjob();

	setjob();

	jobbody();

    endjob();

    return 0;
}

#include "comm.h"
#include "chemdata.h"

int main(int argc, char* argv[])
{
	int grid_main(int argc, char* argv[]);

	if (argc == 1)
	{
		char cmd[512];
		sprintf(cmd, "mpiexec -n %d \"%s\" extra_arg", grid_main(argc, argv), *argv);
		return system(cmd);
	}
	return scarf_main(argc, argv);
}