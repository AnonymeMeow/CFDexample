/*
 *  main.c
 *
 *  A Structure-based Cfd Application for Reactive Flow (SCARF).
 *  Constructed in SWL lab. on Jan 07, 2014
 *  Last updated on Oct. 27, 2014
 */

#include "comm.h"

int main(int argc, char* argv[])
{
	grid();
	readjob();
	setjob();
	jobbody();
	endjob();
	return 0;
}