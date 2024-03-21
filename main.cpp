/*
 *  main.c
 *
 *  A Structure-based Cfd Application for Reactive Flow (SCARF).
 *  Constructed in SWL lab. on Jan 07, 2014
 *  Last updated on Oct. 27, 2014
 */

#define SCARF_DECLSPEC

#include "chemdata.hpp"
#include "comm.hpp"

void grid();
void readjob();
void setjob();
void jobbody();
void endjob();

int main(int argc, char *argv[])
{
	grid();
	readjob();
	setjob();
	jobbody();
	endjob();

    return 0;
}