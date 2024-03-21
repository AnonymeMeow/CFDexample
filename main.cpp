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

extern "C"
{
	void grid();
	void readjob();
	void setjob();
	void jobbody();
	void endjob();
}

#include <thread>

int main(int argc, char *argv[])
{
	grid();
	readjob();
	setjob();

	std::thread *threads = new std::thread[nproc];
	for (int i = 0; i < nproc; i++)
	{
		threads[i] = std::thread(
			[&] (int thread_id) {
				MyID = thread_id;
				jobbody();
			},
			i
		);
	}

	for (int i = 0; i < nproc; i++)
	{
		threads[i].join();
	}

	endjob();

    return 0;
}