
/*
 *  main.c
 *
 *  A Structure-based Cfd Application for Reactive Flow (SCARF).
 *  Constructed in SWL lab. on Jan 07, 2014
 *  Last updated on Oct. 27, 2014
 */

int main(int argc, char *argv[])
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

