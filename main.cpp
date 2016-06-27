/**************************************************************************
To start the program using simply in cmd:

$ ./start

This will write down a SAW in the file "FinalWalk".
Maximum 4 input parameters are possible:

$ ./start 1000 "0" 10 1000

1st length of the SAW
2nd the name of the initial walk, "0" means a new walk
3rd the number of pivot per step in MC chain, which will carry out at the same time at the end
4st number of step to go in MC chain
**************************************************************************/
#include "global.h"
#include "point.h"
#include "walk.h"

// Using "Walk w; w.run();" to generate a SAW and write it down in "FinalWalk"
int main(int argc,char *argv[])
{
	int length=100;
	const char* init_name="0";
	int innter_loop = 1000;
	int outer_loop = 10;
	
	if(argc>=2) length=atoi(argv[1]);
	if(argc>=3) init_name=argv[2];
   if(argc>=4) innter_loop=atoi(argv[3]);
	if(argc>=5) outer_loop=atoi(argv[4]);
   
 	Walk w(length, init_name, innter_loop);
 	w.run(outer_loop);

#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}