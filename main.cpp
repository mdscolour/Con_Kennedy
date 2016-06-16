#include "global.h"
#include "point.h"
#include "walk.h"

// Using "Walk w; w.run();" to generate a SAW and write it down in "FinalWalk"
int main(int argc,char *argv[])
{
	int length=100;
	char* init_name="0";
	int MCStep = 10;
	int stepnum = 1;
	
	if(argc>=2) length=atoi(argv[1]);
	if(argc>=3) init_name=argv[2];
    if(argc>=4) MCStep=atoi(argv[3]);
	if(argc>=5) stepnum=atoi(argv[4]);
   
	Walk w(length, init_name, MCStep);
	for(int i=0;i<stepnum;i++)
	{
	      w.run();
	}
#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}