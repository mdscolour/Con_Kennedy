#include "global.h"
#include "point.h"
#include "walk.h"

int main()
{
	Walk w(1000, "0", 3);
	for (int i = 0; i < 10;i++)
	{
		w.run();
	}


#ifdef _WIN32
	printf("\n");
	system("pause");
#endif
	return 0;
}

