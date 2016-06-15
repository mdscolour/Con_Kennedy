#pragma once
#define _CRT_SECURE_NO_DEPRECATE
#define M_PI       3.14159265358979323846
#include <stdlib.h>
#include<math.h>
#include<stdio.h>
#include <iostream>
using namespace std;

#ifdef linux
#include <sys/time.h>
#define RNG_NAME drand48
#endif
#ifdef _WIN32
#include <time.h>
inline double GetRand()
{
	//srand(time(NULL));
	return (rand() % 10000) / 10000.0; // ".0"
}
#define RNG_NAME GetRand
#endif


#define MODELNAME Sphere
#define RADIUS 0.4
#define RADIANNUM 100
//#define NUM_SYM 47
#define NUM_SYM 119*RADIANNUM*RADIANNUM*RADIANNUM
#define ROTNUM 15  // number of the rotation matrix 0-14