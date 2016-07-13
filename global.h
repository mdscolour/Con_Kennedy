/**************************************************************************
Global definition everywhere, mainly macro and library
and will:
// Define "RNG_NAME" for random number, in linux using drand48(), in windows using rand(),
// the initialization of the seed is in the constructor of class Walk
**************************************************************************/
#pragma once
#define _CRT_SECURE_NO_DEPRECATE
#define M_PI       3.14159265358979323846
#include <stdlib.h>
#include<math.h>
#include<stdio.h>
#include <iostream>
#include <vector>
#include<fstream>
#include<string>
#include<iterator>
#include <sstream>
#include <time.h>
using namespace std;

// Define "RNG_NAME" for random number, in linux using drand48(), in windows using rand(),
// the initialization of the seed is in the constructor of class Walk
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

// Other parameters
//#define MODEL_NAME Sphere
#define RADIUS 0.4
#define RIGID 20
#define OPERATION_NAME OpMatrix
