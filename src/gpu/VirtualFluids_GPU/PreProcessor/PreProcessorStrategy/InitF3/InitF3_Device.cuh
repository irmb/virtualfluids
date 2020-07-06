#ifndef LB_INIT_F3_H
#define LB_INIT_F3_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Init_F3(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* rho,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* G6,
	bool EvenOrOdd);

#endif