#ifndef LB_INIT_SP27_H
#define LB_INIT_SP27_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Init_SP_27(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* rho,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD,
	bool EvenOrOdd);

#endif