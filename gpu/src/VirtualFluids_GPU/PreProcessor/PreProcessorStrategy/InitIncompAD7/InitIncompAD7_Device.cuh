#ifndef LB_INIT_INCOMP_AD7_H
#define LB_INIT_INCOMP_AD7_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Init_Incomp_AD_7(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* Conc,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD7,
	bool EvenOrOdd);

#endif