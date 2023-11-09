#ifndef K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_Device_H
#define K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_Device(	
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
    real* quadricLimiters,
	bool EvenOrOdd);
#endif 