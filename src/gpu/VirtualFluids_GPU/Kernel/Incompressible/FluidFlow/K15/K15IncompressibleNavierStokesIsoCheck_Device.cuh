#ifndef K15IncompressibleNavierStokesIsoCheck_Device_27
#define K15IncompressibleNavierStokesIsoCheck_Device_27

#include <DataTypes.h>
#include <curand.h>

__global__ void K15IncompressibleNavierStokesIsoCheck_Device(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* dxxUx,
	real* dyyUy,
	real* dzzUz,
	int size_Mat,
	bool EvenOrOdd);

#endif