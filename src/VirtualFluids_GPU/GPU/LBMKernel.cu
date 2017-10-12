#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <helper_cuda.h>

#include "LBM/LB.h"
#include "GPU/GPU_Kernels.cuh"

#include "GPU_Interface.h"


void KernelBGKSPSimple27(unsigned int numberOfThreads, 
                                    doubflo s9,
                                    unsigned int* bcMatD,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    doubflo* DD,
                                    int size_Mat,
                                    bool EvenOrOdd)
{
   int Grid = (size_Mat / numberOfThreads) + 1;
   int Grid1, Grid2;
   if (Grid > 512)
   {
      Grid1 = 512;
      Grid2 = (Grid / Grid1) + 1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 grid(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      LB_Kernel_BGK_SP_Simple_27<<< grid, threads >>>(  s9,
                                                         bcMatD,
                                                         neighborX,
                                                         neighborY,
                                                         neighborZ,
                                                         DD,
                                                         size_Mat,
                                                         EvenOrOdd); 
      getLastCudaError("LB_Kernel_BGK_SP_Simple_27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void InitSP27(   unsigned int numberOfThreads,
                            unsigned int* neighborX,
                            unsigned int* neighborY,
                            unsigned int* neighborZ,
                            unsigned int* geoD,
                            doubflo* rho,
                            doubflo* ux,
                            doubflo* uy,
                            doubflo* uz,
                            unsigned int size_Mat,
                            doubflo* DD,
                            bool EvenOrOdd)
{
   int Grid = (size_Mat / numberOfThreads)+1;
   int Grid1, Grid2;
   if (Grid>512)
   {
      Grid1 = 512;
      Grid2 = (Grid/Grid1)+1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 grid(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      LBInitSP27<<< grid, threads >>>( neighborX,
                                       neighborY,
                                       neighborZ,
                                       geoD,
                                       rho,
                                       ux,
                                       uy,
                                       uz,
                                       size_Mat,
                                       DD,
                                       EvenOrOdd);
      getLastCudaError("LBInitSP27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void CalcMacCompSP27( doubflo* vxD,
								 doubflo* vyD,
								 doubflo* vzD,
								 doubflo* rhoD,
								 doubflo* pressD,
								 unsigned int* geoD,
								 unsigned int* neighborX,
								 unsigned int* neighborY,
								 unsigned int* neighborZ,
								 unsigned int size_Mat,
								 unsigned int numberOfThreads, 
								 doubflo* DD,
								 bool evenOrOdd)
{ 
   int Grid = (size_Mat / numberOfThreads)+1;
   int Grid1, Grid2;
   if (Grid>512)
   {
      Grid1 = 512;
      Grid2 = (Grid/Grid1)+1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 grid(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      LBCalcMacCompSP27<<< grid, threads >>> (   vxD, 
												 vyD, 
												 vzD, 
												 rhoD, 
												 pressD, 
												 geoD, 
												 neighborX,
												 neighborY,
												 neighborZ,
												 size_Mat, 
												 DD, 
												 evenOrOdd); 
      getLastCudaError("LBCalcMacCompSP27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void QDevComp27( unsigned int numberOfThreads,
							int nx,
							int ny,
							doubflo* DD, 
							int* k_Q, 
							doubflo* QQ,
							unsigned int sizeQ,
							unsigned int kQ, 
							doubflo om1, 
							unsigned int* neighborX,
							unsigned int* neighborY,
							unsigned int* neighborZ,
							unsigned int size_Mat, 
							bool evenOrOdd)
{
   int Grid = (kQ / numberOfThreads)+1;
   int Grid1, Grid2;
   if (Grid>512)
   {
      Grid1 = 512;
      Grid2 = (Grid/Grid1)+1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      QDeviceComp27<<< gridQ, threads >>> (nx,
										   ny,
										   DD, 
										   k_Q, 
										   QQ,
										   sizeQ,
										   kQ, 
										   om1, 
										   neighborX,
										   neighborY,
										   neighborZ,
										   size_Mat, 
										   evenOrOdd);
      getLastCudaError("QDeviceComp27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void QVelDevComp27(unsigned int numberOfThreads,
							  int nx,
							  int ny,
							  doubflo* vx,
							  doubflo* vy,
							  doubflo* vz,
							  doubflo* DD, 
							  int* k_Q, 
							  doubflo* QQ,
							  unsigned int sizeQ,
							  unsigned int kQ, 
							  doubflo om1, 
							  unsigned int* neighborX,
							  unsigned int* neighborY,
							  unsigned int* neighborZ,
							  unsigned int size_Mat, 
							  bool evenOrOdd)
{
   int Grid = (kQ / numberOfThreads)+1;
   int Grid1, Grid2;
   if (Grid>512)
   {
      Grid1 = 512;
      Grid2 = (Grid/Grid1)+1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      QVelDeviceComp27<<< gridQ, threads >>> (nx,
											  ny,
											  vx,
											  vy,
											  vz,
											  DD, 
											  k_Q, 
											  QQ,
											  sizeQ,
											  kQ, 
											  om1, 
											  neighborX,
											  neighborY,
											  neighborZ,
											  size_Mat, 
											  evenOrOdd);
      getLastCudaError("QVelDeviceComp27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void QPressDevOld27(  unsigned int numberOfThreads,
                                     doubflo* rhoBC,
                                     doubflo* DD, 
                                     int* k_Q, 
                                     int* k_N, 
                                     unsigned int kQ, 
                                     doubflo om1, 
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int size_Mat, 
                                     bool evenOrOdd)
{
   int Grid = (kQ / numberOfThreads)+1;
   int Grid1, Grid2;
   if (Grid>512)
   {
      Grid1 = 512;
      Grid2 = (Grid/Grid1)+1;
   } 
   else
   {
      Grid1 = 1;
      Grid2 = Grid;
   }
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      QPressDeviceOld27<<< gridQ, threads >>> ( rhoBC,
                                                DD, 
                                                k_Q, 
                                                k_N, 
                                                kQ, 
                                                om1, 
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                size_Mat, 
                                                evenOrOdd);
      getLastCudaError("QPressDeviceOld27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void GetSendFsPreDev27(doubflo* DD,
								  doubflo* bufferFs,
								  int* sendIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat, 
								  bool evenOrOdd,
								  unsigned int numberOfThreads)
{
	int Grid = (buffmax / numberOfThreads)+1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid/Grid1)+1;
	} 
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1 );

	getSendFsPre27<<< grid, threads >>>(DD, 
										bufferFs, 
										sendIndex, 
										buffmax,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat, 
										evenOrOdd);
	getLastCudaError("getSendFsPre27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void GetSendFsPostDev27(doubflo* DD,
								   doubflo* bufferFs,
								   int* sendIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat, 
								   bool evenOrOdd,
								   unsigned int numberOfThreads)
{
	int Grid = (buffmax / numberOfThreads)+1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid/Grid1)+1;
	} 
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1 );

	getSendFsPost27<<< grid, threads >>>(DD, 
										 bufferFs, 
										 sendIndex, 
										 buffmax,
										 neighborX,
										 neighborY,
										 neighborZ,
										 size_Mat, 
										 evenOrOdd);
	getLastCudaError("getSendFsPost27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void SetRecvFsPreDev27(doubflo* DD,
								  doubflo* bufferFs,
								  int* recvIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat, 
								  bool evenOrOdd,
								  unsigned int numberOfThreads)
{
	int Grid = (buffmax / numberOfThreads)+1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid/Grid1)+1;
	} 
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1 );

	setRecvFsPre27<<< grid, threads >>>(DD, 
										bufferFs, 
										recvIndex, 
										buffmax,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat, 
										evenOrOdd);
	getLastCudaError("setRecvFsPre27 execution failed"); 
}
//////////////////////////////////////////////////////////////////////////
void SetRecvFsPostDev27(doubflo* DD,
								   doubflo* bufferFs,
								   int* recvIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat, 
								   bool evenOrOdd,
								   unsigned int numberOfThreads)
{
	int Grid = (buffmax / numberOfThreads)+1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid/Grid1)+1;
	} 
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1 );

	setRecvFsPost27<<< grid, threads >>>(DD, 
										 bufferFs, 
										 recvIndex, 
										 buffmax,
										 neighborX,
										 neighborY,
										 neighborZ,
										 size_Mat, 
										 evenOrOdd);
	getLastCudaError("setRecvFsPost27 execution failed"); 
}
