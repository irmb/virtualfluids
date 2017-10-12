#ifndef D3Q27_KERNELS_H
#define D3Q27_KERNELS_H

#include "cuda_runtime.h"
#include "LBM/LB.h"

__global__ void LB_Kernel_BGK_SP_Simple_27(doubflo s9,
                                                      unsigned int* bcMatD,
                                                      unsigned int* neighborX,
                                                      unsigned int* neighborY,
                                                      unsigned int* neighborZ,
                                                      doubflo* DDStart,
                                                      int size_Mat,
                                                      bool EvenOrOdd);

__global__ void LBInitSP27( unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int* geoD,
                                       doubflo* rho,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       unsigned int size_Mat,
                                       doubflo* DD,
                                       bool EvenOrOdd);

__global__ void LBCalcMacCompSP27( doubflo* vxD,
											  doubflo* vyD,
											  doubflo* vzD,
											  doubflo* rhoD,
											  doubflo* pressD,
											  unsigned int* geoD,
											  unsigned int* neighborX,
											  unsigned int* neighborY,
											  unsigned int* neighborZ,
											  unsigned int size_Mat,
											  doubflo* DD,
											  bool evenOrOdd);

__global__ void QDeviceComp27(int inx,
										 int iny,
										 doubflo* DD, 
										 int* k_Q, 
										 doubflo* QQ,
										 unsigned int sizeQ,
										 int kQ, 
										 doubflo minusomega, 
										 unsigned int* neighborX,
										 unsigned int* neighborY,
										 unsigned int* neighborZ,
										 unsigned int size_Mat, 
										 bool evenOrOdd);

__global__ void QVelDeviceComp27(int inx,
											int iny,
											doubflo* vx,
											doubflo* vy,
											doubflo* vz,
											doubflo* DD, 
											int* k_Q, 
											doubflo* QQ,
											unsigned int sizeQ,
											int kQ, 
											doubflo om1, 
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
											bool evenOrOdd);

__global__ void QPressDeviceOld27(doubflo* rhoBC,
                                             doubflo* DD, 
                                             int* k_Q, 
                                             int* k_N, 
                                             int kQ, 
                                             doubflo om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool evenOrOdd);

__global__ void getSendFsPre27(doubflo* DD,
										  doubflo* bufferFs,
										  int* sendIndex,
                                          int buffmax,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat, 
                                          bool evenOrOdd);

__global__ void getSendFsPost27(doubflo* DD,
										   doubflo* bufferFs,
										   int* sendIndex,
                                           int buffmax,
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned int size_Mat, 
                                           bool evenOrOdd);

__global__ void setRecvFsPre27(doubflo* DD,
										  doubflo* bufferFs,
										  int* recvIndex,
                                          int buffmax,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat, 
                                          bool evenOrOdd);

__global__ void setRecvFsPost27(doubflo* DD,
										   doubflo* bufferFs,
										   int* recvIndex,
                                           int buffmax,
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned int size_Mat, 
                                           bool evenOrOdd);


#endif
							 