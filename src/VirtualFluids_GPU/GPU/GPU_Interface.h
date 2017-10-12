#ifndef GPU_Interface_H
#define GPU_Interface_H


int devCheck(int gpudevice);

//////////////////////////////////////////////////////////////////////////
//Kernel
//////////////////////////////////////////////////////////////////////////
void KernelBGKSPSimple27(unsigned int numberOfThreads, 
                                    doubflo s9,
                                    unsigned int* bcMatD,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    doubflo* DD,
                                    int size_Mat,
                                    bool EvenOrOdd);

void InitSP27(unsigned int numberOfThreads,
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
                         bool EvenOrOdd);

void CalcMacCompSP27(doubflo* vxD,
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
								bool evenOrOdd);

void QDevComp27(unsigned int numberOfThreads,
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
						   bool evenOrOdd);

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
							  bool evenOrOdd);

void QPressDevOld27(unsigned int numberOfThreads,
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
                               bool evenOrOdd);

void GetSendFsPreDev27(doubflo* DD,
								  doubflo* bufferFs,
								  int* sendIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat, 
								  bool evenOrOdd,
								  unsigned int numberOfThreads);

void GetSendFsPostDev27(doubflo* DD,
								   doubflo* bufferFs,
								   int* sendIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat, 
								   bool evenOrOdd,
								   unsigned int numberOfThreads);

void SetRecvFsPreDev27(doubflo* DD,
								  doubflo* bufferFs,
								  int* recvIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat, 
								  bool evenOrOdd,
								  unsigned int numberOfThreads);

void SetRecvFsPostDev27(doubflo* DD,
								   doubflo* bufferFs,
								   int* recvIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat, 
								   bool evenOrOdd,
								   unsigned int numberOfThreads);


#endif
