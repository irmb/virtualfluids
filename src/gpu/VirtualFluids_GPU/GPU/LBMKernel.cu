//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ /
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
// includes, cuda
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include "cuda/CudaGrid.h"

// includes, kernels
#include "GPU/GPU_Kernels.cuh"

#include "Parameter/Parameter.h"
//////////////////////////////////////////////////////////////////////////
void KernelCas27( unsigned int grid_nx,
                             unsigned int grid_ny,
                             unsigned int grid_nz,
                             real s9,
                             unsigned int* bcMatD,
                             unsigned int* neighborX,
                             unsigned int* neighborY,
                             unsigned int* neighborZ,
                             real* DD,
                             int size_Mat,
                             bool EvenOrOdd)
{
   dim3 threads       ( grid_nx, 1, 1 );
   dim3 grid          ( grid_ny, grid_nz );   // Gitter fuer Kollision und Propagation

      LB_Kernel_Casc27<<< grid, threads >>>( s9,
                                             bcMatD,
                                             neighborX,
                                             neighborY,
                                             neighborZ,
                                             DD,
                                             size_Mat,
                                             EvenOrOdd);
     getLastCudaError("LB_Kernel_Casc27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelCasSP27( unsigned int numberOfThreads,
                               real s9,
                               unsigned int* bcMatD,
                               unsigned int* neighborX,
                               unsigned int* neighborY,
                               unsigned int* neighborZ,
                               real* DD,
                               int size_Mat,
                               bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

      LB_Kernel_Casc_SP_27<<< grid.grid, grid.threads >>>(s9,
                                                bcMatD,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                DD,
                                                size_Mat,
                                                EvenOrOdd);
      getLastCudaError("LB_Kernel_Casc_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelCasSPMS27( unsigned int numberOfThreads,
                                 real s9,
                                 unsigned int* bcMatD,
                                 unsigned int* neighborX,
                                 unsigned int* neighborY,
                                 unsigned int* neighborZ,
                                 real* DD,
                                 int size_Mat,
                                 bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

      LB_Kernel_Casc_SP_MS_27<<< grid.grid, grid.threads >>>(s9,
                                                   bcMatD,
                                                   neighborX,
                                                   neighborY,
                                                   neighborZ,
                                                   DD,
                                                   size_Mat,
                                                   EvenOrOdd);
      getLastCudaError("LB_Kernel_Casc_SP_MS_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelCasSPMSOHM27( unsigned int numberOfThreads,
                                    real s9,
                                    unsigned int* bcMatD,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    real* DD,
                                    int size_Mat,
                                    bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

      LB_Kernel_Casc_SP_MS_OHM_27<<< grid.grid, grid.threads >>>(  s9,
                                                         bcMatD,
                                                         neighborX,
                                                         neighborY,
                                                         neighborZ,
                                                         DD,
                                                         size_Mat,
                                                         EvenOrOdd);
      getLastCudaError("LB_Kernel_Casc_SP_MS_OHM_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKumCompSRTSP27(
	unsigned int numberOfThreads,
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

   LB_Kernel_Kum_New_Comp_SRT_SP_27 <<< grid.grid, grid.threads >>>(
	   omega,
	   bcMatD,
	   neighborX,
	   neighborY,
	   neighborZ,
	   DDStart,
	   size_Mat,
	   level,
	   forces,
	   EvenOrOdd);
      getLastCudaError("LB_Kernel_Kum_New_Comp_SRT_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKum1hSP27(    unsigned int numberOfThreads,
									real omega,
									real deltaPhi,
									real angularVelocity,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* coordX,
									real* coordY,
									real* coordZ,
									real* DDStart,
									int size_Mat,
									bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

		LB_Kernel_Kum_1h_SP_27<<< grid.grid, grid.threads >>>(omega,
													deltaPhi,
													angularVelocity,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													coordX,
													coordY,
													coordZ,
													DDStart,
													size_Mat,
													EvenOrOdd);
		getLastCudaError("LB_Kernel_Kum_New_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelCascadeSP27(  unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
									bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

		LB_Kernel_Cascade_SP_27<<< grid.grid, grid.threads >>>(s9,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													DD,
													size_Mat,
													EvenOrOdd);
		getLastCudaError("LB_Kernel_Cascade_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKumNewSP27(   unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
									bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


		LB_Kernel_Kum_New_SP_27<<< grid.grid, grid.threads >>>(s9,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													DD,
													size_Mat,
													EvenOrOdd);
		getLastCudaError("LB_Kernel_Kum_New_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKumNewCompSP27(unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
									int size_Array,
									int level,
									real* forces,
									bool EvenOrOdd)
{
	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

		//LB_Kernel_Kum_New_Comp_SP_27<<< grid.grid, grid.threads >>>(	s9,
		//													bcMatD,
		//													neighborX,
		//													neighborY,
		//													neighborZ,
		//													DD,
		//													size_Mat,
		//													level,
		//													forces,
		//													EvenOrOdd);
		//getLastCudaError("LB_Kernel_Kum_New_Comp_SP_27 execution failed");
}

//////////////////////////////////////////////////////////////////////////
void CumulantOnePreconditionedErrorDiffusionChimCompSP27(unsigned int numberOfThreads,
																	real s9,
																	unsigned int* bcMatD,
																	unsigned int* neighborX,
																	unsigned int* neighborY,
																	unsigned int* neighborZ,
																	real* DD,
																	int size_Mat,
																	int size_Array,
																	int level,
																	real* forces,
																	bool EvenOrOdd)
{
	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	Cumulant_One_preconditioned_errorDiffusion_chim_Comp_SP_27 <<< grid.grid, grid.threads >>>(	s9,
																						bcMatD,
																						neighborX,
																						neighborY,
																						neighborZ,
																						DD,
																						size_Mat,
																						level,
																						forces,
																						EvenOrOdd);
		getLastCudaError("Cumulant_One_preconditioned_chim_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CumulantOnePreconditionedChimCompSP27(  unsigned int numberOfThreads,
														real s9,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DD,
														int size_Mat,
														int size_Array,
														int level,
														real* forces,
														bool EvenOrOdd)
{
	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	Cumulant_One_preconditioned_chim_Comp_SP_27 <<< grid.grid, grid.threads >>>(	s9,
																		bcMatD,
																		neighborX,
																		neighborY,
																		neighborZ,
																		DD,
																		size_Mat,
																		level,
																		forces,
																		EvenOrOdd);
		getLastCudaError("Cumulant_One_preconditioned_chim_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CumulantOneChimCompSP27(unsigned int numberOfThreads,
										real s9,
										unsigned int* bcMatD,
										unsigned int* neighborX,
										unsigned int* neighborY,
										unsigned int* neighborZ,
										real* DD,
										int size_Mat,
										int size_Array,
										int level,
										real* forces,
										bool EvenOrOdd)
{
	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	Cumulant_One_chim_Comp_SP_27 <<< grid.grid, grid.threads >>>(	s9,
														bcMatD,
														neighborX,
														neighborY,
														neighborZ,
														DD,
														size_Mat,
														level,
														forces,
														EvenOrOdd);
		getLastCudaError("Cumulant_One_chim_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKumIsoTestSP27(unsigned int numberOfThreads,
									 real s9,
									 unsigned int* bcMatD,
									 unsigned int* neighborX,
									 unsigned int* neighborY,
									 unsigned int* neighborZ,
									 real* DD,
									 real* dxxUx,
									 real* dyyUy,
									 real* dzzUz,
									 int size_Mat,
									 bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	LB_Kernel_Kum_IsoTest_SP_27<<< grid.grid, grid.threads >>>(s9,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													DD,
													dxxUx,
													dyyUy,
													dzzUz,
													size_Mat,
													EvenOrOdd);
	getLastCudaError("LB_Kernel_Kum_IsoTest_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelKumCompSP27(  unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
									bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


		LB_Kernel_Kum_Comp_SP_27<<< grid.grid, grid.threads >>>(s9,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													DD,
													size_Mat,
													EvenOrOdd);
		getLastCudaError("LB_Kernel_Kum_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelPMCumOneCompSP27(unsigned int numberOfThreads,
									   real omega,
									   unsigned int* neighborX,
									   unsigned int* neighborY,
									   unsigned int* neighborZ,
									   real* DD,
									   int size_Mat,
									   int level,
									   real* forces,
									   real porosity,
									   real darcy,
									   real forchheimer,
									   unsigned int sizeOfPorousMedia,
									   unsigned int* nodeIdsPorousMedia,
									   bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	LB_Kernel_PM_Cum_One_Comp_SP_27 <<< grid.grid, grid.threads >>>(omega,
														  neighborX,
														  neighborY,
														  neighborZ,
														  DD,
														  size_Mat,
														  level,
														  forces,
														  porosity,
														  darcy,
														  forchheimer,
														  sizeOfPorousMedia,
														  nodeIdsPorousMedia,
														  EvenOrOdd);
	getLastCudaError("LB_Kernel_PM_Cum_One_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelWaleBySoniMalavCumAA2016CompSP27(
	unsigned int numberOfThreads,
	real s9,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* neighborWSB,
	real* veloX,
	real* veloY,
	real* veloZ,
	real* DD,
	real* turbulentViscosity,
	int size_Mat,
	int size_Array,
	int level,
	real* forces,
	bool EvenOrOdd)
{
	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	LB_Kernel_WaleBySoniMalav_Cum_AA2016_Comp_SP_27 << < grid.grid, grid.threads >> >(
		s9,
		bcMatD,
		neighborX,
		neighborY,
		neighborZ,
		neighborWSB,
		veloX,
		veloY,
		veloZ,
		DD,
		turbulentViscosity,
		size_Mat,
		level,
		forces,
		EvenOrOdd);
	getLastCudaError("LB_Kernel_WaleBySoniMalav_Cum_AA2016_Comp_SP_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelADincomp7(   unsigned int numberOfThreads,
								   real diffusivity,
								   unsigned int* bcMatD,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   real* DD,
								   real* DD7,
								   int size_Mat,
								   bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

      LB_Kernel_AD_Incomp_7<<< grid.grid, grid.threads >>>( diffusivity,
												  bcMatD,
												  neighborX,
												  neighborY,
												  neighborZ,
												  DD,
												  DD7,
												  size_Mat,
												  EvenOrOdd);
      getLastCudaError("LB_Kernel_AD_Incomp_7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void KernelADincomp27( unsigned int numberOfThreads,
								  real diffusivity,
								  unsigned int* bcMatD,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  real* DD,
								  real* DD27,
								  int size_Mat,
								  bool EvenOrOdd)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LB_Kernel_AD_Incomp_27<<< grid.grid, grid.threads >>>( diffusivity,
													bcMatD,
													neighborX,
													neighborY,
													neighborZ,
													DD,
													DD27,
													size_Mat,
													EvenOrOdd);
	getLastCudaError("LB_Kernel_AD_Incomp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Init27( int myid,
                        int numprocs,
                        real u0,
                        unsigned int* geoD,
                        unsigned int* neighborX,
                        unsigned int* neighborY,
                        unsigned int* neighborZ,
                        real* vParab,
                        unsigned int size_Mat,
                        unsigned int grid_nx,
                        unsigned int grid_ny,
                        unsigned int grid_nz,
                        real* DD,
                        int level,
                        int maxlevel)
{
   dim3 threads       ( grid_nx, 1, 1 );
   dim3 grid          ( grid_ny, grid_nz );   // Gitter fuer Kollision und Propagation

	LBInit27<<< grid, threads >>> (  myid,
                                       numprocs,
                                       u0,
                                       geoD,
                                       neighborX,
                                       neighborY,
                                       neighborZ,
                                       vParab,
                                       size_Mat,
                                       grid_nx,
                                       grid_ny,
                                       grid_nz,
                                       DD,
                                       level,
                                       maxlevel);
	getLastCudaError("LBInit27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void InitNonEqPartSP27( unsigned int numberOfThreads,
                                   unsigned int* neighborX,
                                   unsigned int* neighborY,
                                   unsigned int* neighborZ,
                                   unsigned int* neighborWSB,
                                   unsigned int* geoD,
                                   real* rho,
                                   real* ux,
                                   real* uy,
                                   real* uz,
                                   unsigned int size_Mat,
                                   real* DD,
                                   real omega,
                                   bool EvenOrOdd)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBInitNonEqPartSP27<<< grid.grid, grid.threads >>>( neighborX,
                                                neighborY,
                                                neighborZ,
                                                neighborWSB,
                                                geoD,
                                                rho,
                                                ux,
                                                uy,
                                                uz,
                                                size_Mat,
                                                DD,
                                                omega,
                                                EvenOrOdd);
	getLastCudaError("LBInitNonEqPartSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void InitThS7(     unsigned int numberOfThreads,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int* geoD,
                              real* Conc,
                              real* ux,
                              real* uy,
                              real* uz,
                              unsigned int size_Mat,
                              real* DD7,
                              bool EvenOrOdd)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	InitAD7<<< grid.grid, grid.threads >>>( neighborX,
                                       neighborY,
                                       neighborZ,
                                       geoD,
                                       Conc,
                                       ux,
                                       uy,
                                       uz,
                                       size_Mat,
                                       DD7,
                                       EvenOrOdd);
	getLastCudaError("InitAD7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void InitADDev27( unsigned int numberOfThreads,
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           unsigned int* geoD,
                           real* Conc,
                           real* ux,
                           real* uy,
                           real* uz,
                           unsigned int size_Mat,
                           real* DD27,
                           bool EvenOrOdd)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	InitAD27<<< grid.grid, grid.threads >>>(neighborX,
                                       neighborY,
                                       neighborZ,
                                       geoD,
                                       Conc,
                                       ux,
                                       uy,
                                       uz,
                                       size_Mat,
                                       DD27,
                                       EvenOrOdd);
	getLastCudaError("InitAD27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PostProcessorF3_2018Fehlberg(
	unsigned int numberOfThreads,
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* rhoOut,
	real* vxOut,
	real* vyOut,
	real* vzOut,
	real* DDStart,
	real* G6,
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LB_PostProcessor_F3_2018_Fehlberg <<< grid.grid, grid.threads >>> (   omega,
																  bcMatD,
																  neighborX,
																  neighborY,
																  neighborZ,
																  rhoOut,
																  vxOut,
																  vyOut,
																  vzOut,
																  DDStart,
																  G6,
																  size_Mat,
																  level,
																  forces,
																  EvenOrOdd);
	getLastCudaError("LB_PostProcessor_F3_2018_Fehlberg execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMac27( real* vxD,
                           real* vyD,
                           real* vzD,
                           real* rhoD,
                           unsigned int* geoD,
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           unsigned int size_Mat,
                           unsigned int grid_nx,
                           unsigned int grid_ny,
                           unsigned int grid_nz,
                           real* DD,
                           bool isEvenTimestep)
{
   dim3 threads       ( grid_nx, 1, 1 );
   dim3 grid          ( grid_ny, grid_nz );

	LBCalcMac27<<< grid, threads >>> (  vxD,
                                          vyD,
                                          vzD,
                                          rhoD,
                                          geoD,
                                          neighborX,
                                          neighborY,
                                          neighborZ,
                                          size_Mat,
                                          DD,
                                          isEvenTimestep);
	getLastCudaError("LBCalcMac27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacSP27( real* vxD,
                             real* vyD,
                             real* vzD,
                             real* rhoD,
                             real* pressD,
                             unsigned int* geoD,
                             unsigned int* neighborX,
                             unsigned int* neighborY,
                             unsigned int* neighborZ,
                             unsigned int size_Mat,
                             unsigned int numberOfThreads,
                             real* DD,
                             bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMacSP27<<< grid.grid, grid.threads >>> (   vxD,
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
                                             isEvenTimestep);
	getLastCudaError("LBCalcMacSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacCompSP27( real* vxD,
								 real* vyD,
								 real* vzD,
								 real* rhoD,
								 real* pressD,
								 unsigned int* geoD,
								 unsigned int* neighborX,
								 unsigned int* neighborY,
								 unsigned int* neighborZ,
								 unsigned int size_Mat,
								 unsigned int numberOfThreads,
								 real* DD,
								 bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMacCompSP27<<< grid.grid, grid.threads >>> (   vxD,
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
												 isEvenTimestep);
	getLastCudaError("LBCalcMacSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacThS7(  real* Conc,
                              unsigned int* geoD,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              unsigned int numberOfThreads,
                              real* DD7,
                              bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	CalcConc7<<< grid.grid, grid.threads >>> (Conc,
                                          geoD,
                                          neighborX,
                                          neighborY,
                                          neighborZ,
                                          size_Mat,
                                          DD7,
                                          isEvenTimestep);
	getLastCudaError("CalcConc7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PlaneConcThS7(real* Conc,
							  int* kPC,
							  unsigned int numberOfPointskPC,
							  unsigned int* geoD,
							  unsigned int* neighborX,
							  unsigned int* neighborY,
							  unsigned int* neighborZ,
							  unsigned int size_Mat,
                              unsigned int numberOfThreads,
							  real* DD7,
							  bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

	GetPlaneConc7<<< grid.grid, grid.threads >>> (	Conc,
												kPC,
												numberOfPointskPC,
												geoD,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												DD7,
												isEvenTimestep);
	getLastCudaError("GetPlaneConc7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PlaneConcThS27(real* Conc,
							   int* kPC,
							   unsigned int numberOfPointskPC,
							   unsigned int* geoD,
							   unsigned int* neighborX,
							   unsigned int* neighborY,
							   unsigned int* neighborZ,
							   unsigned int size_Mat,
                               unsigned int numberOfThreads,
							   real* DD27,
							   bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

	GetPlaneConc27<<< grid.grid, grid.threads >>> (	Conc,
												kPC,
												numberOfPointskPC,
												geoD,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												DD27,
												isEvenTimestep);
	getLastCudaError("GetPlaneConc27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcConcentration27( unsigned int numberOfThreads,
                                     real* Conc,
                                     unsigned int* geoD,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int size_Mat,
                                     real* DD27,
                                     bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	CalcConc27<<< grid.grid, grid.threads >>> (  Conc,
                                             geoD,
                                             neighborX,
                                             neighborY,
                                             neighborZ,
                                             size_Mat,
                                             DD27,
                                             isEvenTimestep);
	getLastCudaError("CalcConc27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMedSP27(  real* vxD,
                              real* vyD,
                              real* vzD,
                              real* rhoD,
                              real* pressD,
                              unsigned int* geoD,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              unsigned int numberOfThreads,
                              real* DD,
                              bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMedSP27<<< grid.grid, grid.threads >>> (   vxD,
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
                                             isEvenTimestep);
	getLastCudaError("LBCalcMedSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMedCompSP27(  real* vxD,
								  real* vyD,
								  real* vzD,
								  real* rhoD,
								  real* pressD,
								  unsigned int* geoD,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat,
								  unsigned int numberOfThreads,
								  real* DD,
								  bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMedCompSP27<<< grid.grid, grid.threads >>> (   vxD,
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
												 isEvenTimestep);
	getLastCudaError("LBCalcMedSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMedCompAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int* geoD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	unsigned int numberOfThreads,
	real* DD,
	real* DD_AD,
	bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMedCompAD27 <<< grid.grid, grid.threads >>> (
		vxD,
		vyD,
		vzD,
		rhoD,
		pressD,
		concD,
		geoD,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		DD,
		DD_AD,
		isEvenTimestep);
	getLastCudaError("LBCalcMedAD27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacMedSP27(  real* vxD,
                                 real* vyD,
                                 real* vzD,
                                 real* rhoD,
                                 real* pressD,
                                 unsigned int* geoD,
                                 unsigned int* neighborX,
                                 unsigned int* neighborY,
                                 unsigned int* neighborZ,
                                 unsigned int tdiff,
                                 unsigned int size_Mat,
                                 unsigned int numberOfThreads,
                                 bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcMacMedSP27<<< grid.grid, grid.threads >>> (   vxD,
                                                vyD,
                                                vzD,
                                                rhoD,
                                                pressD,
                                                geoD,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                tdiff,
                                                size_Mat,
                                                isEvenTimestep);
	getLastCudaError("LBCalcMacMedSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMedianValuesSP27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	unsigned int size_Mat,
	unsigned int numberOfThreads,
	bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);


	LBResetMedianValuesSP27 << < grid.grid, grid.threads >> > (
		vxD,
		vyD,
		vzD,
		rhoD,
		pressD,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("LBResetMedianValuesSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMedianValuesAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int size_Mat,
	unsigned int numberOfThreads,
	bool isEvenTimestep)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBResetMedianValuesAD27 << < grid.grid, grid.threads >> > (
		vxD,
		vyD,
		vzD,
		rhoD,
		pressD,
		concD,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("LBResetMedianValuesAD27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc2ndMomentsIncompSP27(real* kxyFromfcNEQ,
										 real* kyzFromfcNEQ,
										 real* kxzFromfcNEQ,
										 real* kxxMyyFromfcNEQ,
										 real* kxxMzzFromfcNEQ,
										 unsigned int* geoD,
										 unsigned int* neighborX,
										 unsigned int* neighborY,
										 unsigned int* neighborZ,
										 unsigned int size_Mat,
										 unsigned int numberOfThreads,
										 real* DD,
										 bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalc2ndMomentsIncompSP27<<< grid.grid, grid.threads >>> (  kxyFromfcNEQ,
														 kyzFromfcNEQ,
														 kxzFromfcNEQ,
														 kxxMyyFromfcNEQ,
														 kxxMzzFromfcNEQ,
														 geoD,
														 neighborX,
														 neighborY,
														 neighborZ,
														 size_Mat,
														 DD,
														 isEvenTimestep);
	getLastCudaError("LBCalc2ndMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc2ndMomentsCompSP27( real* kxyFromfcNEQ,
										real* kyzFromfcNEQ,
										real* kxzFromfcNEQ,
										real* kxxMyyFromfcNEQ,
										real* kxxMzzFromfcNEQ,
										unsigned int* geoD,
										unsigned int* neighborX,
										unsigned int* neighborY,
										unsigned int* neighborZ,
										unsigned int size_Mat,
										unsigned int numberOfThreads,
										real* DD,
										bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalc2ndMomentsCompSP27<<< grid.grid, grid.threads >>> (kxyFromfcNEQ,
													 kyzFromfcNEQ,
													 kxzFromfcNEQ,
													 kxxMyyFromfcNEQ,
													 kxxMzzFromfcNEQ,
													 geoD,
													 neighborX,
													 neighborY,
													 neighborZ,
													 size_Mat,
													 DD,
													 isEvenTimestep);
	getLastCudaError("LBCalc2ndMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc3rdMomentsIncompSP27(real* CUMbbb,
										 real* CUMabc,
										 real* CUMbac,
										 real* CUMbca,
										 real* CUMcba,
										 real* CUMacb,
										 real* CUMcab,
										 unsigned int* geoD,
										 unsigned int* neighborX,
										 unsigned int* neighborY,
										 unsigned int* neighborZ,
										 unsigned int size_Mat,
										 unsigned int numberOfThreads,
										 real* DD,
										 bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalc3rdMomentsIncompSP27<<< grid.grid, grid.threads >>> (  CUMbbb,
														 CUMabc,
														 CUMbac,
														 CUMbca,
														 CUMcba,
														 CUMacb,
														 CUMcab,
														 geoD,
														 neighborX,
														 neighborY,
														 neighborZ,
														 DD,
														 size_Mat,
														 isEvenTimestep);
	getLastCudaError("LBCalc3rdMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc3rdMomentsCompSP27( real* CUMbbb,
										real* CUMabc,
										real* CUMbac,
										real* CUMbca,
										real* CUMcba,
										real* CUMacb,
										real* CUMcab,
										unsigned int* geoD,
										unsigned int* neighborX,
										unsigned int* neighborY,
										unsigned int* neighborZ,
										unsigned int size_Mat,
										unsigned int numberOfThreads,
										real* DD,
										bool isEvenTimestep)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalc3rdMomentsCompSP27<<< grid.grid, grid.threads >>> (CUMbbb,
													 CUMabc,
													 CUMbac,
													 CUMbca,
													 CUMcba,
													 CUMacb,
													 CUMcab,
													 geoD,
													 neighborX,
													 neighborY,
													 neighborZ,
													 DD,
													 size_Mat,
													 isEvenTimestep);
	getLastCudaError("LBCalc3rdMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcHigherMomentsIncompSP27(real* CUMcbb,
											real* CUMbcb,
											real* CUMbbc,
											real* CUMcca,
											real* CUMcac,
											real* CUMacc,
											real* CUMbcc,
											real* CUMcbc,
											real* CUMccb,
											real* CUMccc,
											unsigned int* geoD,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat,
											unsigned int numberOfThreads,
											real* DD,
											bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcHigherMomentsIncompSP27<<< grid.grid, grid.threads >>> (CUMcbb,
														  CUMbcb,
														  CUMbbc,
														  CUMcca,
														  CUMcac,
														  CUMacc,
														  CUMbcc,
														  CUMcbc,
														  CUMccb,
														  CUMccc,
														  geoD,
														  neighborX,
														  neighborY,
														  neighborZ,
														  DD,
														  size_Mat,
														  isEvenTimestep);
	getLastCudaError("LBCalcHigherMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcHigherMomentsCompSP27(  real* CUMcbb,
											real* CUMbcb,
											real* CUMbbc,
											real* CUMcca,
											real* CUMcac,
											real* CUMacc,
											real* CUMbcc,
											real* CUMcbc,
											real* CUMccb,
											real* CUMccc,
											unsigned int* geoD,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat,
											unsigned int numberOfThreads,
											real* DD,
											bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);

	LBCalcHigherMomentsCompSP27<<< grid.grid, grid.threads >>> (  CUMcbb,
														  CUMbcb,
														  CUMbbc,
														  CUMcca,
														  CUMcac,
														  CUMacc,
														  CUMbcc,
														  CUMcbc,
														  CUMccb,
														  CUMccc,
														  geoD,
														  neighborX,
														  neighborY,
														  neighborZ,
														  DD,
														  size_Mat,
														  isEvenTimestep);
	getLastCudaError("LBCalcHigherMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void LBCalcMeasurePoints27(real* vxMP,
                                      real* vyMP,
                                      real* vzMP,
                                      real* rhoMP,
                                      unsigned int* kMP,
                                      unsigned int numberOfPointskMP,
                                      unsigned int MPClockCycle,
                                      unsigned int t,
                                      unsigned int* geoD,
                                      unsigned int* neighborX,
                                      unsigned int* neighborY,
                                      unsigned int* neighborZ,
                                      unsigned int size_Mat,
                                      real* DD,
                                      unsigned int numberOfThreads,
                                      bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskMP);

	LBCalcMeasurePoints<<< grid.grid, grid.threads >>> (vxMP,
                                                vyMP,
                                                vzMP,
                                                rhoMP,
                                                kMP,
                                                numberOfPointskMP,
                                                MPClockCycle,
                                                t,
                                                geoD,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                size_Mat,
                                                DD,
                                                isEvenTimestep);
	getLastCudaError("LBCalcMeasurePoints execution failed");
}
//////////////////////////////////////////////////////////////////////////
void BcPress27( int nx,
                           int ny,
                           int tz,
                           unsigned int grid_nx,
                           unsigned int grid_ny,
                           unsigned int* bcMatD,
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           real* DD,
                           unsigned int size_Mat,
                           bool isEvenTimestep)
{
	dim3 threads       ( grid_nx, 1, 1 );
	dim3 grid          ( grid_ny, 1 );

	LB_BC_Press_East27<<< grid, threads >>> ( nx,
                                                ny,
                                                tz,
                                                bcMatD,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                DD,
                                                size_Mat,
                                                isEvenTimestep);
	getLastCudaError("LB_BC_Press_East27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void BcVel27(int nx,
                        int ny,
                        int nz,
                        int itz,
                        unsigned int grid_nx,
                        unsigned int grid_ny,
                        unsigned int* bcMatD,
                        unsigned int* neighborX,
                        unsigned int* neighborY,
                        unsigned int* neighborZ,
                        real* DD,
                        unsigned int size_Mat,
                        bool isEvenTimestep,
                        real u0x,
                        real om)
{
	dim3 threads       ( grid_nx, 1, 1 );
	dim3 grid          ( grid_ny, 1 );

	LB_BC_Vel_West_27<<< grid, threads >>> (  nx,
                                                ny,
                                                nz,
                                                itz,
                                                bcMatD,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                DD,
                                                size_Mat,
                                                isEvenTimestep,
                                                u0x,
                                                grid_nx,
                                                grid_ny,
                                                om);
	getLastCudaError("LB_BC_Vel_West_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressDev7( unsigned int numberOfThreads,
                              real* DD,
                              real* DD7,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q,
                              real* QQ,
                              unsigned int numberOfBCnodes,
                              real om1,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADPress7<<< grid.grid, grid.threads >>>( DD,
                                       DD7,
                                       temp,
                                       velo,
                                       diffusivity,
                                       k_Q,
                                       QQ,
                                       numberOfBCnodes,
                                       om1,
                                       neighborX,
                                       neighborY,
                                       neighborZ,
                                       size_Mat,
                                       isEvenTimestep);
	getLastCudaError("QADPress7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressDev27(unsigned int numberOfThreads,
                              real* DD,
                              real* DD27,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q,
                              real* QQ,
                              unsigned int numberOfBCnodes,
                              real om1,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADPress27<<< grid.grid, grid.threads >>>(   DD,
                                          DD27,
                                          temp,
                                          velo,
                                          diffusivity,
                                          k_Q,
                                          QQ,
                                          numberOfBCnodes,
                                          om1,
                                          neighborX,
                                          neighborY,
                                          neighborZ,
                                          size_Mat,
                                          isEvenTimestep);
	getLastCudaError("QADPress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressNEQNeighborDev27(
											unsigned int numberOfThreads,
											real* DD,
											real* DD27,
											int* k_Q,
											int* k_N,
											int numberOfBCnodes,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat,
											bool isEvenTimestep
										)
{

	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

	QADPressNEQNeighbor27<<< grid.grid, grid.threads >>>(
												DD,
												DD27,
												k_Q,
												k_N,
												numberOfBCnodes,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												isEvenTimestep
											  );
   	getLastCudaError("QADPressNEQNeighbor27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVelDev7(unsigned int numberOfThreads,
                           real* DD,
                           real* DD7,
                           real* temp,
                           real* velo,
                           real diffusivity,
                           int* k_Q,
                           real* QQ,
                           unsigned int numberOfBCnodes,
                           real om1,
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           unsigned int size_Mat,
                           bool isEvenTimestep)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADVel7<<< grid.grid, grid.threads >>> (  
                                       DD,
                                       DD7,
                                       temp,
                                       velo,
                                       diffusivity,
                                       k_Q,
                                       QQ,
                                       numberOfBCnodes,
                                       om1,
                                       neighborX,
                                       neighborY,
                                       neighborZ,
                                       size_Mat,
                                       isEvenTimestep);
	getLastCudaError("QADVel7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVelDev27(  unsigned int numberOfThreads,
                              real* DD,
                              real* DD27,
                              real* temp,
                              real* velo,
                              real diffusivity,
                              int* k_Q,
                              real* QQ,
                              unsigned int numberOfBCnodes,
                              real om1,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADVel27<<< grid.grid, grid.threads >>> ( DD,
                                      DD27,
                                      temp,
                                      velo,
                                      diffusivity,
                                      k_Q,
                                      QQ,
                                      numberOfBCnodes,
                                      om1,
                                      neighborX,
                                      neighborY,
                                      neighborZ,
                                      size_Mat,
                                      isEvenTimestep);
      getLastCudaError("QADVel27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADDev7(unsigned int numberOfThreads,
                        real* DD,
                        real* DD7,
                        real* temp,
                        real diffusivity,
                        int* k_Q,
                        real* QQ,
                        unsigned int numberOfBCnodes,
                        real om1,
                        unsigned int* neighborX,
                        unsigned int* neighborY,
                        unsigned int* neighborZ,
                        unsigned int size_Mat,
                        bool isEvenTimestep)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QAD7<<< grid.grid, grid.threads >>> (     DD,
                                       DD7,
                                       temp,
                                       diffusivity,
                                       k_Q,
                                       QQ,
                                       numberOfBCnodes,
                                       om1,
                                       neighborX,
                                       neighborY,
                                       neighborZ,
                                       size_Mat,
                                       isEvenTimestep);
      getLastCudaError("QAD7 execution failed");
}


//////////////////////////////////////////////////////////////////////////
// Other advection diffusion kernels are in kernel factory :(
void FactorizedCentralMomentsAdvectionDiffusionDeviceKernel(
   uint numberOfThreads,
   real omegaDiffusivity,
   uint* typeOfGridNode,
   uint* neighborX,
   uint* neighborY,
   uint* neighborZ,
   real* distributions,
   real* distributionsAD,
   int size_Mat,
   real* forces,
   bool isEvenTimestep)
{
   int Grid = (size_Mat / numberOfThreads) + 1;
   dim3 grid(Grid, 1, 1);
   dim3 threads(numberOfThreads, 1, 1);

   Factorized_Central_Moments_Advection_Diffusion_Device_Kernel <<< grid, threads >>> (
      omegaDiffusivity,
      typeOfGridNode,
      neighborX,
      neighborY,
      neighborZ,
      distributions,
      distributionsAD,
      size_Mat,
      forces,
      isEvenTimestep);
   getLastCudaError("Factorized_Central_Moments_Advection_Diffusion_Device_Kernel execution failed");
}

//////////////////////////////////////////////////////////////////////////
void ADSlipVelDevComp(
	uint numberOfThreads,
	real * normalX,
	real * normalY,
	real * normalZ,
	real * distributions,
	real * distributionsAD,
	int* QindexArray,
	real * Qarrays,
	uint numberOfBCnodes,
	real omegaDiffusivity,
	uint * neighborX,
	uint * neighborY,
	uint * neighborZ,
	uint size_Mat,
	bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

	AD_SlipVelDeviceComp << < grid.grid, grid.threads >> > (
		normalX,
		normalY,
		normalZ,
		distributions,
		distributionsAD,
		QindexArray,
		Qarrays,
		numberOfBCnodes,
		omegaDiffusivity,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("AD_SlipVelDeviceComp execution failed");
}
//////////////////////////////////////////////////////////////////////////

void QADDirichletDev27( unsigned int numberOfThreads,
								   real* DD,
								   real* DD27,
								   real* temp,
								   real diffusivity,
								   int* k_Q,
								   real* QQ,
								   unsigned int numberOfBCnodes,
								   real om1,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep)
{
   	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADDirichlet27<<< grid.grid, grid.threads >>> (
											   DD,
											   DD27,
											   temp,
											   diffusivity,
											   k_Q,
											   QQ,
											   numberOfBCnodes,
											   om1,
											   neighborX,
											   neighborY,
											   neighborZ,
											   size_Mat,
											   isEvenTimestep);
      getLastCudaError("QADDirichletDev27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADBBDev27(unsigned int numberOfThreads,
                           real* DD,
                           real* DD27,
                           real* temp,
                           real diffusivity,
                           int* k_Q,
                           real* QQ,
                           unsigned int numberOfBCnodes,
                           real om1,
                           unsigned int* neighborX,
                           unsigned int* neighborY,
                           unsigned int* neighborZ,
                           unsigned int size_Mat,
                           bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADBB27<<< grid.grid, grid.threads >>> (  DD,
                                       DD27,
                                       temp,
                                       diffusivity,
                                       k_Q,
                                       QQ,
                                       numberOfBCnodes,
                                       om1,
                                       neighborX,
                                       neighborY,
                                       neighborZ,
                                       size_Mat,
                                       isEvenTimestep);
      getLastCudaError("QADBB27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QNoSlipADincompDev7(unsigned int numberOfThreads,
									real* DD,
									real* DD7,
									real* temp,
									real diffusivity,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QNoSlipADincomp7<<< grid.grid, grid.threads >>> (
											   DD,
											   DD7,
											   temp,
											   diffusivity,
											   k_Q,
											   QQ,
											   numberOfBCnodes,
											   om1,
											   neighborX,
											   neighborY,
											   neighborZ,
											   size_Mat,
											   isEvenTimestep);
      getLastCudaError("QNoSlipADincomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QNoSlipADincompDev27(  unsigned int numberOfThreads,
									   real* DD,
									   real* DD27,
									   real* temp,
									   real diffusivity,
									   int* k_Q,
									   real* QQ,
									   unsigned int numberOfBCnodes,
									   real om1,
									   unsigned int* neighborX,
									   unsigned int* neighborY,
									   unsigned int* neighborZ,
									   unsigned int size_Mat,
									   bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QNoSlipADincomp27<<< grid.grid, grid.threads >>> (
											   DD,
											   DD27,
											   temp,
											   diffusivity,
											   k_Q,
											   QQ,
											   numberOfBCnodes,
											   om1,
											   neighborX,
											   neighborY,
											   neighborZ,
											   size_Mat,
											   isEvenTimestep);
      getLastCudaError("QNoSlipADincomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVeloIncompDev7( unsigned int numberOfThreads,
								   real* DD,
								   real* DD7,
								   real* temp,
								   real* velo,
								   real diffusivity,
								   int* k_Q,
								   real* QQ,
								   unsigned int numberOfBCnodes,
								   real om1,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADVeloIncomp7<<< grid.grid, grid.threads >>> ( DD,
	  										   DD7,
											   temp,
											   velo,
											   diffusivity,
											   k_Q,
											   QQ,
											   numberOfBCnodes,
											   om1,
											   neighborX,
											   neighborY,
											   neighborZ,
											   size_Mat,
											   isEvenTimestep);
      getLastCudaError("QADVeloIncomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADVeloIncompDev27(   unsigned int numberOfThreads,
									  real* DD,
									  real* DD27,
									  real* temp,
									  real* velo,
									  real diffusivity,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADVeloIncomp27<<< grid.grid, grid.threads >>> (
											  DD,
											  DD27,
											  temp,
											  velo,
											  diffusivity,
											  k_Q,
											  QQ,
											  numberOfBCnodes,
											  om1,
											  neighborX,
											  neighborY,
											  neighborZ,
											  size_Mat,
											  isEvenTimestep);
      getLastCudaError("QADVeloIncomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressIncompDev7( unsigned int numberOfThreads,
									  real* DD,
									  real* DD7,
									  real* temp,
									  real* velo,
									  real diffusivity,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADPressIncomp7<<< grid.grid, grid.threads >>>(
											   DD,
											   DD7,
											   temp,
											   velo,
											   diffusivity,
											   k_Q,
											   QQ,
											   numberOfBCnodes,
											   om1,
											   neighborX,
											   neighborY,
											   neighborZ,
											   size_Mat,
											   isEvenTimestep);
      getLastCudaError("QADPressIncomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QADPressIncompDev27(  unsigned int numberOfThreads,
									  real* DD,
									  real* DD27,
									  real* temp,
									  real* velo,
									  real diffusivity,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QADPressIncomp27<<< grid.grid, grid.threads >>>(DD, 
	  										  DD27, 
											  temp,
											  velo,
											  diffusivity,
											  k_Q,
											  QQ,
											  numberOfBCnodes,
											  om1,
											  neighborX,
											  neighborY,
											  neighborZ,
											  size_Mat,
											  isEvenTimestep);
      getLastCudaError("QADPressIncomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

      QDevice27<<< grid, threads >>> (
            parameterDevice->distributions.f[0],
            boundaryCondition->k,
            boundaryCondition->q27[0],
            boundaryCondition->numberOfBCnodes,
            parameterDevice->omega,
            parameterDevice->neighborX,
            parameterDevice->neighborY,
            parameterDevice->neighborZ,
            parameterDevice->numberOfNodes,
            parameterDevice->isEvenTimestep);

      getLastCudaError("QDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

      QDeviceComp27<<< grid, threads >>> (
           parameterDevice->distributions.f[0],
           boundaryCondition->k,
           boundaryCondition->q27[0],
           boundaryCondition->numberOfBCnodes,
           parameterDevice->omega,
           parameterDevice->neighborX,
           parameterDevice->neighborY,
           parameterDevice->neighborZ,
           parameterDevice->numberOfNodes,
           parameterDevice->isEvenTimestep);
      getLastCudaError("QDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QDevCompThinWalls27(unsigned int numberOfThreads,
									real* DD,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* geom,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int* neighborWSB,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QDeviceCompThinWallsPartOne27 <<< grid.grid, grid.threads >>> (DD,
														 k_Q,
														 QQ,
														 numberOfBCnodes,
														 om1,
														 neighborX,
														 neighborY,
														 neighborZ,
														 size_Mat,
														 isEvenTimestep);
   getLastCudaError("QDeviceCompThinWallsPartOne27 execution failed");

   QThinWallsPartTwo27 <<< grid.grid, grid.threads >>> ( DD,
												k_Q,
												QQ,
												numberOfBCnodes,
												geom,
												neighborX,
												neighborY,
												neighborZ,
												neighborWSB,
												size_Mat,
												isEvenTimestep);
   getLastCudaError("QThinWallsPartTwo27 execution failed");

}
//////////////////////////////////////////////////////////////////////////
void QDev3rdMomentsComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1);

   QDevice3rdMomentsComp27<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QDevice3rdMomentsComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QDevIncompHighNu27( unsigned int numberOfThreads,
									real* DD,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QDeviceIncompHighNu27<<< grid.grid, grid.threads >>> (
												   DD,
												   k_Q,
												   QQ,
												   numberOfBCnodes,
												   om1,
												   neighborX,
												   neighborY,
												   neighborZ,
												   size_Mat,
												   isEvenTimestep);
      getLastCudaError("QDeviceIncompHighNu27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QDevCompHighNu27(   unsigned int numberOfThreads,
									real* DD,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QDeviceCompHighNu27<<< grid.grid, grid.threads >>> (
												   DD,
												   k_Q,
												   QQ,
												   numberOfBCnodes,
												   om1,
												   neighborX,
												   neighborY,
												   neighborZ,
												   size_Mat,
												   isEvenTimestep);
      getLastCudaError("QDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevicePlainBB27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QVelDevPlainBB27<<< grid, threads >>> (
         boundaryCondition->Vx,
         boundaryCondition->Vy,
         boundaryCondition->Vz,
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QVelDevicePlainBB27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDeviceCouette27(unsigned int numberOfThreads,
									real* vx,
									real* vy,
									real* vz,
									real* DD,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QVelDevCouette27<<< grid.grid, grid.threads >>> ( vx,
												vy,
												vz,
												DD,
												k_Q,
												QQ,
												numberOfBCnodes,
												om1,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												isEvenTimestep);
      getLastCudaError("QVelDevicePlainBB27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevice1h27(   unsigned int numberOfThreads,
								  int nx,
								  int ny,
								  real* vx,
								  real* vy,
								  real* vz,
								  real* DD,
								  int* k_Q,
								  real* QQ,
								  unsigned int numberOfBCnodes,
								  real om1,
								  real Phi,
								  real angularVelocity,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  real* coordX,
								  real* coordY,
								  real* coordZ,
								  unsigned int size_Mat,
								  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QVelDev1h27<<< grid.grid, grid.threads >>> (nx,
                                          ny,
                                          vx,
                                          vy,
                                          vz,
                                          DD,
                                          k_Q,
                                          QQ,
                                          numberOfBCnodes,
                                          om1,
										  Phi,
										  angularVelocity,
                                          neighborX,
                                          neighborY,
                                          neighborZ,
										  coordX,
										  coordY,
										  coordZ,
                                          size_Mat,
                                          isEvenTimestep);
      getLastCudaError("QVelDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

      QVelDevice27<<< grid, threads >>> (
            parameterDevice->nx,
            parameterDevice->ny,
            boundaryCondition->Vx,
            boundaryCondition->Vy,
            boundaryCondition->Vz,
            parameterDevice->distributions.f[0],
            boundaryCondition->k,
            boundaryCondition->q27[0],
            boundaryCondition->numberOfBCnodes,
            parameterDevice->omega,
            parameterDevice->neighborX,
            parameterDevice->neighborY,
            parameterDevice->neighborZ,
            parameterDevice->numberOfNodes,
            parameterDevice->isEvenTimestep);
      getLastCudaError("QVelDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevCompPlusSlip27(unsigned int numberOfThreads,
									  real* vx,
									  real* vy,
									  real* vz,
									  real* DD,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QVelDeviceCompPlusSlip27<<< grid.grid, grid.threads >>> (
													  vx,
													  vy,
													  vz,
													  DD,
													  k_Q,
													  QQ,
													  numberOfBCnodes,
													  om1,
													  neighborX,
													  neighborY,
													  neighborZ,
													  size_Mat,
													  isEvenTimestep);
      getLastCudaError("QVelDeviceCompPlusSlip27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid(parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QVelDeviceComp27<<< grid, threads >>> (
            boundaryCondition->Vx,
            boundaryCondition->Vy,
            boundaryCondition->Vz,
            parameterDevice->distributions.f[0],
            boundaryCondition->k,        
            boundaryCondition->q27[0],
            boundaryCondition->numberOfBCnodes,
            parameterDevice->omega,
            parameterDevice->neighborX,
            parameterDevice->neighborY,
            parameterDevice->neighborZ,
            parameterDevice->numberOfNodes,
            parameterDevice->isEvenTimestep);
   getLastCudaError("QVelDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevCompThinWalls27(unsigned int numberOfThreads,
							           real* vx,
							           real* vy,
							           real* vz,
							           real* DD,
							           int* k_Q,
							           real* QQ,
							           unsigned int numberOfBCnodes,
							           real om1,
									     unsigned int* geom,
							           unsigned int* neighborX,
							           unsigned int* neighborY,
							           unsigned int* neighborZ,
									     unsigned int* neighborWSB,
							           unsigned int size_Mat,
							           bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QVelDeviceCompThinWallsPartOne27<<< grid.grid, grid.threads >>> (vx,
											                  vy,
											                  vz,
											                  DD,
											                  k_Q,
											                  QQ,
											                  numberOfBCnodes,
											                  om1,
											                  neighborX,
											                  neighborY,
											                  neighborZ,
											                  size_Mat,
											                  isEvenTimestep);
   getLastCudaError("QVelDeviceCompThinWallsPartOne27 execution failed");

	QThinWallsPartTwo27 <<< grid.grid, grid.threads >>> (
       DD,
       k_Q,
       QQ,
       numberOfBCnodes,
       geom,
       neighborX,
       neighborY,
       neighborZ,
       neighborWSB,
       size_Mat,
       isEvenTimestep);
   getLastCudaError("QThinWallsPartTwo27 execution failed");
}

void QVelDevCompZeroPress27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QVelDeviceCompZeroPress27<<< grid, threads >>> (
            boundaryCondition->Vx,
            boundaryCondition->Vy,
            boundaryCondition->Vz,
            parameterDevice->distributions.f[0],
            boundaryCondition->k,
            boundaryCondition->q27[0],
            boundaryCondition->numberOfBCnodes,
            parameterDevice->omega,
            parameterDevice->neighborX,
            parameterDevice->neighborY,
            parameterDevice->neighborZ,
            parameterDevice->numberOfNodes,
            parameterDevice->isEvenTimestep);
   getLastCudaError("QVelDeviceCompZeroPress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevIncompHighNu27(unsigned int numberOfThreads,
									  real* vx,
									  real* vy,
									  real* vz,
									  real* DD,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QVelDeviceIncompHighNu27<<< grid.grid, grid.threads >>> (
													  vx,
													  vy,
													  vz,
													  DD,
													  k_Q,
													  QQ,
													  numberOfBCnodes,
													  om1,
													  neighborX,
													  neighborY,
													  neighborZ,
													  size_Mat,
													  isEvenTimestep);
      getLastCudaError("QVelDeviceIncompHighNu27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVelDevCompHighNu27(  unsigned int numberOfThreads,
									  real* vx,
									  real* vy,
									  real* vz,
									  real* DD,
									  int* k_Q,
									  real* QQ,
									  unsigned int numberOfBCnodes,
									  real om1,
									  unsigned int* neighborX,
									  unsigned int* neighborY,
									  unsigned int* neighborZ,
									  unsigned int size_Mat,
									  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      QVelDeviceCompHighNu27<<< grid.grid, grid.threads >>> (
													  vx,
													  vy,
													  vz,
													  DD,
													  k_Q,
													  QQ,
													  numberOfBCnodes,
													  om1,
													  neighborX,
													  neighborY,
													  neighborZ,
													  size_Mat,
													  isEvenTimestep);
      getLastCudaError("QVelDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVeloDevEQ27(unsigned int numberOfThreads,
							 real* VeloX,
							 real* VeloY,
							 real* VeloZ,
							 real* DD,
							 int* k_Q,
							 int numberOfBCnodes,
							 real om1,
							 unsigned int* neighborX,
							 unsigned int* neighborY,
							 unsigned int* neighborZ,
							 unsigned int size_Mat,
							 bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QVeloDeviceEQ27<<< grid.grid, grid.threads >>> (VeloX,
											 VeloY,
											 VeloZ,
											 DD,
											 k_Q,
											 numberOfBCnodes,
											 om1,
											 neighborX,
											 neighborY,
											 neighborZ,
											 size_Mat,
											 isEvenTimestep);
      getLastCudaError("QVeloDeviceEQ27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QVeloStreetDevEQ27(
	uint  numberOfThreads,
	real* veloXfraction,
	real* veloYfraction,
	int*  naschVelo,
	real* DD,
	int*  naschIndex,
	int   numberOfStreetNodes,
	real  velocityRatio,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint  size_Mat,
	bool  isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfStreetNodes);

	QVeloStreetDeviceEQ27 << < grid.grid, grid.threads >> > (
		veloXfraction,
		veloYfraction,
		naschVelo,
		DD,
		naschIndex,
		numberOfStreetNodes,
		velocityRatio,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("QVeloStreetDeviceEQ27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QSlipDevice27<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QSlipDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipDevCompTurbulentViscosity27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );
   
   QSlipDeviceComp27TurbViscosity<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->turbViscosity,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QSlipDeviceComp27TurbViscosity execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipPressureDevCompTurbulentViscosity27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QSlipPressureDeviceComp27TurbViscosity<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->turbViscosity,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QSlipDeviceComp27TurbViscosity execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );
   
   QSlipDeviceComp27<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QSlipDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void BBSlipDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QSlipDeviceComp27<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("BBSlipDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipGeomDevComp27(unsigned int numberOfThreads,
								   real* DD,
								   int* k_Q,
								   real* QQ,
								   unsigned int numberOfBCnodes,
								   real om1,
								   real* NormalX,
								   real* NormalY,
								   real* NormalZ,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep)
{
	vf::cuda::CudaGrid grid(numberOfThreads, numberOfBCnodes);

   QSlipGeomDeviceComp27<<< grid.grid, grid.threads >>> (DD,
												   k_Q,
												   QQ,
												   numberOfBCnodes,
												   om1,
												   NormalX,
												   NormalY,
												   NormalZ,
												   neighborX,
												   neighborY,
												   neighborZ,
												   size_Mat,
												   isEvenTimestep);
   getLastCudaError("QSlipGeomDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QSlipNormDevComp27(unsigned int numberOfThreads,
								   real* DD,
								   int* k_Q,
								   real* QQ,
								   unsigned int numberOfBCnodes,
								   real om1,
								   real* NormalX,
								   real* NormalY,
								   real* NormalZ,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QSlipNormDeviceComp27<<< grid.grid, grid.threads >>> (DD,
												   k_Q,
												   QQ,
												   numberOfBCnodes,
												   om1,
												   NormalX,
												   NormalY,
												   NormalZ,
												   neighborX,
												   neighborY,
												   neighborZ,
												   size_Mat,
												   isEvenTimestep);
      getLastCudaError("QSlipGeomDeviceComp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QStressDevComp27(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
   dim3 grid = vf::cuda::getCudaGrid(  para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

      QStressDeviceComp27<<< grid, threads >>> (
         para->getParD(level)->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->kN,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         para->getParD(level)->omega,
         para->getParD(level)->turbViscosity,
         para->getParD(level)->velocityX,
         para->getParD(level)->velocityY,
         para->getParD(level)->velocityY,
         boundaryCondition->normalX,
         boundaryCondition->normalY,
         boundaryCondition->normalZ,
         boundaryCondition->Vx,
         boundaryCondition->Vy,
         boundaryCondition->Vz,
         boundaryCondition->Vx1,
         boundaryCondition->Vy1,
         boundaryCondition->Vz1,
         para->getParD(level)->wallModel.samplingOffset,
         para->getParD(level)->wallModel.z0,
         para->getHasWallModelMonitor(),
         para->getParD(level)->wallModel.u_star,
         para->getParD(level)->wallModel.Fx,
         para->getParD(level)->wallModel.Fy,
         para->getParD(level)->wallModel.Fz,
         para->getParD(level)->neighborX,
         para->getParD(level)->neighborY,
         para->getParD(level)->neighborZ,
         para->getParD(level)->numberOfNodes,
         para->getParD(level)->isEvenTimestep);
      getLastCudaError("QSlipDeviceComp27 execution failed");
}

//////////////////////////////////////////////////////////////////////////
void BBStressDev27(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
   dim3 grid = vf::cuda::getCudaGrid( para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

   BBStressDevice27<<< grid, threads >>> (
      para->getParD(level)->distributions.f[0],
      boundaryCondition->k,
      boundaryCondition->kN,
      boundaryCondition->q27[0],
      boundaryCondition->numberOfBCnodes,
      para->getParD(level)->velocityX,
      para->getParD(level)->velocityY,
      para->getParD(level)->velocityY,
      boundaryCondition->normalX,
      boundaryCondition->normalY,
      boundaryCondition->normalZ,
      boundaryCondition->Vx,
      boundaryCondition->Vy,
      boundaryCondition->Vz,
      boundaryCondition->Vx1,
      boundaryCondition->Vy1,
      boundaryCondition->Vz1,
      para->getParD(level)->wallModel.samplingOffset,
      para->getParD(level)->wallModel.z0,
      para->getHasWallModelMonitor(),
      para->getParD(level)->wallModel.u_star,
      para->getParD(level)->wallModel.Fx,
      para->getParD(level)->wallModel.Fy,
      para->getParD(level)->wallModel.Fz,
      para->getParD(level)->neighborX,
      para->getParD(level)->neighborY,
      para->getParD(level)->neighborZ,
      para->getParD(level)->numberOfNodes,
      para->getParD(level)->isEvenTimestep);
      getLastCudaError("BBStressDevice27 execution failed");
}

//////////////////////////////////////////////////////////////////////////
void BBStressPressureDev27(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
{
   dim3 grid = vf::cuda::getCudaGrid( para->getParD(level)->numberofthreads, boundaryCondition->numberOfBCnodes);
   dim3 threads(para->getParD(level)->numberofthreads, 1, 1 );

   BBStressPressureDevice27<<< grid, threads >>> (
      para->getParD(level)->distributions.f[0],
      boundaryCondition->k,
      boundaryCondition->kN,
      boundaryCondition->q27[0],
      boundaryCondition->numberOfBCnodes,
      para->getParD(level)->velocityX,
      para->getParD(level)->velocityY,
      para->getParD(level)->velocityY,
      boundaryCondition->normalX,
      boundaryCondition->normalY,
      boundaryCondition->normalZ,
      boundaryCondition->Vx,
      boundaryCondition->Vy,
      boundaryCondition->Vz,
      boundaryCondition->Vx1,
      boundaryCondition->Vy1,
      boundaryCondition->Vz1,
      para->getParD(level)->wallModel.samplingOffset,
      para->getParD(level)->wallModel.z0,
      para->getHasWallModelMonitor(),
      para->getParD(level)->wallModel.u_star,
      para->getParD(level)->wallModel.Fx,
      para->getParD(level)->wallModel.Fy,
      para->getParD(level)->wallModel.Fz,
      para->getParD(level)->neighborX,
      para->getParD(level)->neighborY,
      para->getParD(level)->neighborZ,
      para->getParD(level)->numberOfNodes,
      para->getParD(level)->isEvenTimestep);
      getLastCudaError("BBStressDevice27 execution failed");
}

//////////////////////////////////////////////////////////////////////////
void QPressDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QPressDevice27<<< grid, threads >>> (
      boundaryCondition->RhoBC,
      parameterDevice->distributions.f[0],
      boundaryCondition->k,
      boundaryCondition->q27[0],
      boundaryCondition->numberOfBCnodes,
      parameterDevice->omega,
      parameterDevice->neighborX,
      parameterDevice->neighborY,
      parameterDevice->neighborZ,
      parameterDevice->numberOfNodes,
      parameterDevice->isEvenTimestep);
   getLastCudaError("QPressDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevAntiBB27(  unsigned int numberOfThreads,
                                    real* rhoBC,
									real* vx,
									real* vy,
									real* vz,
									real* DD,
									int* k_Q,
									real* QQ,
									int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPressDeviceAntiBB27<<< grid.grid, grid.threads >>>( rhoBC,
												vx,
												vy,
												vz,
												DD,
												k_Q,
												QQ,
												numberOfBCnodes,
												om1,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												isEvenTimestep);
   getLastCudaError("QPressDeviceAntiBB27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevFixBackflow27( unsigned int numberOfThreads,
                                        real* rhoBC,
                                        real* DD,
                                        int* k_Q,
                                        unsigned int numberOfBCnodes,
                                        real om1,
                                        unsigned int* neighborX,
                                        unsigned int* neighborY,
                                        unsigned int* neighborZ,
                                        unsigned int size_Mat,
                                        bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPressDeviceFixBackflow27<<< grid.grid, grid.threads >>> (  rhoBC,
                                                         DD,
                                                         k_Q,
                                                         numberOfBCnodes,
                                                         om1,
                                                         neighborX,
                                                         neighborY,
                                                         neighborZ,
                                                         size_Mat,
                                                         isEvenTimestep);
   getLastCudaError("QPressDeviceFixBackflow27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevDirDepBot27(  unsigned int numberOfThreads,
                                       real* rhoBC,
                                       real* DD,
                                       int* k_Q,
                                       unsigned int numberOfBCnodes,
                                       real om1,
                                       unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int size_Mat,
                                       bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPressDeviceDirDepBot27<<< grid.grid, grid.threads >>> ( rhoBC,
                                                      DD,
                                                      k_Q,
                                                      numberOfBCnodes,
                                                      om1,
                                                      neighborX,
                                                      neighborY,
                                                      neighborZ,
                                                      size_Mat,
                                                      isEvenTimestep);
   getLastCudaError("QPressDeviceDirDepBot27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressNoRhoDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QPressNoRhoDevice27<<< grid, threads >>> (
         boundaryCondition->RhoBC,
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->kN,
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep,
         vf::lbm::dir::DIR_P00);
   getLastCudaError("QPressNoRhoDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressZeroRhoOutflowDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QPressZeroRhoOutflowDevice27<<< grid, threads >>> (
         boundaryCondition->RhoBC,
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->kN,
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep,
         vf::lbm::dir::DIR_P00,
         parameterDevice->outflowPressureCorrectionFactor);
   getLastCudaError("QPressZeroRhoOutflowDev27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QInflowScaleByPressDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QInflowScaleByPressDevice27<<< grid, threads >>> (
           boundaryCondition->RhoBC,
           parameterDevice->distributions.f[0],
           boundaryCondition->k,
           boundaryCondition->kN,
           boundaryCondition->numberOfBCnodes,
           parameterDevice->omega,
           parameterDevice->neighborX,
           parameterDevice->neighborY,
           parameterDevice->neighborZ,
           parameterDevice->numberOfNodes,
           parameterDevice->isEvenTimestep);
   getLastCudaError("QInflowScaleByPressDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevOld27(  unsigned int numberOfThreads,
                                     real* rhoBC,
                                     real* DD,
                                     int* k_Q,
                                     int* k_N,
                                     unsigned int numberOfBCnodes,
                                     real om1,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int size_Mat,
                                     bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPressDeviceOld27<<< grid.grid, grid.threads >>> ( rhoBC,
                                                DD,
                                                k_Q,
                                                k_N,
                                                numberOfBCnodes,
                                                om1,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                size_Mat,
                                                isEvenTimestep);
   getLastCudaError("QPressDeviceOld27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevIncompNEQ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QPressDeviceIncompNEQ27<<< grid, threads >>> (
         boundaryCondition->RhoBC,
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->kN,
         boundaryCondition->numberOfBCnodes,
         parameterDevice->omega,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("QPressDeviceIncompNEQ27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevNEQ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   QPressDeviceNEQ27<<< grid, threads >>> (
        boundaryCondition->RhoBC,
        parameterDevice->distributions.f[0],
        boundaryCondition->k,
        boundaryCondition->kN,
        boundaryCondition->numberOfBCnodes,
        parameterDevice->omega,
        parameterDevice->neighborX,
        parameterDevice->neighborY,
        parameterDevice->neighborZ,
        parameterDevice->numberOfNodes,
        parameterDevice->isEvenTimestep);
   getLastCudaError("QPressDevNEQ27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevEQZ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

      QPressDeviceEQZ27<<< grid, threads >>> (
            boundaryCondition->RhoBC,
            parameterDevice->distributions.f[0],
            boundaryCondition->k,
            boundaryCondition->kN,
            parameterDevice->kDistTestRE.f[0],
            boundaryCondition->numberOfBCnodes,
            parameterDevice->omega,
            parameterDevice->neighborX,
            parameterDevice->neighborY,
            parameterDevice->neighborZ,
            parameterDevice->numberOfNodes,
            parameterDevice->isEvenTimestep);
      getLastCudaError("QPressDeviceEQZ27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevZero27(unsigned int numberOfThreads,
                                real* DD,
                                int* k_Q,
                                unsigned int numberOfBCnodes,
                                unsigned int* neighborX,
                                unsigned int* neighborY,
                                unsigned int* neighborZ,
                                unsigned int size_Mat,
                                bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPressDeviceZero27<<< grid.grid, grid.threads >>> (DD,
                                                k_Q,
                                                numberOfBCnodes,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                size_Mat,
                                                isEvenTimestep);
   getLastCudaError("QPressDeviceOld27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDevFake27(     unsigned int numberOfThreads,
                                     real* rhoBC,
                                     real* DD,
                                     int* k_Q,
                                     int* k_N,
                                     unsigned int numberOfBCnodes,
                                     real om1,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int size_Mat,
                                     bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);


      QPressDeviceFake27<<< grid.grid, grid.threads >>> (rhoBC,
                                                DD,
                                                k_Q,
                                                k_N,
                                                numberOfBCnodes,
                                                om1,
                                                neighborX,
                                                neighborY,
                                                neighborZ,
                                                size_Mat,
                                                isEvenTimestep);
      getLastCudaError("QPressDeviceFake27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void BBDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
   dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
   dim3 threads(parameterDevice->numberofthreads, 1, 1 );

   BBDevice27<<< grid, threads >>> (
         parameterDevice->distributions.f[0],
         boundaryCondition->k,
         boundaryCondition->q27[0],
         boundaryCondition->numberOfBCnodes,
         parameterDevice->neighborX,
         parameterDevice->neighborY,
         parameterDevice->neighborZ,
         parameterDevice->numberOfNodes,
         parameterDevice->isEvenTimestep);
   getLastCudaError("BBDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPressDev27_IntBB(  unsigned int numberOfThreads,
									real* rho,
									real* DD,
									int* k_Q,
									real* QQ,
									unsigned int numberOfBCnodes,
									real om1,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int size_Mat,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

	QPressDevice27_IntBB<<< grid.grid, grid.threads >>> (rho,
													DD,
													k_Q,
													QQ,
													numberOfBCnodes,
													om1,
													neighborX,
													neighborY,
													neighborZ,
													size_Mat,
													isEvenTimestep);
	getLastCudaError("QPressDevice27_IntBB execution failed");
}
// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
//////////////////////////////////////////////////////////////////////////
void PressSchlaffer27(unsigned int numberOfThreads,
                                 real* rhoBC,
                                 real* DD,
                                 real* vx0,
                                 real* vy0,
                                 real* vz0,
                                 real* deltaVz0,
                                 int* k_Q,
                                 int* k_N,
                                 int numberOfBCnodes,
                                 real om1,
                                 unsigned int* neighborX,
                                 unsigned int* neighborY,
                                 unsigned int* neighborZ,
                                 unsigned int size_Mat,
                                 bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   PressSchlaff27<<< grid.grid, grid.threads >>>(  rhoBC,
                                             DD,
                                             vx0,
                                             vy0,
                                             vz0,
                                             deltaVz0,
                                             k_Q,
                                             k_N,
                                             numberOfBCnodes,
                                             om1,
                                             neighborX,
                                             neighborY,
                                             neighborZ,
                                             size_Mat,
                                             isEvenTimestep);
   getLastCudaError("PressSchlaff27 execution failed");
}
// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
//////////////////////////////////////////////////////////////////////////
void VelSchlaffer27(  unsigned int numberOfThreads,
                                 int t,
                                 real* DD,
                                 real* vz0,
                                 real* deltaVz0,
                                 int* k_Q,
                                 int* k_N,
                                 int numberOfBCnodes,
                                 real om1,
                                 unsigned int* neighborX,
                                 unsigned int* neighborY,
                                 unsigned int* neighborZ,
                                 unsigned int size_Mat,
                                 bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   VelSchlaff27<<< grid.grid, grid.threads >>>( t,
                                          DD,
                                          vz0,
                                          deltaVz0,
                                          k_Q,
                                          k_N,
                                          numberOfBCnodes,
                                          om1,
                                          neighborX,
                                          neighborY,
                                          neighborZ,
                                          size_Mat,
                                          isEvenTimestep);
      getLastCudaError("VelSchlaff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void QPrecursorDevCompZeroPress(  uint numberOfThreads, real tRatio,
                                             real* DD, real* QQ, int* k_Q, 
                                             uint sizeQ, uint numberOfBCnodes,
                                             real omega, real velocityRatio,
                                             uint* neighborX, uint* neighborY, uint* neighborZ,
                                             uint* neighborsNT, uint* neighborsNB, uint* neighborsST, uint* neighborsSB,
                                             real* weightsNT, real* weightsNB, real* weightsST, real* weightsSB,
                                             real* vxLast, real* vyLast, real* vzLast,
                                             real* vxCurrent, real* vyCurrent, real* vzCurrent,
                                             real velocityX, real velocityY, real velocityZ, 
                                             unsigned long long size_Mat, bool evenOrOdd)
{

   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

   QPrecursorDeviceCompZeroPress<<< grid.grid, grid.threads >>>(k_Q, numberOfBCnodes, sizeQ, omega, DD, QQ,
                                                               neighborX, neighborY, neighborZ,
                                                               neighborsNT, neighborsNB, neighborsST, neighborsSB,
                                                               weightsNT, weightsNB, weightsST, weightsSB,
                                                               vxLast, vyLast, vzLast,
                                                               vxCurrent, vyCurrent, vzCurrent, 
                                                               velocityX, velocityY, velocityZ, 
                                                               tRatio, velocityRatio, size_Mat, evenOrOdd);
   getLastCudaError("QPrecursorDeviceCompZeroPress execution failed"); 


}
//////////////////////////////////////////////////////////////////////////
extern "C" void PropVelo(   unsigned int numberOfThreads,
                            unsigned int* neighborX,
                            unsigned int* neighborY,
                            unsigned int* neighborZ,
                            real* rho,
                            real* ux,
                            real* uy,
                            real* uz,
                            int* k_Q,
							unsigned int size_Prop,
                            unsigned int size_Mat,
                            unsigned int* bcMatD,
                            real* DD,
                            bool EvenOrOdd)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Prop);

      PropellerBC<<< grid.grid, grid.threads >>>(neighborX,
                                       neighborY,
                                       neighborZ,
                                       rho,
                                       ux,
                                       uy,
                                       uz,
									   k_Q,
									   size_Prop,
                                       size_Mat,
									   bcMatD,
                                       DD,
                                       EvenOrOdd);
      getLastCudaError("PropellerBC execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF27( real* DC,
                        real* DF,
                        unsigned int* neighborCX,
                        unsigned int* neighborCY,
                        unsigned int* neighborCZ,
                        unsigned int* neighborFX,
                        unsigned int* neighborFY,
                        unsigned int* neighborFZ,
                        unsigned int size_MatC,
                        unsigned int size_MatF,
                        bool isEvenTimestep,
                        unsigned int* posCSWB,
                        unsigned int* posFSWB,
                        unsigned int kCF,
                        real omCoarse,
                        real omFine,
                        real nu,
                        unsigned int nxC,
                        unsigned int nyC,
                        unsigned int nxF,
                        unsigned int nyF,
                        unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);
   
      scaleCF27<<< grid.grid, grid.threads >>> ( DC,
                                             DF,
                                             neighborCX,
                                             neighborCY,
                                             neighborCZ,
                                             neighborFX,
                                             neighborFY,
                                             neighborFZ,
                                             size_MatC,
                                             size_MatF,
                                             isEvenTimestep,
                                             posCSWB,
                                             posFSWB,
                                             kCF,
                                             omCoarse,
                                             omFine,
                                             nu,
                                             nxC,
                                             nyC,
                                             nxF,
                                             nyF);
      getLastCudaError("scaleCF27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFEff27(real* DC,
                             real* DF,
                             unsigned int* neighborCX,
                             unsigned int* neighborCY,
                             unsigned int* neighborCZ,
                             unsigned int* neighborFX,
                             unsigned int* neighborFY,
                             unsigned int* neighborFZ,
                             unsigned int size_MatC,
                             unsigned int size_MatF,
                             bool isEvenTimestep,
                             unsigned int* posCSWB,
                             unsigned int* posFSWB,
                             unsigned int kCF,
                             real omCoarse,
                             real omFine,
                             real nu,
                             unsigned int nxC,
                             unsigned int nyC,
                             unsigned int nxF,
                             unsigned int nyF,
                             unsigned int numberOfThreads,
                             OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFEff27<<< grid.grid, grid.threads >>> ( DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offCF);
      getLastCudaError("scaleCFEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFLast27(real* DC,
                              real* DF,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posCSWB,
                              unsigned int* posFSWB,
                              unsigned int kCF,
                              real omCoarse,
                              real omFine,
                              real nu,
                              unsigned int nxC,
                              unsigned int nyC,
                              unsigned int nxF,
                              unsigned int nyF,
                              unsigned int numberOfThreads,
                              OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFLast27<<< grid.grid, grid.threads >>> (DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offCF);
      getLastCudaError("scaleCFLast27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFpress27(  real* DC,
                                 real* DF,
                                 unsigned int* neighborCX,
                                 unsigned int* neighborCY,
                                 unsigned int* neighborCZ,
                                 unsigned int* neighborFX,
                                 unsigned int* neighborFY,
                                 unsigned int* neighborFZ,
                                 unsigned int size_MatC,
                                 unsigned int size_MatF,
                                 bool isEvenTimestep,
                                 unsigned int* posCSWB,
                                 unsigned int* posFSWB,
                                 unsigned int kCF,
                                 real omCoarse,
                                 real omFine,
                                 real nu,
                                 unsigned int nxC,
                                 unsigned int nyC,
                                 unsigned int nxF,
                                 unsigned int nyF,
                                 unsigned int numberOfThreads,
                                 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFpress27<<< grid.grid, grid.threads >>>(DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offCF);
      getLastCudaError("scaleCFpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_27(  real* DC,
                                 real* DF,
                                 unsigned int* neighborCX,
                                 unsigned int* neighborCY,
                                 unsigned int* neighborCZ,
                                 unsigned int* neighborFX,
                                 unsigned int* neighborFY,
                                 unsigned int* neighborFZ,
                                 unsigned int size_MatC,
                                 unsigned int size_MatF,
                                 bool isEvenTimestep,
                                 unsigned int* posCSWB,
                                 unsigned int* posFSWB,
                                 unsigned int kCF,
                                 real omCoarse,
                                 real omFine,
                                 real nu,
                                 unsigned int nxC,
                                 unsigned int nyC,
                                 unsigned int nxF,
                                 unsigned int nyF,
                                 unsigned int numberOfThreads,
                                 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_Fix_27<<< grid.grid, grid.threads >>>(DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offCF);
      getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_comp_27( real* DC,
									 real* DF,
									 unsigned int* neighborCX,
									 unsigned int* neighborCY,
									 unsigned int* neighborCZ,
									 unsigned int* neighborFX,
									 unsigned int* neighborFY,
									 unsigned int* neighborFZ,
									 unsigned int size_MatC,
									 unsigned int size_MatF,
									 bool isEvenTimestep,
									 unsigned int* posCSWB,
									 unsigned int* posFSWB,
									 unsigned int kCF,
									 real omCoarse,
									 real omFine,
									 real nu,
									 unsigned int nxC,
									 unsigned int nyC,
									 unsigned int nxF,
									 unsigned int nyF,
									 unsigned int numberOfThreads,
									 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_Fix_comp_27<<< grid.grid, grid.threads >>>(   DC,
														DF,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posCSWB,
														posFSWB,
														kCF,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offCF);
      getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_0817_comp_27(real* DC,
									 real* DF,
									 unsigned int* neighborCX,
									 unsigned int* neighborCY,
									 unsigned int* neighborCZ,
									 unsigned int* neighborFX,
									 unsigned int* neighborFY,
									 unsigned int* neighborFZ,
									 unsigned int size_MatC,
									 unsigned int size_MatF,
									 bool isEvenTimestep,
									 unsigned int* posCSWB,
									 unsigned int* posFSWB,
									 unsigned int kCF,
									 real omCoarse,
									 real omFine,
									 real nu,
									 unsigned int nxC,
									 unsigned int nyC,
									 unsigned int nxF,
									 unsigned int nyF,
									 unsigned int numberOfThreads,
									 OffCF offCF,
                            CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_0817_comp_27<<< grid.grid, grid.threads, 0, stream >>>(  DC,
														DF,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posCSWB,
														posFSWB,
														kCF,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offCF);
      getLastCudaError("scaleCF_0817_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3_2018(real* DC,
										  real* DF,
										  real* G6,
										  unsigned int* neighborCX,
										  unsigned int* neighborCY,
										  unsigned int* neighborCZ,
										  unsigned int* neighborFX,
										  unsigned int* neighborFY,
										  unsigned int* neighborFZ,
										  unsigned int size_MatC,
										  unsigned int size_MatF,
										  bool isEvenTimestep,
										  unsigned int* posCSWB,
										  unsigned int* posFSWB,
										  unsigned int kCF,
										  real omCoarse,
										  real omFine,
										  real nu,
										  unsigned int nxC,
										  unsigned int nyC,
										  unsigned int nxF,
										  unsigned int nyF,
										  unsigned int numberOfThreads,
										  OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_comp_D3Q27F3_2018 <<< grid.grid, grid.threads >>>(DC,
															DF,
															G6,
															neighborCX,
															neighborCY,
															neighborCZ,
															neighborFX,
															neighborFY,
															neighborFZ,
															size_MatC,
															size_MatF,
															isEvenTimestep,
															posCSWB,
															posFSWB,
															kCF,
															omCoarse,
															omFine,
															nu,
															nxC,
															nyC,
															nxF,
															nyF,
															offCF);
      getLastCudaError("scaleCF_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3(real* DC,
									 real* DF,
									 real* G6,
									 unsigned int* neighborCX,
									 unsigned int* neighborCY,
									 unsigned int* neighborCZ,
									 unsigned int* neighborFX,
									 unsigned int* neighborFY,
									 unsigned int* neighborFZ,
									 unsigned int size_MatC,
									 unsigned int size_MatF,
									 bool isEvenTimestep,
									 unsigned int* posCSWB,
									 unsigned int* posFSWB,
									 unsigned int kCF,
									 real omCoarse,
									 real omFine,
									 real nu,
									 unsigned int nxC,
									 unsigned int nyC,
									 unsigned int nxF,
									 unsigned int nyF,
									 unsigned int numberOfThreads,
									 OffCF offCF,
                            CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_comp_D3Q27F3 <<< grid.grid, grid.threads, 0, stream >>>( DC,
														DF,
														G6,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posCSWB,
														posFSWB,
														kCF,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offCF);
      getLastCudaError("scaleCF_comp_D3Q27F3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_staggered_time_comp_27(  real* DC,
												 real* DF,
												 unsigned int* neighborCX,
												 unsigned int* neighborCY,
												 unsigned int* neighborCZ,
												 unsigned int* neighborFX,
												 unsigned int* neighborFY,
												 unsigned int* neighborFZ,
												 unsigned int size_MatC,
												 unsigned int size_MatF,
												 bool isEvenTimestep,
												 unsigned int* posCSWB,
												 unsigned int* posFSWB,
												 unsigned int kCF,
												 real omCoarse,
												 real omFine,
												 real nu,
												 unsigned int nxC,
												 unsigned int nyC,
												 unsigned int nxF,
												 unsigned int nyF,
												 unsigned int numberOfThreads,
												 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_staggered_time_comp_27<<< grid.grid, grid.threads >>>(    DC,
																	DF,
																	neighborCX,
																	neighborCY,
																	neighborCZ,
																	neighborFX,
																	neighborFY,
																	neighborFZ,
																	size_MatC,
																	size_MatF,
																	isEvenTimestep,
																	posCSWB,
																	posFSWB,
																	kCF,
																	omCoarse,
																	omFine,
																	nu,
																	nxC,
																	nyC,
																	nxF,
																	nyF,
																	offCF);
      getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_comp_27(   real* DC,
										 real* DF,
										 unsigned int* neighborCX,
										 unsigned int* neighborCY,
										 unsigned int* neighborCZ,
										 unsigned int* neighborFX,
										 unsigned int* neighborFY,
										 unsigned int* neighborFZ,
										 unsigned int size_MatC,
										 unsigned int size_MatF,
										 bool isEvenTimestep,
										 unsigned int* posCSWB,
										 unsigned int* posFSWB,
										 unsigned int kCF,
										 real omCoarse,
										 real omFine,
										 real nu,
										 unsigned int nxC,
										 unsigned int nyC,
										 unsigned int nxF,
										 unsigned int nyF,
										 unsigned int numberOfThreads,
										 OffCF offCF,
                               CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_RhoSq_comp_27<<< grid.grid, grid.threads, 0, stream >>>( DC,
														DF,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posCSWB,
														posFSWB,
														kCF,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offCF);
      getLastCudaError("scaleCF_RhoSq_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_3rdMom_comp_27(real* DC,
											 real* DF,
											 unsigned int* neighborCX,
											 unsigned int* neighborCY,
											 unsigned int* neighborCZ,
											 unsigned int* neighborFX,
											 unsigned int* neighborFY,
											 unsigned int* neighborFZ,
											 unsigned int size_MatC,
											 unsigned int size_MatF,
											 bool isEvenTimestep,
											 unsigned int* posCSWB,
											 unsigned int* posFSWB,
											 unsigned int kCF,
											 real omCoarse,
											 real omFine,
											 real nu,
											 unsigned int nxC,
											 unsigned int nyC,
											 unsigned int nxF,
											 unsigned int nyF,
											 unsigned int numberOfThreads,
											 OffCF offCF,
                                  CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_RhoSq_3rdMom_comp_27<<< grid.grid, grid.threads, 0, stream >>>(  DC,
																DF,
																neighborCX,
																neighborCY,
																neighborCZ,
																neighborFX,
																neighborFY,
																neighborFZ,
																size_MatC,
																size_MatF,
																isEvenTimestep,
																posCSWB,
																posFSWB,
																kCF,
																omCoarse,
																omFine,
																nu,
																nxC,
																nyC,
																nxF,
																nyF,
																offCF);
      getLastCudaError("scaleCF_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_AA2016_comp_27(real* DC,
									   real* DF,
									   unsigned int* neighborCX,
									   unsigned int* neighborCY,
									   unsigned int* neighborCZ,
									   unsigned int* neighborFX,
									   unsigned int* neighborFY,
									   unsigned int* neighborFZ,
									   unsigned int size_MatC,
									   unsigned int size_MatF,
									   bool isEvenTimestep,
									   unsigned int* posCSWB,
									   unsigned int* posFSWB,
									   unsigned int kCF,
									   real omCoarse,
									   real omFine,
									   real nu,
									   unsigned int nxC,
									   unsigned int nyC,
									   unsigned int nxF,
									   unsigned int nyF,
									   unsigned int numberOfThreads,
									   OffCF offCF,
                              CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_AA2016_comp_27<<< grid.grid, grid.threads, 0, stream >>>(DC,
														DF,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posCSWB,
														posFSWB,
														kCF,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offCF);
      getLastCudaError("scaleCF_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_NSPress_27(  real* DC,
									 real* DF,
									 unsigned int* neighborCX,
									 unsigned int* neighborCY,
									 unsigned int* neighborCZ,
									 unsigned int* neighborFX,
									 unsigned int* neighborFY,
									 unsigned int* neighborFZ,
									 unsigned int size_MatC,
									 unsigned int size_MatF,
									 bool isEvenTimestep,
									 unsigned int* posCSWB,
									 unsigned int* posFSWB,
									 unsigned int kCF,
									 real omCoarse,
									 real omFine,
									 real nu,
									 unsigned int nxC,
									 unsigned int nyC,
									 unsigned int nxF,
									 unsigned int nyF,
									 unsigned int numberOfThreads,
									 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCF_NSPress_27<<< grid.grid, grid.threads >>>(DC,
													DF,
													neighborCX,
													neighborCY,
													neighborCZ,
													neighborFX,
													neighborFY,
													neighborFZ,
													size_MatC,
													size_MatF,
													isEvenTimestep,
													posCSWB,
													posFSWB,
													kCF,
													omCoarse,
													omFine,
													nu,
													nxC,
													nyC,
													nxF,
													nyF,
													offCF);
      getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThSMG7(   real* DC,
                                 real* DF,
                                 real* DD7C,
                                 real* DD7F,
                                 unsigned int* neighborCX,
                                 unsigned int* neighborCY,
                                 unsigned int* neighborCZ,
                                 unsigned int* neighborFX,
                                 unsigned int* neighborFY,
                                 unsigned int* neighborFZ,
                                 unsigned int size_MatC,
                                 unsigned int size_MatF,
                                 bool isEvenTimestep,
                                 unsigned int* posCSWB,
                                 unsigned int* posFSWB,
                                 unsigned int kCF,
                                 real nu,
                                 real diffusivity_fine,
                                 unsigned int numberOfThreads,
                                 OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFThSMG7<<< grid.grid, grid.threads >>> (DC,
                                                DF,
                                                DD7C,
                                                DD7F,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                nu,
                                                diffusivity_fine,
                                                offCF);
      getLastCudaError("scaleCFThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS7(  real* DC,
                              real* DF,
                              real* DD7C,
                              real* DD7F,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posCSWB,
                              unsigned int* posFSWB,
                              unsigned int kCF,
                              real nu,
                              real diffusivity_fine,
                              unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFThS7<<< grid.grid, grid.threads >>> (  DC,
                                                DF,
                                                DD7C,
                                                DD7F,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                nu,
                                                diffusivity_fine);
      getLastCudaError("scaleCFThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS27( real* DC,
                              real* DF,
                              real* DD27C,
                              real* DD27F,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posCSWB,
                              unsigned int* posFSWB,
                              unsigned int kCF,
                              real nu,
                              real diffusivity_fine,
                              unsigned int numberOfThreads,
							  OffCF offCF)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

      scaleCFThS27<<< grid.grid, grid.threads >>> ( DC,
                                                DF,
                                                DD27C,
                                                DD27F,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posCSWB,
                                                posFSWB,
                                                kCF,
                                                nu,
                                                diffusivity_fine,
										        offCF);
      getLastCudaError("scaleCFThS27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC27( real* DC,
                           real* DF,
                           unsigned int* neighborCX,
                           unsigned int* neighborCY,
                           unsigned int* neighborCZ,
                           unsigned int* neighborFX,
                           unsigned int* neighborFY,
                           unsigned int* neighborFZ,
                           unsigned int size_MatC,
                           unsigned int size_MatF,
                           bool isEvenTimestep,
                           unsigned int* posC,
                           unsigned int* posFSWB,
                           unsigned int kFC,
                           real omCoarse,
                           real omFine,
                           real nu,
                           unsigned int nxC,
                           unsigned int nyC,
                           unsigned int nxF,
                           unsigned int nyF,
                           unsigned int numberOfThreads)
{
   
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC27<<< grid.grid, grid.threads >>> ( DC,
                                             DF,
                                             neighborCX,
                                             neighborCY,
                                             neighborCZ,
                                             neighborFX,
                                             neighborFY,
                                             neighborFZ,
                                             size_MatC,
                                             size_MatF,
                                             isEvenTimestep,
                                             posC,
                                             posFSWB,
                                             kFC,
                                             omCoarse,
                                             omFine,
                                             nu,
                                             nxC,
                                             nyC,
                                             nxF,
                                             nyF);
      getLastCudaError("scaleFC27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCEff27(real* DC,
                             real* DF,
                             unsigned int* neighborCX,
                             unsigned int* neighborCY,
                             unsigned int* neighborCZ,
                             unsigned int* neighborFX,
                             unsigned int* neighborFY,
                             unsigned int* neighborFZ,
                             unsigned int size_MatC,
                             unsigned int size_MatF,
                             bool isEvenTimestep,
                             unsigned int* posC,
                             unsigned int* posFSWB,
                             unsigned int kFC,
                             real omCoarse,
                             real omFine,
                             real nu,
                             unsigned int nxC,
                             unsigned int nyC,
                             unsigned int nxF,
                             unsigned int nyF,
                             unsigned int numberOfThreads,
                             OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCEff27<<< grid.grid, grid.threads >>> ( DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posC,
                                                posFSWB,
                                                kFC,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offFC);
      getLastCudaError("scaleFCEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCLast27(real* DC,
                              real* DF,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real omCoarse,
                              real omFine,
                              real nu,
                              unsigned int nxC,
                              unsigned int nyC,
                              unsigned int nxF,
                              unsigned int nyF,
                              unsigned int numberOfThreads,
                              OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCLast27<<< grid.grid, grid.threads >>> (DC,
                                                DF,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posC,
                                                posFSWB,
                                                kFC,
                                                omCoarse,
                                                omFine,
                                                nu,
                                                nxC,
                                                nyC,
                                                nxF,
                                                nyF,
                                                offFC);
      getLastCudaError("Kernel execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCpress27(real* DC,
                              real* DF,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real omCoarse,
                              real omFine,
                              real nu,
                              unsigned int nxC,
                              unsigned int nyC,
                              unsigned int nxF,
                              unsigned int nyF,
                              unsigned int numberOfThreads,
                              OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCpress27<<< grid.grid, grid.threads >>> (  DC,
                                                   DF,
                                                   neighborCX,
                                                   neighborCY,
                                                   neighborCZ,
                                                   neighborFX,
                                                   neighborFY,
                                                   neighborFZ,
                                                   size_MatC,
                                                   size_MatF,
                                                   isEvenTimestep,
                                                   posC,
                                                   posFSWB,
                                                   kFC,
                                                   omCoarse,
                                                   omFine,
                                                   nu,
                                                   nxC,
                                                   nyC,
                                                   nxF,
                                                   nyF,
                                                   offFC);
      getLastCudaError("scaleFCpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_27(real* DC,
                              real* DF,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real omCoarse,
                              real omFine,
                              real nu,
                              unsigned int nxC,
                              unsigned int nyC,
                              unsigned int nxF,
                              unsigned int nyF,
                              unsigned int numberOfThreads,
                              OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_Fix_27<<< grid.grid, grid.threads >>> (  DC,
                                                   DF,
                                                   neighborCX,
                                                   neighborCY,
                                                   neighborCZ,
                                                   neighborFX,
                                                   neighborFY,
                                                   neighborFZ,
                                                   size_MatC,
                                                   size_MatF,
                                                   isEvenTimestep,
                                                   posC,
                                                   posFSWB,
                                                   kFC,
                                                   omCoarse,
                                                   omFine,
                                                   nu,
                                                   nxC,
                                                   nyC,
                                                   nxF,
                                                   nyF,
                                                   offFC);
      getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_comp_27(  real* DC,
									  real* DF,
									  unsigned int* neighborCX,
									  unsigned int* neighborCY,
									  unsigned int* neighborCZ,
									  unsigned int* neighborFX,
									  unsigned int* neighborFY,
									  unsigned int* neighborFZ,
									  unsigned int size_MatC,
									  unsigned int size_MatF,
									  bool isEvenTimestep,
									  unsigned int* posC,
									  unsigned int* posFSWB,
									  unsigned int kFC,
									  real omCoarse,
									  real omFine,
									  real nu,
									  unsigned int nxC,
									  unsigned int nyC,
									  unsigned int nxF,
									  unsigned int nyF,
									  unsigned int numberOfThreads,
									  OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_Fix_comp_27<<< grid.grid, grid.threads >>> ( DC,
													   DF,
													   neighborCX,
													   neighborCY,
													   neighborCZ,
													   neighborFX,
													   neighborFY,
													   neighborFZ,
													   size_MatC,
													   size_MatF,
													   isEvenTimestep,
													   posC,
													   posFSWB,
													   kFC,
													   omCoarse,
													   omFine,
													   nu,
													   nxC,
													   nyC,
													   nxF,
													   nyF,
													   offFC);
      getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_0817_comp_27( real* DC,
									  real* DF,
									  unsigned int* neighborCX,
									  unsigned int* neighborCY,
									  unsigned int* neighborCZ,
									  unsigned int* neighborFX,
									  unsigned int* neighborFY,
									  unsigned int* neighborFZ,
									  unsigned int size_MatC,
									  unsigned int size_MatF,
									  bool isEvenTimestep,
									  unsigned int* posC,
									  unsigned int* posFSWB,
									  unsigned int kFC,
									  real omCoarse,
									  real omFine,
									  real nu,
									  unsigned int nxC,
									  unsigned int nyC,
									  unsigned int nxF,
									  unsigned int nyF,
									  unsigned int numberOfThreads,
									  OffFC offFC,
                             CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_0817_comp_27<<< grid.grid, grid.threads, 0, stream >>> (DC,
													   DF,
													   neighborCX,
													   neighborCY,
													   neighborCZ,
													   neighborFX,
													   neighborFY,
													   neighborFZ,
													   size_MatC,
													   size_MatF,
													   isEvenTimestep,
													   posC,
													   posFSWB,
													   kFC,
													   omCoarse,
													   omFine,
													   nu,
													   nxC,
													   nyC,
													   nxF,
													   nyF,
													   offFC);
      getLastCudaError("scaleFC_0817_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3_2018( real* DC,
										   real* DF,
										   real* G6,
										   unsigned int* neighborCX,
										   unsigned int* neighborCY,
										   unsigned int* neighborCZ,
										   unsigned int* neighborFX,
										   unsigned int* neighborFY,
										   unsigned int* neighborFZ,
										   unsigned int size_MatC,
										   unsigned int size_MatF,
										   bool isEvenTimestep,
										   unsigned int* posC,
										   unsigned int* posFSWB,
										   unsigned int kFC,
										   real omCoarse,
										   real omFine,
										   real nu,
										   unsigned int nxC,
										   unsigned int nyC,
										   unsigned int nxF,
										   unsigned int nyF,
										   unsigned int numberOfThreads,
										   OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

     scaleFC_comp_D3Q27F3_2018 <<< grid.grid, grid.threads >>> (DC,
															DF,
															G6,
															neighborCX,
															neighborCY,
															neighborCZ,
															neighborFX,
															neighborFY,
															neighborFZ,
															size_MatC,
															size_MatF,
															isEvenTimestep,
															posC,
															posFSWB,
															kFC,
															omCoarse,
															omFine,
															nu,
															nxC,
															nyC,
															nxF,
															nyF,
															offFC);
      getLastCudaError("scaleFC_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3( real* DC,
									  real* DF,
									  real* G6,
									  unsigned int* neighborCX,
									  unsigned int* neighborCY,
									  unsigned int* neighborCZ,
									  unsigned int* neighborFX,
									  unsigned int* neighborFY,
									  unsigned int* neighborFZ,
									  unsigned int size_MatC,
									  unsigned int size_MatF,
									  bool isEvenTimestep,
									  unsigned int* posC,
									  unsigned int* posFSWB,
									  unsigned int kFC,
									  real omCoarse,
									  real omFine,
									  real nu,
									  unsigned int nxC,
									  unsigned int nyC,
									  unsigned int nxF,
									  unsigned int nyF,
									  unsigned int numberOfThreads,
									  OffFC offFC,
                             CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

     scaleFC_comp_D3Q27F3 <<< grid.grid, grid.threads, 0, stream >>> (DC,
													   DF,
													   G6,
													   neighborCX,
													   neighborCY,
													   neighborCZ,
													   neighborFX,
													   neighborFY,
													   neighborFZ,
													   size_MatC,
													   size_MatF,
													   isEvenTimestep,
													   posC,
													   posFSWB,
													   kFC,
													   omCoarse,
													   omFine,
													   nu,
													   nxC,
													   nyC,
													   nxF,
													   nyF,
													   offFC);
      getLastCudaError("scaleFC_0817_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_staggered_time_comp_27(   real* DC,
												  real* DF,
												  unsigned int* neighborCX,
												  unsigned int* neighborCY,
												  unsigned int* neighborCZ,
												  unsigned int* neighborFX,
												  unsigned int* neighborFY,
												  unsigned int* neighborFZ,
												  unsigned int size_MatC,
												  unsigned int size_MatF,
												  bool isEvenTimestep,
												  unsigned int* posC,
												  unsigned int* posFSWB,
												  unsigned int kFC,
												  real omCoarse,
												  real omFine,
												  real nu,
												  unsigned int nxC,
												  unsigned int nyC,
												  unsigned int nxF,
												  unsigned int nyF,
												  unsigned int numberOfThreads,
												  OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_staggered_time_comp_27<<< grid.grid, grid.threads >>> (  DC,
																   DF,
																   neighborCX,
																   neighborCY,
																   neighborCZ,
																   neighborFX,
																   neighborFY,
																   neighborFZ,
																   size_MatC,
																   size_MatF,
																   isEvenTimestep,
																   posC,
																   posFSWB,
																   kFC,
																   omCoarse,
																   omFine,
																   nu,
																   nxC,
																   nyC,
																   nxF,
																   nyF,
																   offFC);
      getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_comp_27(real* DC,
									  real* DF,
									  unsigned int* neighborCX,
									  unsigned int* neighborCY,
									  unsigned int* neighborCZ,
									  unsigned int* neighborFX,
									  unsigned int* neighborFY,
									  unsigned int* neighborFZ,
									  unsigned int size_MatC,
									  unsigned int size_MatF,
									  bool isEvenTimestep,
									  unsigned int* posC,
									  unsigned int* posFSWB,
									  unsigned int kFC,
									  real omCoarse,
									  real omFine,
									  real nu,
									  unsigned int nxC,
									  unsigned int nyC,
									  unsigned int nxF,
									  unsigned int nyF,
									  unsigned int numberOfThreads,
									  OffFC offFC,
                             CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_RhoSq_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
													   DC,
													   DF,
													   neighborCX,
													   neighborCY,
													   neighborCZ,
													   neighborFX,
													   neighborFY,
													   neighborFZ,
													   size_MatC,
													   size_MatF,
													   isEvenTimestep,
													   posC,
													   posFSWB,
													   kFC,
													   omCoarse,
													   omFine,
													   nu,
													   nxC,
													   nyC,
													   nxF,
													   nyF,
													   offFC);
      getLastCudaError("scaleFC_RhoSq_27 execution failed");
}

//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_3rdMom_comp_27( real* DC,
											  real* DF,
											  unsigned int* neighborCX,
											  unsigned int* neighborCY,
											  unsigned int* neighborCZ,
											  unsigned int* neighborFX,
											  unsigned int* neighborFY,
											  unsigned int* neighborFZ,
											  unsigned int size_MatC,
											  unsigned int size_MatF,
											  bool isEvenTimestep,
											  unsigned int* posC,
											  unsigned int* posFSWB,
											  unsigned int kFC,
											  real omCoarse,
											  real omFine,
											  real nu,
											  unsigned int nxC,
											  unsigned int nyC,
											  unsigned int nxF,
											  unsigned int nyF,
											  unsigned int numberOfThreads,
											  OffFC offFC,
                                   CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_RhoSq_3rdMom_comp_27<<< grid.grid, grid.threads, 0, stream >>>(DC,
															  DF,
															  neighborCX,
															  neighborCY,
															  neighborCZ,
															  neighborFX,
															  neighborFY,
															  neighborFZ,
															  size_MatC,
															  size_MatF,
															  isEvenTimestep,
															  posC,
															  posFSWB,
															  kFC,
															  omCoarse,
															  omFine,
															  nu,
															  nxC,
															  nyC,
															  nxF,
															  nyF,
															  offFC);
      getLastCudaError("scaleFC_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_AA2016_comp_27( real* DC,
										real* DF,
										unsigned int* neighborCX,
										unsigned int* neighborCY,
										unsigned int* neighborCZ,
										unsigned int* neighborFX,
										unsigned int* neighborFY,
										unsigned int* neighborFZ,
										unsigned int size_MatC,
										unsigned int size_MatF,
										bool isEvenTimestep,
										unsigned int* posC,
										unsigned int* posFSWB,
										unsigned int kFC,
										real omCoarse,
										real omFine,
										real nu,
										unsigned int nxC,
										unsigned int nyC,
										unsigned int nxF,
										unsigned int nyF,
										unsigned int numberOfThreads,
										OffFC offFC,
                              CUstream_st *stream)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_AA2016_comp_27<<< grid.grid, grid.threads, 0, stream >>>(DC,
														DF,
														neighborCX,
														neighborCY,
														neighborCZ,
														neighborFX,
														neighborFY,
														neighborFZ,
														size_MatC,
														size_MatF,
														isEvenTimestep,
														posC,
														posFSWB,
														kFC,
														omCoarse,
														omFine,
														nu,
														nxC,
														nyC,
														nxF,
														nyF,
														offFC);
      getLastCudaError("scaleFC_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_NSPress_27(real* DC,
								  real* DF,
								  unsigned int* neighborCX,
								  unsigned int* neighborCY,
								  unsigned int* neighborCZ,
								  unsigned int* neighborFX,
								  unsigned int* neighborFY,
								  unsigned int* neighborFZ,
								  unsigned int size_MatC,
								  unsigned int size_MatF,
								  bool isEvenTimestep,
								  unsigned int* posC,
								  unsigned int* posFSWB,
								  unsigned int kFC,
								  real omCoarse,
								  real omFine,
								  real nu,
								  unsigned int nxC,
								  unsigned int nyC,
								  unsigned int nxF,
								  unsigned int nyF,
								  unsigned int numberOfThreads,
								  OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFC_NSPress_27<<< grid.grid, grid.threads >>> (  DC,
													   DF,
													   neighborCX,
													   neighborCY,
													   neighborCZ,
													   neighborFX,
													   neighborFY,
													   neighborFZ,
													   size_MatC,
													   size_MatF,
													   isEvenTimestep,
													   posC,
													   posFSWB,
													   kFC,
													   omCoarse,
													   omFine,
													   nu,
													   nxC,
													   nyC,
													   nxF,
													   nyF,
													   offFC);
      getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThSMG7(real* DC,
                              real* DF,
                              real* DD7C,
                              real* DD7F,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real nu,
                              real diffusivity_coarse,
                              unsigned int numberOfThreads,
                              OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCThSMG7<<< grid.grid, grid.threads >>>( DC,
                                                DF,
                                                DD7C,
                                                DD7F,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posC,
                                                posFSWB,
                                                kFC,
                                                nu,
                                                diffusivity_coarse,
                                                offFC);
      getLastCudaError("scaleFCThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS7(  real* DC,
                              real* DF,
                              real* DD7C,
                              real* DD7F,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real nu,
                              real diffusivity_coarse,
                              unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCThS7<<< grid.grid, grid.threads >>>(DC,
                                             DF,
                                             DD7C,
                                             DD7F,
                                             neighborCX,
                                             neighborCY,
                                             neighborCZ,
                                             neighborFX,
                                             neighborFY,
                                             neighborFZ,
                                             size_MatC,
                                             size_MatF,
                                             isEvenTimestep,
                                             posC,
                                             posFSWB,
                                             kFC,
                                             nu,
                                             diffusivity_coarse);
      getLastCudaError("scaleFCThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS27( real* DC,
                              real* DF,
                              real* DD27C,
                              real* DD27F,
                              unsigned int* neighborCX,
                              unsigned int* neighborCY,
                              unsigned int* neighborCZ,
                              unsigned int* neighborFX,
                              unsigned int* neighborFY,
                              unsigned int* neighborFZ,
                              unsigned int size_MatC,
                              unsigned int size_MatF,
                              bool isEvenTimestep,
                              unsigned int* posC,
                              unsigned int* posFSWB,
                              unsigned int kFC,
                              real nu,
                              real diffusivity_coarse,
                              unsigned int numberOfThreads,
							  OffFC offFC)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

      scaleFCThS27<<< grid.grid, grid.threads >>>(  DC,
                                                DF,
                                                DD27C,
                                                DD27F,
                                                neighborCX,
                                                neighborCY,
                                                neighborCZ,
                                                neighborFX,
                                                neighborFY,
                                                neighborFZ,
                                                size_MatC,
                                                size_MatF,
                                                isEvenTimestep,
                                                posC,
                                                posFSWB,
                                                kFC,
                                                nu,
                                                diffusivity_coarse,
												offFC);
      getLastCudaError("scaleFCThS27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void DragLiftPostD27(real* DD,
								int* k_Q,
								real* QQ,
								int numberOfBCnodes,
								double *DragX,
								double *DragY,
								double *DragZ,
								unsigned int* neighborX,
								unsigned int* neighborY,
								unsigned int* neighborZ,
								unsigned int size_Mat,
								bool isEvenTimestep,
								unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

	DragLiftPost27<<< grid.grid, grid.threads >>>(DD,
										k_Q,
										QQ,
										numberOfBCnodes,
										DragX,
										DragY,
										DragZ,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat,
										isEvenTimestep);
	getLastCudaError("DragLift27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void DragLiftPreD27( real* DD,
								int* k_Q,
								real* QQ,
								int numberOfBCnodes,
								double *DragX,
								double *DragY,
								double *DragZ,
								unsigned int* neighborX,
								unsigned int* neighborY,
								unsigned int* neighborZ,
								unsigned int size_Mat,
								bool isEvenTimestep,
								unsigned int numberOfThreads)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

	DragLiftPre27<<< grid.grid, grid.threads >>>( DD,
										k_Q,
										QQ,
										numberOfBCnodes,
										DragX,
										DragY,
										DragZ,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat,
										isEvenTimestep);
	getLastCudaError("DragLift27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcCPtop27(real* DD,
							int* cpIndex,
							int nonCp,
							double *cpPress,
							unsigned int* neighborX,
							unsigned int* neighborY,
							unsigned int* neighborZ,
							unsigned int size_Mat,
							bool isEvenTimestep,
							unsigned int numberOfThreads)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

	CalcCP27<<< grid.grid, grid.threads >>>(DD,
								  cpIndex,
								  nonCp,
								  cpPress,
								  neighborX,
								  neighborY,
								  neighborZ,
								  size_Mat,
								  isEvenTimestep);
	getLastCudaError("CalcCP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcCPbottom27( real* DD,
								int* cpIndex,
								int nonCp,
								double *cpPress,
								unsigned int* neighborX,
								unsigned int* neighborY,
								unsigned int* neighborZ,
								unsigned int size_Mat,
								bool isEvenTimestep,
								unsigned int numberOfThreads)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

	CalcCP27<<< grid.grid, grid.threads >>>(DD,
								  cpIndex,
								  nonCp,
								  cpPress,
								  neighborX,
								  neighborY,
								  neighborZ,
								  size_Mat,
								  isEvenTimestep);
	getLastCudaError("CalcCP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void GetSendFsPreDev27(real* DD,
								  real* bufferFs,
								  int* sendIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat,
								  bool isEvenTimestep,
								  unsigned int numberOfThreads,
								  cudaStream_t stream)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	getSendFsPre27<<< grid.grid, grid.threads, 0, stream >>>(DD,
										bufferFs,
										sendIndex,
										buffmax,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat,
										isEvenTimestep);
	getLastCudaError("getSendFsPre27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void GetSendFsPostDev27(real* DD,
								   real* bufferFs,
								   int* sendIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep,
								   unsigned int numberOfThreads,
								   cudaStream_t stream)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	getSendFsPost27<<< grid.grid, grid.threads, 0, stream >>>(DD,
										 bufferFs,
										 sendIndex,
										 buffmax,
										 neighborX,
										 neighborY,
										 neighborZ,
										 size_Mat,
										 isEvenTimestep);
	getLastCudaError("getSendFsPost27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void SetRecvFsPreDev27(real* DD,
								  real* bufferFs,
								  int* recvIndex,
								  int buffmax,
								  unsigned int* neighborX,
								  unsigned int* neighborY,
								  unsigned int* neighborZ,
								  unsigned int size_Mat,
								  bool isEvenTimestep,
								  unsigned int numberOfThreads,
	                              cudaStream_t stream)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	setRecvFsPre27<<< grid.grid, grid.threads, 0, stream >>>(DD,
										bufferFs,
										recvIndex,
										buffmax,
										neighborX,
										neighborY,
										neighborZ,
										size_Mat,
										isEvenTimestep);
	getLastCudaError("setRecvFsPre27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void SetRecvFsPostDev27(real* DD,
								   real* bufferFs,
								   int* recvIndex,
								   int buffmax,
								   unsigned int* neighborX,
								   unsigned int* neighborY,
								   unsigned int* neighborZ,
								   unsigned int size_Mat,
								   bool isEvenTimestep,
	                               unsigned int numberOfThreads,
	                               cudaStream_t stream)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	setRecvFsPost27<<< grid.grid, grid.threads, 0, stream >>>(DD,
										 bufferFs,
										 recvIndex,
										 buffmax,
										 neighborX,
										 neighborY,
										 neighborZ,
										 size_Mat,
										 isEvenTimestep);
	getLastCudaError("setRecvFsPost27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void getSendGsDevF3(
	real* G6,
	real* bufferGs,
	int* sendIndex,
	int buffmax,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	bool isEvenTimestep,
	unsigned int numberOfThreads)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	getSendGsF3 <<< grid.grid, grid.threads >>> (
		G6,
		bufferGs,
		sendIndex,
		buffmax,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("getSendGsF3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void setRecvGsDevF3(
	real* G6,
	real* bufferGs,
	int* recvIndex,
	int buffmax,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	bool isEvenTimestep,
	unsigned int numberOfThreads)
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

	setRecvGsF3 <<< grid.grid, grid.threads >>> (
		G6,
		bufferGs,
		recvIndex,
		buffmax,
		neighborX,
		neighborY,
		neighborZ,
		size_Mat,
		isEvenTimestep);
	getLastCudaError("setRecvGsF3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void WallFuncDev27(unsigned int numberOfThreads,
							  real* vx,
							  real* vy,
							  real* vz,
							  real* DD,
							  int* k_Q,
							  real* QQ,
							  unsigned int numberOfBCnodes,
							  real om1,
							  unsigned int* neighborX,
							  unsigned int* neighborY,
							  unsigned int* neighborZ,
							  unsigned int size_Mat,
							  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

      WallFunction27<<< grid.grid, grid.threads >>> (
											  vx,
											  vy,
											  vz,
											  DD,
											  k_Q,
											  QQ,
											  numberOfBCnodes,
											  om1,
											  neighborX,
											  neighborY,
											  neighborZ,
											  size_Mat,
											  isEvenTimestep);
      getLastCudaError("WallFunction27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void SetOutputWallVelocitySP27(unsigned int numberOfThreads,
										  real* vxD,
										  real* vyD,
										  real* vzD,
										  real* vxWall,
										  real* vyWall,
										  real* vzWall,
										  int numberOfWallNodes,
										  int* kWallNodes,
										  real* rhoD,
										  real* pressD,
										  unsigned int* geoD,
										  unsigned int* neighborX,
										  unsigned int* neighborY,
										  unsigned int* neighborZ,
										  unsigned int size_Mat,
										  real* DD,
										  bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfWallNodes);

      LBSetOutputWallVelocitySP27<<< grid.grid, grid.threads >>> (	vxD,
															vyD,
															vzD,
															vxWall,
															vyWall,
															vzWall,
															numberOfWallNodes,
															kWallNodes,
															rhoD,
															pressD,
															geoD,
															neighborX,
															neighborY,
															neighborZ,
															size_Mat,
															DD,
															isEvenTimestep);
      getLastCudaError("LBSetOutputWallVelocitySP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void GetVelotoForce27(unsigned int numberOfThreads,
								 real* DD,
								 int* bcIndex,
								 int nonAtBC,
								 real* Vx,
								 real* Vy,
								 real* Vz,
								 unsigned int* neighborX,
								 unsigned int* neighborY,
								 unsigned int* neighborZ,
								 unsigned int size_Mat,
								 bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonAtBC);

      GetVeloforForcing27<<< grid.grid, grid.threads >>> (DD,
												bcIndex,
												nonAtBC,
												Vx,
												Vy,
												Vz,
												neighborX,
												neighborY,
												neighborZ,
												size_Mat,
												isEvenTimestep);
      getLastCudaError("GetVeloforForcing27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void InitParticlesDevice(real* coordX,
									real* coordY,
									real* coordZ,
									real* coordParticleXlocal,
									real* coordParticleYlocal,
									real* coordParticleZlocal,
									real* coordParticleXglobal,
									real* coordParticleYglobal,
									real* coordParticleZglobal,
									real* veloParticleX,
									real* veloParticleY,
									real* veloParticleZ,
									real* randArray,
									unsigned int* particleID,
									unsigned int* cellBaseID,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int* neighborWSB,
									int level,
									unsigned int numberOfParticles,
									unsigned int size_Mat,
									unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfParticles);

   InitParticles<<< grid.grid, grid.threads >>> (coordX,
										coordY,
										coordZ,
										coordParticleXlocal,
										coordParticleYlocal,
										coordParticleZlocal,
										coordParticleXglobal,
										coordParticleYglobal,
										coordParticleZglobal,
										veloParticleX,
										veloParticleY,
										veloParticleZ,
										randArray,
										particleID,
										cellBaseID,
										bcMatD,
										neighborX,
										neighborY,
										neighborZ,
										neighborWSB,
										level,
										numberOfParticles,
										size_Mat);
      getLastCudaError("InitParticles execution failed");
}
//////////////////////////////////////////////////////////////////////////
void MoveParticlesDevice(real* coordX,
									real* coordY,
									real* coordZ,
									real* coordParticleXlocal,
									real* coordParticleYlocal,
									real* coordParticleZlocal,
									real* coordParticleXglobal,
									real* coordParticleYglobal,
									real* coordParticleZglobal,
									real* veloParticleX,
									real* veloParticleY,
									real* veloParticleZ,
									real* DD,
									real  omega,
									unsigned int* particleID,
									unsigned int* cellBaseID,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int* neighborWSB,
							        int level,
									unsigned int timestep,
									unsigned int numberOfTimesteps,
									unsigned int numberOfParticles,
									unsigned int size_Mat,
									unsigned int numberOfThreads,
									bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfParticles);

   MoveParticles<<< grid.grid, grid.threads >>> (coordX,
										coordY,
										coordZ,
										coordParticleXlocal,
										coordParticleYlocal,
										coordParticleZlocal,
										coordParticleXglobal,
										coordParticleYglobal,
										coordParticleZglobal,
										veloParticleX,
										veloParticleY,
										veloParticleZ,
										DD,
										omega,
										particleID,
										cellBaseID,
										bcMatD,
										neighborX,
										neighborY,
										neighborZ,
										neighborWSB,
										level,
										timestep,
										numberOfTimesteps,
										numberOfParticles,
										size_Mat,
										isEvenTimestep);
      getLastCudaError("MoveParticles execution failed");
}
//////////////////////////////////////////////////////////////////////////
void initRandomDevice(curandState* state,
								 unsigned int size_Mat,
								 unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);
   initRandom<<< grid.grid, grid.threads >>> (state);
   getLastCudaError("initRandom execution failed");
}
//////////////////////////////////////////////////////////////////////////
void generateRandomValuesDevice( curandState* state,
											unsigned int size_Mat,
											real* randArray,
											unsigned int numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);
   generateRandomValues<<< grid.grid, grid.threads >>> (state,randArray);
   getLastCudaError("generateRandomValues execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcTurbulenceIntensityDevice(
   real* vxx,
   real* vyy,
   real* vzz,
   real* vxy,
   real* vxz,
   real* vyz,
   real* vx_mean,
   real* vy_mean,
   real* vz_mean,
   real* DD,
   uint* typeOfGridNode,
   unsigned int* neighborX,
   unsigned int* neighborY,
   unsigned int* neighborZ,
   unsigned int size_Mat,
   bool isEvenTimestep,
   uint numberOfThreads)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, size_Mat);
   CalcTurbulenceIntensity<<<grid.grid, grid.threads>>>(
     vxx,
     vyy,
     vzz,
	 vxy,
     vxz,
     vyz,
     vx_mean,
     vy_mean,
     vz_mean,
     DD,
     typeOfGridNode,
     neighborX,
     neighborY,
     neighborZ,
     size_Mat,
     isEvenTimestep);

   getLastCudaError("CalcTurbulenceIntensity execution failed");
}













