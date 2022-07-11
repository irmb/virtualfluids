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
extern "C" void KernelCas27( unsigned int grid_nx,
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
extern "C" void KernelCasSP27( unsigned int numberOfThreads,
                               real s9,
                               unsigned int* bcMatD,
                               unsigned int* neighborX,
                               unsigned int* neighborY,
                               unsigned int* neighborZ,
                               real* DD,
                               int size_Mat,
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

      LB_Kernel_Casc_SP_27<<< grid, threads >>>(s9,
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
extern "C" void KernelCasSPMS27( unsigned int numberOfThreads,
                                 real s9,
                                 unsigned int* bcMatD,
                                 unsigned int* neighborX,
                                 unsigned int* neighborY,
                                 unsigned int* neighborZ,
                                 real* DD,
                                 int size_Mat,
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

      LB_Kernel_Casc_SP_MS_27<<< grid, threads >>>(s9,
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
extern "C" void KernelCasSPMSOHM27( unsigned int numberOfThreads,
                                    real s9,
                                    unsigned int* bcMatD,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    real* DD,
                                    int size_Mat,
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

      LB_Kernel_Casc_SP_MS_OHM_27<<< grid, threads >>>(  s9,
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
extern "C" void KernelKumCompSRTSP27(
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

   LB_Kernel_Kum_New_Comp_SRT_SP_27 <<< grid, threads >>>(
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
extern "C" void KernelKum1hSP27(    unsigned int numberOfThreads,
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

		LB_Kernel_Kum_1h_SP_27<<< grid, threads >>>(omega,
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
extern "C" void KernelCascadeSP27(  unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
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

		LB_Kernel_Cascade_SP_27<<< grid, threads >>>(s9,
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
extern "C" void KernelKumNewSP27(   unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
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

		LB_Kernel_Kum_New_SP_27<<< grid, threads >>>(s9,
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
extern "C" void KernelKumNewCompSP27(unsigned int numberOfThreads,
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

		//LB_Kernel_Kum_New_Comp_SP_27<<< grid, threads >>>(	s9,
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
extern "C" void CumulantOnePreconditionedErrorDiffusionChimCompSP27(unsigned int numberOfThreads,
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	Cumulant_One_preconditioned_errorDiffusion_chim_Comp_SP_27 <<< grid, threads >>>(	s9,
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
extern "C" void CumulantOnePreconditionedChimCompSP27(  unsigned int numberOfThreads,
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	Cumulant_One_preconditioned_chim_Comp_SP_27 <<< grid, threads >>>(	s9,
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
extern "C" void CumulantOneChimCompSP27(unsigned int numberOfThreads,
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	Cumulant_One_chim_Comp_SP_27 <<< grid, threads >>>(	s9,
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
extern "C" void KernelKumIsoTestSP27(unsigned int numberOfThreads,
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

	LB_Kernel_Kum_IsoTest_SP_27<<< grid, threads >>>(s9,
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
extern "C" void KernelKumCompSP27(  unsigned int numberOfThreads,
									real s9,
									unsigned int* bcMatD,
									unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									real* DD,
									int size_Mat,
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

		LB_Kernel_Kum_Comp_SP_27<<< grid, threads >>>(s9,
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
extern "C" void KernelPMCumOneCompSP27(unsigned int numberOfThreads,
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_PM_Cum_One_Comp_SP_27 <<< grid, threads >>>(omega,
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
extern "C" void KernelWaleBySoniMalavCumAA2016CompSP27(
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
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Kernel_WaleBySoniMalav_Cum_AA2016_Comp_SP_27 << < grid, threads >> >(
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
extern "C" void KernelADincomp7(   unsigned int numberOfThreads,
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

      LB_Kernel_AD_Incomp_7<<< grid, threads >>>( diffusivity,
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
extern "C" void KernelADincomp27( unsigned int numberOfThreads,
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

      LB_Kernel_AD_Incomp_27<<< grid, threads >>>( diffusivity,
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
extern "C" void Init27( int myid,
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
extern "C" void InitNonEqPartSP27( unsigned int numberOfThreads,
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

      LBInitNonEqPartSP27<<< grid, threads >>>( neighborX,
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
extern "C" void InitThS7(     unsigned int numberOfThreads,
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

      InitAD7<<< grid, threads >>>( neighborX,
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
extern "C" void InitADDev27( unsigned int numberOfThreads,
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

      InitAD27<<< grid, threads >>>(neighborX,
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
extern "C" void PostProcessorF3_2018Fehlberg(
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
	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid>512)
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
	dim3 threads(numberOfThreads, 1, 1);

	  LB_PostProcessor_F3_2018_Fehlberg <<< grid, threads >>> (   omega,
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
extern "C" void CalcMac27( real* vxD,
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
extern "C" void CalcMacSP27( real* vxD,
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

      LBCalcMacSP27<<< grid, threads >>> (   vxD,
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
extern "C" void CalcMacCompSP27( real* vxD,
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
												 isEvenTimestep);
      getLastCudaError("LBCalcMacSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void CalcMacThS7(  real* Conc,
                              unsigned int* geoD,
                              unsigned int* neighborX,
                              unsigned int* neighborY,
                              unsigned int* neighborZ,
                              unsigned int size_Mat,
                              unsigned int numberOfThreads,
                              real* DD7,
                              bool isEvenTimestep)
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

      CalcConc7<<< grid, threads >>> (Conc,
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
extern "C" void PlaneConcThS7(real* Conc,
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
   int Grid = (numberOfPointskPC / numberOfThreads)+1;
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

      GetPlaneConc7<<< grid, threads >>> (	Conc,
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
extern "C" void PlaneConcThS27(real* Conc,
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
   int Grid = (numberOfPointskPC / numberOfThreads)+1;
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

      GetPlaneConc27<<< grid, threads >>> (	Conc,
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
extern "C" void CalcConcentration27( unsigned int numberOfThreads,
                                     real* Conc,
                                     unsigned int* geoD,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     unsigned int size_Mat,
                                     real* DD27,
                                     bool isEvenTimestep)
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

      CalcConc27<<< grid, threads >>> (  Conc,
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
extern "C" void CalcMedSP27(  real* vxD,
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

      LBCalcMedSP27<<< grid, threads >>> (   vxD,
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
extern "C" void CalcMedCompSP27(  real* vxD,
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

      LBCalcMedCompSP27<<< grid, threads >>> (   vxD,
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
extern "C" void CalcMedCompAD27(
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
	dim3 threads(numberOfThreads, 1, 1);

	LBCalcMedCompAD27 <<< grid, threads >>> (
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
extern "C" void CalcMacMedSP27(  real* vxD,
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

      LBCalcMacMedSP27<<< grid, threads >>> (   vxD,
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
extern "C" void ResetMedianValuesSP27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	unsigned int size_Mat,
	unsigned int numberOfThreads,
	bool isEvenTimestep)
{
	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid>512)
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
	dim3 threads(numberOfThreads, 1, 1);

	LBResetMedianValuesSP27 << < grid, threads >> > (
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
extern "C" void ResetMedianValuesAD27(
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
	dim3 threads(numberOfThreads, 1, 1);

	LBResetMedianValuesAD27 << < grid, threads >> > (
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
extern "C" void Calc2ndMomentsIncompSP27(real* kxyFromfcNEQ,
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

      LBCalc2ndMomentsIncompSP27<<< grid, threads >>> (  kxyFromfcNEQ,
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
extern "C" void Calc2ndMomentsCompSP27( real* kxyFromfcNEQ,
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

      LBCalc2ndMomentsCompSP27<<< grid, threads >>> (kxyFromfcNEQ,
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
extern "C" void Calc3rdMomentsIncompSP27(real* CUMbbb,
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

      LBCalc3rdMomentsIncompSP27<<< grid, threads >>> (  CUMbbb,
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
extern "C" void Calc3rdMomentsCompSP27( real* CUMbbb,
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

      LBCalc3rdMomentsCompSP27<<< grid, threads >>> (CUMbbb,
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
extern "C" void CalcHigherMomentsIncompSP27(real* CUMcbb,
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

      LBCalcHigherMomentsIncompSP27<<< grid, threads >>> (CUMcbb,
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
extern "C" void CalcHigherMomentsCompSP27(  real* CUMcbb,
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

      LBCalcHigherMomentsCompSP27<<< grid, threads >>> (  CUMcbb,
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
extern "C" void LBCalcMeasurePoints27(real* vxMP,
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
   int Grid = (numberOfPointskMP / numberOfThreads)+1;
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

      LBCalcMeasurePoints<<< grid, threads >>> (vxMP,
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
extern "C" void BcPress27( int nx,
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
extern "C" void BcVel27(int nx,
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
extern "C" void QADPressDev7( unsigned int numberOfThreads,
                              int nx,
                              int ny,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADPress7<<< gridQ, threads >>>( DD,
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
extern "C" void QADPressDev27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADPress27<<< gridQ, threads >>>(   DD,
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
extern "C" void QADPressNEQNeighborDev27(
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

   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

   QADPressNEQNeighbor27<<< gridQ, threads >>>(
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
extern "C" void QADVelDev7(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADVel7<<< gridQ, threads >>> (  
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
extern "C" void QADVelDev27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADVel27<<< gridQ, threads >>> ( DD,
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
extern "C" void QADDev7(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QAD7<<< gridQ, threads >>> (     DD,
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
extern "C" void FactorizedCentralMomentsAdvectionDiffusionDeviceKernel(
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
extern "C" void ADSlipVelDevComp(
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
	int Grid = (numberOfBCnodes / numberOfThreads) + 1;
	dim3 gridQ(Grid, 1, 1);
	dim3 threads(numberOfThreads, 1, 1);

	AD_SlipVelDeviceComp << < gridQ, threads >> > (
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

extern "C" void QADDirichletDev27( unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADDirichlet27<<< gridQ, threads >>> (
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
extern "C" void QADBBDev27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADBB27<<< gridQ, threads >>> (  DD,
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
extern "C" void QNoSlipADincompDev7(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QNoSlipADincomp7<<< gridQ, threads >>> (
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
extern "C" void QNoSlipADincompDev27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QNoSlipADincomp27<<< gridQ, threads >>> (
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
extern "C" void QADVeloIncompDev7( unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADVeloIncomp7<<< gridQ, threads >>> ( 
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
      getLastCudaError("QADVeloIncomp7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void QADVeloIncompDev27(   unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADVeloIncomp27<<< gridQ, threads >>> (
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
extern "C" void QADPressIncompDev7( unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADPressIncomp7<<< gridQ, threads >>>(
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
extern "C" void QADPressIncompDev27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QADPressIncomp27<<< gridQ, threads >>>(
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
      getLastCudaError("QADPressIncomp27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void QDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QDevCompThinWalls27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

   QDeviceCompThinWallsPartOne27 <<< gridQ, threads >>> (DD,
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

   QThinWallsPartTwo27 <<< gridQ, threads >>> ( DD,
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
extern "C" void QDev3rdMomentsComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QDevIncompHighNu27( unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QDeviceIncompHighNu27<<< gridQ, threads >>> (
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
extern "C" void QDevCompHighNu27(   unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QDeviceCompHighNu27<<< gridQ, threads >>> (
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
extern "C" void QVelDevicePlainBB27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QVelDeviceCouette27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVelDevCouette27<<< gridQ, threads >>> ( vx,
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
extern "C" void QVelDevice1h27(   unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVelDev1h27<<< gridQ, threads >>> (nx,
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
extern "C" void QVelDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QVelDevCompPlusSlip27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVelDeviceCompPlusSlip27<<< gridQ, threads >>> (
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
extern "C" void QVelDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QVelDevCompThinWalls27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

   QVelDeviceCompThinWallsPartOne27<<< gridQ, threads >>> (vx,
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

	QThinWallsPartTwo27 <<< gridQ, threads >>> (
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

extern "C" void QVelDevCompZeroPress27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QVelDevIncompHighNu27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVelDeviceIncompHighNu27<<< gridQ, threads >>> (
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
extern "C" void QVelDevCompHighNu27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVelDeviceCompHighNu27<<< gridQ, threads >>> (
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
extern "C" void QVeloDevEQ27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QVeloDeviceEQ27<<< gridQ, threads >>> (VeloX,
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
extern "C" void QVeloStreetDevEQ27(
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
	int Grid = (numberOfStreetNodes / numberOfThreads) + 1;
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
	dim3 gridQ(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	QVeloStreetDeviceEQ27 << < gridQ, threads >> > (
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
extern "C" void QSlipDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QSlipDevCompTurbulentViscosity27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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

extern "C" void QSlipDevComp27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QSlipGeomDevComp27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QSlipGeomDeviceComp27<<< gridQ, threads >>> (DD,
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
extern "C" void QSlipNormDevComp27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QSlipNormDeviceComp27<<< gridQ, threads >>> (DD,
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
extern "C" void QStressDevComp27(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
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
extern "C" void BBStressDev27(Parameter *para,  QforBoundaryConditions* boundaryCondition, const int level)
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
extern "C" void QPressDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDevAntiBB27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

    QPressDeviceAntiBB27<<< gridQ, threads >>>( rhoBC,
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
extern "C" void QPressDevFixBackflow27( unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QPressDeviceFixBackflow27<<< gridQ, threads >>> (  rhoBC,
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
extern "C" void QPressDevDirDepBot27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QPressDeviceDirDepBot27<<< gridQ, threads >>> ( rhoBC,
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
extern "C" void QPressNoRhoDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
         parameterDevice->isEvenTimestep);
   getLastCudaError("QPressNoRhoDevice27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void QInflowScaleByPressDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDevOld27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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
extern "C" void QPressDevIncompNEQ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDevNEQ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDevEQZ27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDevZero27(unsigned int numberOfThreads,
                                real* DD,
                                int* k_Q,
                                unsigned int numberOfBCnodes,
                                unsigned int* neighborX,
                                unsigned int* neighborY,
                                unsigned int* neighborZ,
                                unsigned int size_Mat,
                                bool isEvenTimestep)
{
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QPressDeviceZero27<<< gridQ, threads >>> (DD,
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
extern "C" void QPressDevFake27(     unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      QPressDeviceFake27<<< gridQ, threads >>> (rhoBC,
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
extern "C" void BBDev27(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
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
extern "C" void QPressDev27_IntBB(  unsigned int numberOfThreads,
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
	int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

		QPressDevice27_IntBB<<< gridQ, threads >>> (rho,
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
extern "C" void PressSchlaffer27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      PressSchlaff27<<< gridQ, threads >>>(  rhoBC,
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
extern "C" void VelSchlaffer27(  unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      VelSchlaff27<<< gridQ, threads >>>( t,
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
   int Grid = (size_Prop / numberOfThreads)+1;
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

      PropellerBC<<< grid, threads >>>(neighborX,
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
extern "C" void ScaleCF27( real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF27<<< gridINT_CF, threads >>> ( DC,
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
extern "C" void ScaleCFEff27(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFEff27<<< gridINT_CF, threads >>> ( DC,
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
extern "C" void ScaleCFLast27(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFLast27<<< gridINT_CF, threads >>> (DC,
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
extern "C" void ScaleCFpress27(  real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFpress27<<< gridINT_CF, threads >>>(DC,
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
extern "C" void ScaleCF_Fix_27(  real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_Fix_27<<< gridINT_CF, threads >>>(DC,
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
extern "C" void ScaleCF_Fix_comp_27( real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_Fix_comp_27<<< gridINT_CF, threads >>>(   DC,
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
extern "C" void ScaleCF_0817_comp_27(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_0817_comp_27<<< gridINT_CF, threads, 0, stream >>>(  DC,
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
extern "C" void ScaleCF_comp_D3Q27F3_2018(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_comp_D3Q27F3_2018 <<< gridINT_CF, threads >>>(DC,
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
extern "C" void ScaleCF_comp_D3Q27F3(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_comp_D3Q27F3 <<< gridINT_CF, threads, 0, stream >>>( DC,
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
extern "C" void ScaleCF_staggered_time_comp_27(  real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_staggered_time_comp_27<<< gridINT_CF, threads >>>(    DC,
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
extern "C" void ScaleCF_RhoSq_comp_27(   real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_RhoSq_comp_27<<< gridINT_CF, threads, 0, stream >>>( DC,
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
extern "C" void ScaleCF_RhoSq_3rdMom_comp_27(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_RhoSq_3rdMom_comp_27<<< gridINT_CF, threads, 0, stream >>>(  DC,
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
extern "C" void ScaleCF_AA2016_comp_27(real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_AA2016_comp_27<<< gridINT_CF, threads, 0, stream >>>(DC,
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
extern "C" void ScaleCF_NSPress_27(  real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCF_NSPress_27<<< gridINT_CF, threads >>>(DC,
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
extern "C" void ScaleCFThSMG7(   real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFThSMG7<<< gridINT_CF, threads >>> (DC,
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
extern "C" void ScaleCFThS7(  real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFThS7<<< gridINT_CF, threads >>> (  DC,
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
extern "C" void ScaleCFThS27( real* DC,
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
   int Grid = (kCF / numberOfThreads)+1;
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
   dim3 gridINT_CF(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleCFThS27<<< gridINT_CF, threads >>> ( DC,
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
extern "C" void ScaleFC27( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC27<<< gridINT_FC, threads >>> ( DC,
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
extern "C" void ScaleFCEff27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCEff27<<< gridINT_FC, threads >>> ( DC,
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
extern "C" void ScaleFCLast27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCLast27<<< gridINT_FC, threads >>> (DC,
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
extern "C" void ScaleFCpress27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCpress27<<< gridINT_FC, threads >>> (  DC,
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
extern "C" void ScaleFC_Fix_27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_Fix_27<<< gridINT_FC, threads >>> (  DC,
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
extern "C" void ScaleFC_Fix_comp_27(  real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_Fix_comp_27<<< gridINT_FC, threads >>> ( DC,
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
extern "C" void ScaleFC_0817_comp_27( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_0817_comp_27<<< gridINT_FC, threads, 0, stream >>> (DC,
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
extern "C" void ScaleFC_comp_D3Q27F3_2018( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

     scaleFC_comp_D3Q27F3_2018 <<< gridINT_FC, threads >>> (DC,
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
extern "C" void ScaleFC_comp_D3Q27F3( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

     scaleFC_comp_D3Q27F3 <<< gridINT_FC, threads, 0, stream >>> (DC,
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
extern "C" void ScaleFC_staggered_time_comp_27(   real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_staggered_time_comp_27<<< gridINT_FC, threads >>> (  DC,
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
extern "C" void ScaleFC_RhoSq_comp_27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_RhoSq_comp_27<<<gridINT_FC, threads, 0, stream>>>(
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
extern "C" void ScaleFC_RhoSq_3rdMom_comp_27( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_RhoSq_3rdMom_comp_27<<< gridINT_FC, threads, 0, stream >>>(DC,
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
extern "C" void ScaleFC_AA2016_comp_27( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_AA2016_comp_27<<< gridINT_FC, threads, 0, stream >>>(DC,
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
extern "C" void ScaleFC_NSPress_27(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFC_NSPress_27<<< gridINT_FC, threads >>> (  DC,
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
extern "C" void ScaleFCThSMG7(real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCThSMG7<<< gridINT_FC, threads >>>( DC,
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
extern "C" void ScaleFCThS7(  real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCThS7<<< gridINT_FC, threads >>>(DC,
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
extern "C" void ScaleFCThS27( real* DC,
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
   int Grid = (kFC / numberOfThreads)+1;
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
   dim3 gridINT_FC(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

      scaleFCThS27<<< gridINT_FC, threads >>>(  DC,
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
extern "C" void DragLiftPostD27(real* DD,
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
	int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

	DragLiftPost27<<< grid, threads >>>(DD,
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
extern "C" void DragLiftPreD27( real* DD,
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
	int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

	DragLiftPre27<<< grid, threads >>>( DD,
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
extern "C" void CalcCPtop27(real* DD,
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
	int Grid = (nonCp / numberOfThreads)+1;
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

	CalcCP27<<< grid, threads >>>(DD,
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
extern "C" void CalcCPbottom27( real* DD,
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
	int Grid = (nonCp / numberOfThreads)+1;
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

	CalcCP27<<< grid, threads >>>(DD,
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
extern "C" void GetSendFsPreDev27(real* DD,
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

	getSendFsPre27<<< grid, threads, 0, stream >>>(DD,
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
extern "C" void GetSendFsPostDev27(real* DD,
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

	getSendFsPost27<<< grid, threads, 0, stream >>>(DD,
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
extern "C" void SetRecvFsPreDev27(real* DD,
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

	setRecvFsPre27<<< grid, threads, 0, stream >>>(DD,
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
extern "C" void SetRecvFsPostDev27(real* DD,
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

	setRecvFsPost27<<< grid, threads, 0, stream >>>(DD,
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
extern "C" void getSendGsDevF3(
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
	int Grid = (buffmax / numberOfThreads) + 1;
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
	dim3 threads(numberOfThreads, 1, 1);

	getSendGsF3 <<< grid, threads >>> (
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
extern "C" void setRecvGsDevF3(
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
	int Grid = (buffmax / numberOfThreads) + 1;
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
	dim3 threads(numberOfThreads, 1, 1);

	setRecvGsF3 <<< grid, threads >>> (
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
extern "C" void WallFuncDev27(unsigned int numberOfThreads,
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
   int Grid = (numberOfBCnodes / numberOfThreads)+1;
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

      WallFunction27<<< gridQ, threads >>> (
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
extern "C" void SetOutputWallVelocitySP27(unsigned int numberOfThreads,
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
   int Grid = (numberOfWallNodes / numberOfThreads)+1;
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

      LBSetOutputWallVelocitySP27<<< gridQ, threads >>> (	vxD,
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
extern "C" void GetVelotoForce27(unsigned int numberOfThreads,
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
   int Grid = (nonAtBC / numberOfThreads)+1;
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

      GetVeloforForcing27<<< gridQ, threads >>> (DD,
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
extern "C" void InitParticlesDevice(real* coordX,
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
   int Grid = (numberOfParticles / numberOfThreads)+1;
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

   InitParticles<<< gridQ, threads >>> (coordX,
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
extern "C" void MoveParticlesDevice(real* coordX,
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
   int Grid = (numberOfParticles / numberOfThreads)+1;
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

   MoveParticles<<< gridQ, threads >>> (coordX,
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
extern "C" void initRandomDevice(curandState* state,
								 unsigned int size_Mat,
								 unsigned int numberOfThreads)
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
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

   initRandom<<< gridQ, threads >>> (state);
   getLastCudaError("initRandom execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void generateRandomValuesDevice( curandState* state,
											unsigned int size_Mat,
											real* randArray,
											unsigned int numberOfThreads)
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
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

   generateRandomValues<<< gridQ, threads >>> (state,randArray);
   getLastCudaError("generateRandomValues execution failed");
}
//////////////////////////////////////////////////////////////////////////
extern "C" void CalcTurbulenceIntensityDevice(
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
   dim3 gridQ(Grid1, Grid2);
   dim3 threads(numberOfThreads, 1, 1 );

   CalcTurbulenceIntensity<<<gridQ, threads>>>(
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













