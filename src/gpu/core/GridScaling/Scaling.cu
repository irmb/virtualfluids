//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================
#include "Scaling.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "LBM/LB.h"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void scaleCF27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine,
                          real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF);

__global__ void scaleCFEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                             unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                             unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                             unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu, unsigned int nxC,
                             unsigned int nyC, unsigned int nxF, unsigned int nyF, ICellNeigh neighborCoarseToFine);

__global__ void scaleCFLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                              unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                              unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                              unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                              unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu, unsigned int nxC,
                              unsigned int nyC, unsigned int nxF, unsigned int nyF, ICellNeigh neighborCoarseToFine);

__global__ void scaleCFpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                               unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                               unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                    unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                     unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                     unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                     unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                     unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                     unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                     ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                          unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                          unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                          unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                          ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                                     unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                     unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                     unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                     unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                     unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                     ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                               unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                                               unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse,
                                               real omFine, real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF,
                                               unsigned int nyF, ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_RhoSq_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                      unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                      unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                      unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                      unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                      unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                      ICellNeigh neighborCoarseToFine);

template <bool hasTurbulentViscosity>
__global__ void scaleCF_compressible(real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse,
                                     unsigned int* neighborYcoarse, unsigned int* neighborZcoarse,
                                     unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine,
                                     unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine,
                                     bool isEvenTimestep, unsigned int* indicesCoarseMMM, unsigned int* indicesFineMMM,
                                     unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine,
                                     real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetCF);

__global__ void scaleCF_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                             unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                             unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                             unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                             unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                             ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                       unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                       unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                       unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                       unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                       unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                       ICellNeigh neighborCoarseToFine);

__global__ void scaleCF_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                   unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                   unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                   unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                   ICellNeigh neighborCoarseToFine);

__global__ void scaleCFThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                              unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                              unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                              unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                              unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                              ICellNeigh neighborCoarseToFine);

__global__ void scaleCFThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                            unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                            unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                            unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                            unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine);

__global__ void scaleCFThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX,
                             unsigned int* neighborCY, unsigned int* neighborCZ, unsigned int* neighborFX,
                             unsigned int* neighborFY, unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                             unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                             ICellNeigh neighborCoarseToFine);

// fine to coarse
__global__ void scaleFC27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                          unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF);

__global__ void scaleFCEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                             unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                             unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                             unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu, unsigned int nxC,
                             unsigned int nyC, unsigned int nxF, unsigned int nyF, ICellNeigh neighborFineToCoarse);

__global__ void scaleFCLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                              unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                              unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                              unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                              unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu, unsigned int nxC,
                              unsigned int nyC, unsigned int nxF, unsigned int nyF, ICellNeigh neighborFineToCoarse);

__global__ void scaleFCpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                               unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                               unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                    unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                     unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                     unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                     unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                     unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                     unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                     ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                          unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                          unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                          unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                          ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                                     unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                     unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                     unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                     unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                     unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                     ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                               unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                               ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_RhoSq_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                      unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                      unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                      unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                      unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                      unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                      ICellNeigh neighborFineToCoarse);

template <bool hasTurbulentViscosity>
__global__ void scaleFC_compressible(real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse,
                                     unsigned int* neighborYcoarse, unsigned int* neighborZcoarse,
                                     unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine,
                                     unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine,
                                     bool isEvenTimestep, unsigned int* indicesCoarse000, unsigned int* indicesFineMMM,
                                     unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine,
                                     real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetFC);

__global__ void scaleFC_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                             unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                             unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                             unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                             unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                             ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                       unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                       unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                       unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                       unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                       unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                       ICellNeigh neighborFineToCoarse);

__global__ void scaleFC_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                   unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                   unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                   unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                   ICellNeigh neighborFineToCoarse);

__global__ void scaleFCThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                              unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                              unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                              unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                              unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                              ICellNeigh neighborFineToCoarse);

__global__ void scaleFCThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                            unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                            unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                            unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                            unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse);

__global__ void scaleFCThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX,
                             unsigned int* neighborCY, unsigned int* neighborCZ, unsigned int* neighborFX,
                             unsigned int* neighborFY, unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                             unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                             unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                             ICellNeigh neighborFineToCoarse);

void ScaleCF27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
               unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
               unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
               unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ,
                                           numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posCSWB, posFSWB, kCF,
                                           omCoarse, omFine, nu, nxC, nyC, nxF, nyF);
    getLastCudaError("scaleCF27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                  unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                  ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFEff27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ,
                                              numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posCSWB, posFSWB, kCF,
                                              omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCFEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                   unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                   ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFLast27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY,
                                               neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posCSWB,
                                               posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCFLast27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFpress27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCFpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_Fix_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                         unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                         unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                         unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                         unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                         unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_Fix_comp_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_Fix_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine,
                          real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_0817_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_0817_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                               unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_comp_D3Q27F3_2018<<<grid.grid, grid.threads>>>(DC, DF, G6, neighborCX, neighborCY, neighborCZ, neighborFX,
                                                           neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
                                                           isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC,
                                                           nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF,
                          bool isEvenTimestep, unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse,
                          real omFine, real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_comp_D3Q27F3<<<grid.grid, grid.threads, 0, stream>>>(DC, DF, G6, neighborCX, neighborCY, neighborCZ, neighborFX,
                                                                 neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
                                                                 isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu,
                                                                 nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_comp_D3Q27F3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                    unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_staggered_time_comp_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_staggered_time_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_comp_27(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                           ICells* coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, coarseToFine->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleCF_RhoSq_comp_27<<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0], parameterDeviceF->distributions.f[0], parameterDeviceC->neighborX,
        parameterDeviceC->neighborY, parameterDeviceC->neighborZ, parameterDeviceF->neighborX, parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ, parameterDeviceC->numberOfNodes, parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep, coarseToFine->coarseCellIndices, coarseToFine->fineCellIndices,
        coarseToFine->numberOfCells, parameterDeviceC->omega, parameterDeviceF->omega, parameterDeviceC->viscosity,
        parameterDeviceC->nx, parameterDeviceC->ny, parameterDeviceF->nx, parameterDeviceF->ny, neighborCoarseToFine);
    getLastCudaError("scaleCF_RhoSq_27 execution failed");
}

template <bool hasTurbulentViscosity>
void ScaleCF_compressible(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                          ICells* coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, coarseToFine->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleCF_compressible<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0], parameterDeviceF->distributions.f[0], parameterDeviceC->neighborX,
        parameterDeviceC->neighborY, parameterDeviceC->neighborZ, parameterDeviceF->neighborX, parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ, parameterDeviceC->numberOfNodes, parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep, coarseToFine->coarseCellIndices, coarseToFine->fineCellIndices,
        coarseToFine->numberOfCells, parameterDeviceC->omega, parameterDeviceF->omega, parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity, neighborCoarseToFine);

    getLastCudaError("scaleCF_compressible execution failed");
}
template void ScaleCF_compressible<true>(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                                         ICells* coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream);
template void ScaleCF_compressible<false>(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                                          ICells* coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream);

//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                  unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                  unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                  unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                  unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_RhoSq_3rdMom_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                            unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                            unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                            unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine,
                            real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                            unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_AA2016_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                        unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                        unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                        unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                        unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                        ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_NSPress_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posCSWB, posFSWB, kCF, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborCoarseToFine);
    getLastCudaError("scaleCF_NSPress_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                   unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThSMG7<<<grid.grid, grid.threads>>>(DC, DF, DD7C, DD7F, neighborCX, neighborCY, neighborCZ, neighborFX,
                                               neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep,
                                               posCSWB, posFSWB, kCF, nu, diffusivity_fine, neighborCoarseToFine);
    getLastCudaError("scaleCFThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                 unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                 unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                 unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                 unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThS7<<<grid.grid, grid.threads>>>(DC, DF, DD7C, DD7F, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY,
                                             neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posCSWB,
                                             posFSWB, kCF, nu, diffusivity_fine);
    getLastCudaError("scaleCFThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX, unsigned int* neighborCY,
                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                  unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThS27<<<grid.grid, grid.threads>>>(DC, DF, DD27C, DD27F, neighborCX, neighborCY, neighborCZ, neighborFX,
                                              neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep,
                                              posCSWB, posFSWB, kCF, nu, diffusivity_fine, neighborCoarseToFine);
    getLastCudaError("scaleCFThS27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
               unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
               unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
               unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ,
                                           numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posC, posFSWB, kFC, omCoarse,
                                           omFine, nu, nxC, nyC, nxF, nyF);
    getLastCudaError("scaleFC27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                  unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                  ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCEff27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ,
                                              numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posC, posFSWB, kFC,
                                              omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFCEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                   unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                   ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCLast27<<<grid.grid, grid.threads>>>(DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY,
                                               neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posC, posFSWB,
                                               kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("Kernel execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCpress27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFCpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_Fix_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                         unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                         unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                         unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                         unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                         unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_Fix_comp_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_Fix_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                          unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_0817_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_0817_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                               unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_comp_D3Q27F3_2018<<<grid.grid, grid.threads>>>(DC, DF, G6, neighborCX, neighborCY, neighborCZ, neighborFX,
                                                           neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
                                                           isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC,
                                                           nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF,
                          bool isEvenTimestep, unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse,
                          real omFine, real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_comp_D3Q27F3<<<grid.grid, grid.threads, 0, stream>>>(DC, DF, G6, neighborCX, neighborCY, neighborCZ, neighborFX,
                                                                 neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
                                                                 isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu,
                                                                 nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_comp_D3Q27F3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                    unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_staggered_time_comp_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_staggered_time_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_comp_27(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                           ICells* fineToCoarse, ICellNeigh& neighborFineToCoarse, CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, fineToCoarse->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleFC_RhoSq_comp_27<<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0], parameterDeviceF->distributions.f[0], parameterDeviceC->neighborX,
        parameterDeviceC->neighborY, parameterDeviceC->neighborZ, parameterDeviceF->neighborX, parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ, parameterDeviceC->numberOfNodes, parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep, fineToCoarse->coarseCellIndices, fineToCoarse->fineCellIndices,
        fineToCoarse->numberOfCells, parameterDeviceC->omega, parameterDeviceF->omega, parameterDeviceC->viscosity,
        parameterDeviceC->nx, parameterDeviceC->ny, parameterDeviceF->nx, parameterDeviceF->ny, neighborFineToCoarse);
    getLastCudaError("scaleFC_RhoSq_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
template <bool hasTurbulentViscosity>
void ScaleFC_compressible(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                          ICells* fineToCoarse, ICellNeigh& neighborFineToCoarse, CUstream_st* stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads, fineToCoarse->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1);

    scaleFC_compressible<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0], parameterDeviceF->distributions.f[0], parameterDeviceC->neighborX,
        parameterDeviceC->neighborY, parameterDeviceC->neighborZ, parameterDeviceF->neighborX, parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ, parameterDeviceC->numberOfNodes, parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep, fineToCoarse->coarseCellIndices, fineToCoarse->fineCellIndices,
        fineToCoarse->numberOfCells, parameterDeviceC->omega, parameterDeviceF->omega, parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity, neighborFineToCoarse);

    getLastCudaError("scaleFC_compressible execution failed");
}
template void ScaleFC_compressible<true>(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                                         ICells* fineToCoarse, ICellNeigh& neighborFineToCoarse, CUstream_st* stream);
template void ScaleFC_compressible<false>(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                                          ICells* fineToCoarse, ICellNeigh& neighborFineToCoarse, CUstream_st* stream);

//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                  unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                  unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                  unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                  unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_RhoSq_3rdMom_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                            unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                            unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                            unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                            unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                            unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_AA2016_comp_27<<<grid.grid, grid.threads, 0, stream>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                        unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                        unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                        unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                        unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                        ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_NSPress_27<<<grid.grid, grid.threads>>>(
        DC, DF, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF,
        isEvenTimestep, posC, posFSWB, kFC, omCoarse, omFine, nu, nxC, nyC, nxF, nyF, neighborFineToCoarse);
    getLastCudaError("scaleFC_NSPress_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                   unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThSMG7<<<grid.grid, grid.threads>>>(DC, DF, DD7C, DD7F, neighborCX, neighborCY, neighborCZ, neighborFX,
                                               neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep,
                                               posC, posFSWB, kFC, nu, diffusivity_coarse, neighborFineToCoarse);
    getLastCudaError("scaleFCThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                 unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                 unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                 unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                 unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThS7<<<grid.grid, grid.threads>>>(DC, DF, DD7C, DD7F, neighborCX, neighborCY, neighborCZ, neighborFX, neighborFY,
                                             neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep, posC, posFSWB,
                                             kFC, nu, diffusivity_coarse);
    getLastCudaError("scaleFCThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX, unsigned int* neighborCY,
                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                  unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThS27<<<grid.grid, grid.threads>>>(DC, DF, DD27C, DD27F, neighborCX, neighborCY, neighborCZ, neighborFX,
                                              neighborFY, neighborFZ, numberOfLBnodesC, numberOfLBnodesF, isEvenTimestep,
                                              posC, posFSWB, kFC, nu, diffusivity_coarse, neighborFineToCoarse);
    getLastCudaError("scaleFCThS27 execution failed");
}
