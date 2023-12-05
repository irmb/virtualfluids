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

#ifndef SCALING_CUH
#define SCALING_CUH

#include <basics/DataTypes.h>

#include "Calculation/Calculation.h"

struct LBMSimulationParameter;

void ScaleCF27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
               unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
               unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
               unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads);

void ScaleFC27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
               unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
               unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
               unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads);

void ScaleCFEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                  unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                  ICellNeigh neighborCoarseToFine);

void ScaleFCEff27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                  unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                  ICellNeigh neighborFineToCoarse);

void ScaleCFLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                   unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                   ICellNeigh neighborCoarseToFine);

void ScaleFCLast27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                   unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                   unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                   ICellNeigh neighborFineToCoarse);

void ScaleCFpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborCoarseToFine);

void ScaleFCpress27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborFineToCoarse);

void ScaleCF_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborCoarseToFine);

void ScaleCF_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                         unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                         unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                         unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                         unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                         unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine);

void ScaleCF_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine,
                          real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream);

void ScaleCF_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                               unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine);

void ScaleCF_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF,
                          bool isEvenTimestep, unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse,
                          real omFine, real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream);

void ScaleCF_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                    unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine);

void ScaleCF_RhoSq_comp_27(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                           ICells* interpolationCellsCoarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream);

template <bool hasTurbulentViscosity>
void ScaleCF_compressible(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                          ICells* interpolationCellsCoarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st* stream);

void ScaleCF_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                  unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                  unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posCSWB,
                                  unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                  unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream);

void ScaleCF_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                            unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                            unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                            unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine,
                            real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                            unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine, CUstream_st* stream);

void ScaleCF_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                        unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                        unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                        unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real omCoarse, real omFine, real nu,
                        unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                        ICellNeigh neighborCoarseToFine);

void ScaleFC_Fix_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                    unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                    unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                    unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                    ICellNeigh neighborFineToCoarse);

void ScaleFC_Fix_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                         unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                         unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                         unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                         unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                         unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse);

void ScaleFC_0817_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                          unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                          unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                          unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                          unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream);

void ScaleFC_comp_D3Q27F3_2018(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                               unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                               unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                               unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                               unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                               unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                               unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse);

void ScaleFC_comp_D3Q27F3(real* DC, real* DF, real* G6, unsigned int* neighborCX, unsigned int* neighborCY,
                          unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                          unsigned int* neighborFZ, unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF,
                          bool isEvenTimestep, unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse,
                          real omFine, real nu, unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                          unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream);

void ScaleFC_staggered_time_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                    unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                    unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                    unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                    unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                    unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                    unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse);

void ScaleFC_RhoSq_comp_27(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                           ICells* interpolationCellsFineToCoarse, ICellNeigh& neighborFineToCoarse, CUstream_st* stream);

template <bool hasTurbulentViscosity>
void ScaleFC_compressible(LBMSimulationParameter* parameterDeviceC, LBMSimulationParameter* parameterDeviceF,
                          ICells* icellFC, ICellNeigh& neighborFineToCoarse, CUstream_st* stream);

void ScaleFC_RhoSq_3rdMom_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY,
                                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY,
                                  unsigned int* neighborFZ, unsigned long long numberOfLBnodesC,
                                  unsigned long long numberOfLBnodesF, bool isEvenTimestep, unsigned int* posC,
                                  unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                                  unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                                  unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream);

void ScaleFC_AA2016_comp_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                            unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                            unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                            unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                            unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF,
                            unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse, CUstream_st* stream);

void ScaleFC_NSPress_27(real* DC, real* DF, unsigned int* neighborCX, unsigned int* neighborCY, unsigned int* neighborCZ,
                        unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                        unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                        unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real omCoarse, real omFine, real nu,
                        unsigned int nxC, unsigned int nyC, unsigned int nxF, unsigned int nyF, unsigned int numberOfThreads,
                        ICellNeigh neighborFineToCoarse);

void ScaleCFThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                 unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                 unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                 unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                 unsigned int numberOfThreads);

void ScaleFCThS7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                 unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                 unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                 unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                 unsigned int numberOfThreads);

void ScaleCFThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                   unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine);

void ScaleFCThSMG7(real* DC, real* DF, real* DD7C, real* DD7F, unsigned int* neighborCX, unsigned int* neighborCY,
                   unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                   unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                   unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                   unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse);

void ScaleCFThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX, unsigned int* neighborCY,
                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posCSWB, unsigned int* posFSWB, unsigned int kCF, real nu, real diffusivity_fine,
                  unsigned int numberOfThreads, ICellNeigh neighborCoarseToFine);

void ScaleFCThS27(real* DC, real* DF, real* DD27C, real* DD27F, unsigned int* neighborCX, unsigned int* neighborCY,
                  unsigned int* neighborCZ, unsigned int* neighborFX, unsigned int* neighborFY, unsigned int* neighborFZ,
                  unsigned long long numberOfLBnodesC, unsigned long long numberOfLBnodesF, bool isEvenTimestep,
                  unsigned int* posC, unsigned int* posFSWB, unsigned int kFC, real nu, real diffusivity_coarse,
                  unsigned int numberOfThreads, ICellNeigh neighborFineToCoarse);

#endif
