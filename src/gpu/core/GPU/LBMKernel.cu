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
//! \author Martin Schoenherr
//=======================================================================================
// includes, cuda
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include <cuda_helper/CudaGrid.h>

// includes, kernels
#include "GPU/GPU_Kernels.cuh"

#include "Parameter/Parameter.h"
//////////////////////////////////////////////////////////////////////////
void CalcMac27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int grid_nx,
    unsigned int grid_ny,
    unsigned int grid_nz,
    real* DD,
    bool isEvenTimestep)
{
   dim3 threads       ( grid_nx, 1, 1 );
   dim3 grid          ( grid_ny, grid_nz );

    LBCalcMac27<<< grid, threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMac27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacSP27<<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMacSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacCompSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacCompSP27<<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMacCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacThS7(
    real* Conc,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD7,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    CalcConc7<<< grid.grid, grid.threads >>> (
        Conc,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD7,
        isEvenTimestep);
    getLastCudaError("CalcConc7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PlaneConcThS7(
    real* Conc,
    int* kPC,
    unsigned int numberOfPointskPC,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD7,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

    GetPlaneConc7<<< grid.grid, grid.threads >>> (
        Conc,
        kPC,
        numberOfPointskPC,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD7,
        isEvenTimestep);
    getLastCudaError("GetPlaneConc7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PlaneConcThS27(
    real* Conc,
    int* kPC,
    unsigned int numberOfPointskPC,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD27,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

    GetPlaneConc27<<< grid.grid, grid.threads >>> (
        Conc,
        kPC,
        numberOfPointskPC,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD27,
        isEvenTimestep);
    getLastCudaError("GetPlaneConc27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcConcentration27(
    unsigned int numberOfThreads,
    real* Conc,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* DD27,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    CalcConc27<<< grid.grid, grid.threads >>> (
        Conc,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD27,
        isEvenTimestep);
    getLastCudaError("CalcConc27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMedSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMedSP27<<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMedSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMedCompSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMedCompSP27<<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMedCompSP27 execution failed");
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
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    real* DD_AD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

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
        numberOfLBnodes,
        DD,
        DD_AD,
        isEvenTimestep);
    getLastCudaError("LBCalcMedCompAD27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcMacMedSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int tdiff,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacMedSP27<<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        tdiff,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBCalcMacMedSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMeanValuesSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesSP27 <<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBResetMeanValuesSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMeanValuesAD27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    real* concD,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesAD27 <<< grid.grid, grid.threads >>> (
        vxD,
        vyD,
        vzD,
        rhoD,
        pressD,
        concD,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBResetMeanValuesAD27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc2ndMomentsIncompSP27(
    real* kxyFromfcNEQ,
    real* kyzFromfcNEQ,
    real* kxzFromfcNEQ,
    real* kxxMyyFromfcNEQ,
    real* kxxMzzFromfcNEQ,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc2ndMomentsIncompSP27<<< grid.grid, grid.threads >>> (
        kxyFromfcNEQ,
        kyzFromfcNEQ,
        kxzFromfcNEQ,
        kxxMyyFromfcNEQ,
        kxxMzzFromfcNEQ,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalc2ndMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc2ndMomentsCompSP27(
    real* kxyFromfcNEQ,
    real* kyzFromfcNEQ,
    real* kxzFromfcNEQ,
    real* kxxMyyFromfcNEQ,
    real* kxxMzzFromfcNEQ,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc2ndMomentsCompSP27<<< grid.grid, grid.threads >>> (
        kxyFromfcNEQ,
        kyzFromfcNEQ,
        kxzFromfcNEQ,
        kxxMyyFromfcNEQ,
        kxxMzzFromfcNEQ,
        geoD,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalc2ndMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc3rdMomentsIncompSP27(
    real* CUMbbb,
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
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc3rdMomentsIncompSP27<<< grid.grid, grid.threads >>> (
        CUMbbb,
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
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBCalc3rdMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void Calc3rdMomentsCompSP27(
    real* CUMbbb,
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
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc3rdMomentsCompSP27<<< grid.grid, grid.threads >>> (
        CUMbbb,
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
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBCalc3rdMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcHigherMomentsIncompSP27(
    real* CUMcbb,
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
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcHigherMomentsIncompSP27<<< grid.grid, grid.threads >>> (
        CUMcbb,
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
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBCalcHigherMomentsIncompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcHigherMomentsCompSP27(
    real* CUMcbb,
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
    unsigned long long numberOfLBnodes,
    unsigned int numberOfThreads,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcHigherMomentsCompSP27<<< grid.grid, grid.threads >>> (
        CUMcbb,
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
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("LBCalcHigherMomentsCompSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void LBCalcMeasurePoints27(
    real* vxMP,
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
    unsigned long long numberOfLBnodes,
    real* DD,
    unsigned int numberOfThreads,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskMP);

    LBCalcMeasurePoints<<< grid.grid, grid.threads >>> (
        vxMP,
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
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBCalcMeasurePoints execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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

    scaleCF27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
void ScaleCFEff27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFEff27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCFEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFLast27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFLast27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCFLast27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFpress27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFpress27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCFpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_Fix_27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_Fix_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_Fix_comp_27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_Fix_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_0817_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_0817_comp_27<<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_0817_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3_2018(
    real* DC,
    real* DF,
    real* G6,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_comp_D3Q27F3_2018 <<< grid.grid, grid.threads >>>(
        DC,
        DF,
        G6,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_comp_D3Q27F3(
    real* DC,
    real* DF,
    real* G6,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_comp_D3Q27F3 <<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        G6,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_comp_D3Q27F3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_staggered_time_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_staggered_time_comp_27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_staggered_time_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_comp_27(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st *stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads,  coarseToFine->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1 );

    scaleCF_RhoSq_comp_27<<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        coarseToFine->coarseCellIndices,
        coarseToFine->fineCellIndices,
        coarseToFine->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->viscosity,
        parameterDeviceC->nx,
        parameterDeviceC->ny,
        parameterDeviceF->nx,
        parameterDeviceF->ny,
        neighborCoarseToFine);
    getLastCudaError("scaleCF_RhoSq_27 execution failed");
}

template<bool hasTurbulentViscosity> void ScaleCF_compressible(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st *stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads,  coarseToFine->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1 );

    scaleCF_compressible<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        coarseToFine->coarseCellIndices,
        coarseToFine->fineCellIndices,
        coarseToFine->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity,
        neighborCoarseToFine);

    getLastCudaError("scaleCF_compressible execution failed");
}
template void ScaleCF_compressible<true>(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st *stream);
template void ScaleCF_compressible<false>(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * coarseToFine, ICellNeigh& neighborCoarseToFine, CUstream_st *stream);

//////////////////////////////////////////////////////////////////////////
void ScaleCF_RhoSq_3rdMom_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_RhoSq_3rdMom_comp_27<<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_AA2016_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_AA2016_comp_27<<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCF_NSPress_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCF_NSPress_27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborCoarseToFine);
    getLastCudaError("scaleCF_NSPress_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThSMG7(
    real* DC,
    real* DF,
    real* DD7C,
    real* DD7F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posCSWB,
    unsigned int* posFSWB,
    unsigned int kCF,
    real nu,
    real diffusivity_fine,
    unsigned int numberOfThreads,
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThSMG7<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        DD7C,
        DD7F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posCSWB,
        posFSWB,
        kCF,
        nu,
        diffusivity_fine,
        neighborCoarseToFine);
    getLastCudaError("scaleCFThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS7(
    real* DC,
    real* DF,
    real* DD7C,
    real* DD7F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posCSWB,
    unsigned int* posFSWB,
    unsigned int kCF,
    real nu,
    real diffusivity_fine,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThS7<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        DD7C,
        DD7F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posCSWB,
        posFSWB,
        kCF,
        nu,
        diffusivity_fine);
    getLastCudaError("scaleCFThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleCFThS27(
    real* DC,
    real* DF,
    real* DD27C,
    real* DD27F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posCSWB,
    unsigned int* posFSWB,
    unsigned int kCF,
    real nu,
    real diffusivity_fine,
    unsigned int numberOfThreads,
    ICellNeigh neighborCoarseToFine)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kCF);

    scaleCFThS27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        DD27C,
        DD27F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posCSWB,
        posFSWB,
        kCF,
        nu,
        diffusivity_fine,
        neighborCoarseToFine);
    getLastCudaError("scaleCFThS27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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

    scaleFC27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
void ScaleFCEff27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCEff27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFCEff27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCLast27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCLast27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("Kernel execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCpress27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCpress27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFCpress27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_Fix_27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_Fix_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_Fix_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_Fix_comp_27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_Fix_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_0817_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_0817_comp_27<<< grid.grid, grid.threads, 0, stream >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_0817_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3_2018(
    real* DC,
    real* DF,
    real* G6,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_comp_D3Q27F3_2018 <<< grid.grid, grid.threads >>> (
        DC,
        DF,
        G6,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_comp_D3Q27F3_2018 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_comp_D3Q27F3(
    real* DC,
    real* DF,
    real* G6,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_comp_D3Q27F3 <<< grid.grid, grid.threads, 0, stream >>> (
        DC,
        DF,
        G6,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_comp_D3Q27F3 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_staggered_time_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_staggered_time_comp_27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_staggered_time_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_comp_27(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * fineToCoarse, ICellNeigh &neighborFineToCoarse, CUstream_st *stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads,  fineToCoarse->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1 );

    scaleFC_RhoSq_comp_27<<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        fineToCoarse->coarseCellIndices,
        fineToCoarse->fineCellIndices,
        fineToCoarse->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->viscosity,
        parameterDeviceC->nx,
        parameterDeviceC->ny,
        parameterDeviceF->nx,
        parameterDeviceF->ny,
        neighborFineToCoarse);
    getLastCudaError("scaleFC_RhoSq_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
template<bool hasTurbulentViscosity> void ScaleFC_compressible(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * fineToCoarse, ICellNeigh &neighborFineToCoarse, CUstream_st *stream)
{
    dim3 grid = vf::cuda::getCudaGrid(parameterDeviceC->numberofthreads,  fineToCoarse->numberOfCells);
    dim3 threads(parameterDeviceC->numberofthreads, 1, 1 );

    scaleFC_compressible<hasTurbulentViscosity><<<grid, threads, 0, stream>>>(
        parameterDeviceC->distributions.f[0],
        parameterDeviceF->distributions.f[0],
        parameterDeviceC->neighborX,
        parameterDeviceC->neighborY,
        parameterDeviceC->neighborZ,
        parameterDeviceF->neighborX,
        parameterDeviceF->neighborY,
        parameterDeviceF->neighborZ,
        parameterDeviceC->numberOfNodes,
        parameterDeviceF->numberOfNodes,
        parameterDeviceC->isEvenTimestep,
        fineToCoarse->coarseCellIndices,
        fineToCoarse->fineCellIndices,
        fineToCoarse->numberOfCells,
        parameterDeviceC->omega,
        parameterDeviceF->omega,
        parameterDeviceC->turbViscosity,
        parameterDeviceF->turbViscosity,
        neighborFineToCoarse);

    getLastCudaError("scaleFC_compressible execution failed");
}
template void ScaleFC_compressible<true>(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * fineToCoarse, ICellNeigh &neighborFineToCoarse, CUstream_st *stream);
template void ScaleFC_compressible<false>(LBMSimulationParameter * parameterDeviceC, LBMSimulationParameter* parameterDeviceF, ICells * fineToCoarse, ICellNeigh &neighborFineToCoarse, CUstream_st *stream);

//////////////////////////////////////////////////////////////////////////
void ScaleFC_RhoSq_3rdMom_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_RhoSq_3rdMom_comp_27<<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_RhoSq_3rdMom_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_AA2016_comp_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse,
    CUstream_st *stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_AA2016_comp_27<<< grid.grid, grid.threads, 0, stream >>>(
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_AA2016_comp_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFC_NSPress_27(
    real* DC,
    real* DF,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
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
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFC_NSPress_27<<< grid.grid, grid.threads >>> (
        DC,
        DF,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
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
        neighborFineToCoarse);
    getLastCudaError("scaleFC_NSPress_27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThSMG7(
    real* DC,
    real* DF,
    real* DD7C,
    real* DD7F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posC,
    unsigned int* posFSWB,
    unsigned int kFC,
    real nu,
    real diffusivity_coarse,
    unsigned int numberOfThreads,
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThSMG7<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        DD7C,
        DD7F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posC,
        posFSWB,
        kFC,
        nu,
        diffusivity_coarse,
        neighborFineToCoarse);
    getLastCudaError("scaleFCThSMG7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS7(
    real* DC,
    real* DF,
    real* DD7C,
    real* DD7F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posC,
    unsigned int* posFSWB,
    unsigned int kFC,
    real nu,
    real diffusivity_coarse,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThS7<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        DD7C,
        DD7F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posC,
        posFSWB,
        kFC,
        nu,
        diffusivity_coarse);
    getLastCudaError("scaleFCThS7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ScaleFCThS27(
    real* DC,
    real* DF,
    real* DD27C,
    real* DD27F,
    unsigned int* neighborCX,
    unsigned int* neighborCY,
    unsigned int* neighborCZ,
    unsigned int* neighborFX,
    unsigned int* neighborFY,
    unsigned int* neighborFZ,
    unsigned long long numberOfLBnodesC,
    unsigned long long numberOfLBnodesF,
    bool isEvenTimestep,
    unsigned int* posC,
    unsigned int* posFSWB,
    unsigned int kFC,
    real nu,
    real diffusivity_coarse,
    unsigned int numberOfThreads,
    ICellNeigh neighborFineToCoarse)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, kFC);

    scaleFCThS27<<< grid.grid, grid.threads >>>(
        DC,
        DF,
        DD27C,
        DD27F,
        neighborCX,
        neighborCY,
        neighborCZ,
        neighborFX,
        neighborFY,
        neighborFZ,
        numberOfLBnodesC,
        numberOfLBnodesF,
        isEvenTimestep,
        posC,
        posFSWB,
        kFC,
        nu,
        diffusivity_coarse,
        neighborFineToCoarse);
    getLastCudaError("scaleFCThS27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void DragLiftPostD27(
    real* DD,
    int* k_Q,
    real* QQ,
    int numberOfBCnodes,
    double *DragX,
    double *DragY,
    double *DragZ,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    DragLiftPost27<<< grid.grid, grid.threads >>>(
        DD,
        k_Q,
        QQ,
        numberOfBCnodes,
        DragX,
        DragY,
        DragZ,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("DragLiftPost27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void DragLiftPreD27(
    real* DD,
    int* k_Q,
    real* QQ,
    int numberOfBCnodes,
    double *DragX,
    double *DragY,
    double *DragZ,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    DragLiftPre27<<< grid.grid, grid.threads >>>(
        DD,
        k_Q,
        QQ,
        numberOfBCnodes,
        DragX,
        DragY,
        DragZ,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("DragLiftPre27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcCPtop27(
    real* DD,
    int* cpIndex,
    int nonCp,
    double *cpPress,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

    CalcCP27<<< grid.grid, grid.threads >>>(
        DD,
        cpIndex,
        nonCp,
        cpPress,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("CalcCP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void CalcCPbottom27(
    real* DD,
    int* cpIndex,
    int nonCp,
    double *cpPress,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

    CalcCP27<<< grid.grid, grid.threads >>>(
        DD,
        cpIndex,
        nonCp,
        cpPress,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("CalcCP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
void SetOutputWallVelocitySP27(
    unsigned int numberOfThreads,
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
    unsigned long long numberOfLBnodes,
    real* DD,
    bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfWallNodes);

    LBSetOutputWallVelocitySP27<<< grid.grid, grid.threads >>> (
        vxD,
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
        numberOfLBnodes,
        DD,
        isEvenTimestep);
    getLastCudaError("LBSetOutputWallVelocitySP27 execution failed");
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
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    uint numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);
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
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("CalcTurbulenceIntensity execution failed");
}
