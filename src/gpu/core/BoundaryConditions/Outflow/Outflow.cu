//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ /
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "LBM/LB.h"
#include <cuda_helper/CudaGrid.h>

#include "BoundaryConditions/Outflow/Outflow_Device.cuh"
#include "Parameter/Parameter.h"

void OutflowNonReflecting(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    OutflowNonReflecting_Device<<< grid, threads >>> (
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
        vf::lbm::dir::dP00);
    getLastCudaError("OutflowNonReflecting_Device execution failed");
}

void OutflowNonReflectingPressureCorrection(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition)
{
    dim3 grid = vf::cuda::getCudaGrid( parameterDevice->numberofthreads,  boundaryCondition->numberOfBCnodes);
    dim3 threads(parameterDevice->numberofthreads, 1, 1 );

    OutflowNonReflectingPressureCorrection_Device<<< grid, threads >>> (
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
        vf::lbm::dir::dP00,
        parameterDevice->outflowPressureCorrectionFactor);
    getLastCudaError("OutflowNonReflectingPressureCorrection_Device execution failed");
}
