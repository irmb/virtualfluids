//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "CalcMean.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "PostProcessor/MacroscopicQuantities.cuh"

void allocMean(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        cudaMemoryManager->cudaAllocMeanOut(lev);
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->vx_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->vy_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->vz_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->rho_SP_Med_Out[pos]   = (real)0.0;
            para->getParH(lev)->press_SP_Med_Out[pos] = (real)0.0;
        }
    }
}





void calcMean(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->vx_SP_Med_Out[pos]    = para->getParH(lev)->vx_SP_Med[pos]   / (real)tdiff;
            para->getParH(lev)->vy_SP_Med_Out[pos]    = para->getParH(lev)->vy_SP_Med[pos]   / (real)tdiff;
            para->getParH(lev)->vz_SP_Med_Out[pos]    = para->getParH(lev)->vz_SP_Med[pos]   / (real)tdiff;
            para->getParH(lev)->rho_SP_Med_Out[pos]   = para->getParH(lev)->rho_SP_Med[pos]  / (real)tdiff;
            para->getParH(lev)->press_SP_Med_Out[pos] = para->getParH(lev)->press_SP_Med[pos]/ (real)tdiff;
        }
    }
}





void resetMean(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        ResetMeanValuesSP27(
            para->getParD(lev)->vx_SP_Med,
            para->getParD(lev)->vy_SP_Med,
            para->getParD(lev)->vz_SP_Med,
            para->getParD(lev)->rho_SP_Med,
            para->getParD(lev)->press_SP_Med,
            para->getParD(lev)->numberOfNodes,
            para->getParD(lev)->numberofthreads,
            para->getParD(lev)->isEvenTimestep);
        getLastCudaError("ResetMeanValuesSP27 execution failed");
    }
}





//Advection-Diffusion
void allocMeanAD(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        cudaMemoryManager->cudaAllocMeanOutAD(lev);
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->vx_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->vy_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->vz_SP_Med_Out[pos]    = (real)0.0;
            para->getParH(lev)->rho_SP_Med_Out[pos]   = (real)0.0;
            para->getParH(lev)->press_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->Conc_Med_Out[pos]     = (real)0.0;
        }
    }
}





void calcMeanAD(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->vx_SP_Med_Out[pos]    = para->getParH(lev)->vx_SP_Med[pos]    / (real)tdiff;
            para->getParH(lev)->vy_SP_Med_Out[pos]    = para->getParH(lev)->vy_SP_Med[pos]    / (real)tdiff;
            para->getParH(lev)->vz_SP_Med_Out[pos]    = para->getParH(lev)->vz_SP_Med[pos]    / (real)tdiff;
            para->getParH(lev)->rho_SP_Med_Out[pos]   = para->getParH(lev)->rho_SP_Med[pos]   / (real)tdiff;
            para->getParH(lev)->press_SP_Med_Out[pos] = para->getParH(lev)->press_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->Conc_Med_Out[pos]     = para->getParH(lev)->Conc_Med[pos]     / (real)tdiff;
        }
    }
}





void resetMeanAD(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        ResetMeanValuesAD27(
            para->getParD(lev)->vx_SP_Med,
            para->getParD(lev)->vy_SP_Med,
            para->getParD(lev)->vz_SP_Med,
            para->getParD(lev)->rho_SP_Med,
            para->getParD(lev)->press_SP_Med,
            para->getParD(lev)->Conc_Med,
            para->getParD(lev)->numberOfNodes,
            para->getParD(lev)->numberofthreads,
            para->getParD(lev)->isEvenTimestep);
        getLastCudaError("ResetMeanValuesAD27 execution failed");
    }
}



