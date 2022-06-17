//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "Calculation/CalcMedian.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

void allocMedian(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		cudaMemoryManager->cudaAllocMedianOut(lev);
		for (unsigned int i = 0; i < para->getParH(lev)->numberOfNodes; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vy_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vz_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->rho_SP_Med_Out[i]   = (real)0.0;
			para->getParH(lev)->press_SP_Med_Out[i] = (real)0.0;
		}
	}
}





void calcMedian(Parameter* para, uint tdiff)
{
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		for (uint i = 0; i < para->getParH(lev)->numberOfNodes; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = para->getParH(lev)->vx_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->vy_SP_Med_Out[i]    = para->getParH(lev)->vy_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->vz_SP_Med_Out[i]    = para->getParH(lev)->vz_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->rho_SP_Med_Out[i]   = para->getParH(lev)->rho_SP_Med[i]  / (real)tdiff;
			para->getParH(lev)->press_SP_Med_Out[i] = para->getParH(lev)->press_SP_Med[i]/ (real)tdiff;
		}
	}
}





void resetMedian(Parameter* para)
{
	for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
		ResetMedianValuesSP27(
			para->getParD(lev)->vx_SP_Med,
			para->getParD(lev)->vy_SP_Med,
			para->getParD(lev)->vz_SP_Med,
			para->getParD(lev)->rho_SP_Med,
			para->getParD(lev)->press_SP_Med,
			para->getParD(lev)->numberOfNodes,
			para->getParD(lev)->numberofthreads,
			para->getParD(lev)->isEvenTimestep);
		getLastCudaError("ResetMedianValuesSP27 execution failed");
	}
}





//Advection-Diffusion
void allocMedianAD(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
	for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
        cudaMemoryManager->cudaAllocMedianOutAD(lev);
		for (unsigned int i = 0; i < para->getParH(lev)->numberOfNodes; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vy_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vz_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->rho_SP_Med_Out[i]   = (real)0.0;
			para->getParH(lev)->press_SP_Med_Out[i] = (real)0.0;
			para->getParH(lev)->Conc_Med_Out[i]     = (real)0.0;
		}
	}
}





void calcMedianAD(Parameter* para, uint tdiff)
{
	for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
		for (uint i = 0; i < para->getParH(lev)->numberOfNodes; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = para->getParH(lev)->vx_SP_Med[i]    / (real)tdiff;
			para->getParH(lev)->vy_SP_Med_Out[i]    = para->getParH(lev)->vy_SP_Med[i]    / (real)tdiff;
			para->getParH(lev)->vz_SP_Med_Out[i]    = para->getParH(lev)->vz_SP_Med[i]    / (real)tdiff;
			para->getParH(lev)->rho_SP_Med_Out[i]   = para->getParH(lev)->rho_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->press_SP_Med_Out[i] = para->getParH(lev)->press_SP_Med[i] / (real)tdiff;
			para->getParH(lev)->Conc_Med_Out[i]     = para->getParH(lev)->Conc_Med[i]     / (real)tdiff;
		}
	}
}





void resetMedianAD(Parameter* para)
{
	for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
		ResetMedianValuesAD27(
			para->getParD(lev)->vx_SP_Med,
			para->getParD(lev)->vy_SP_Med,
			para->getParD(lev)->vz_SP_Med,
			para->getParD(lev)->rho_SP_Med,
			para->getParD(lev)->press_SP_Med,
			para->getParD(lev)->Conc_Med,
			para->getParD(lev)->numberOfNodes,
			para->getParD(lev)->numberofthreads,
			para->getParD(lev)->isEvenTimestep);
		getLastCudaError("ResetMedianValuesSP27 execution failed");
	}
}



