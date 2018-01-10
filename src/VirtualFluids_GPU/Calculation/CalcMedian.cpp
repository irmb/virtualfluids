#include "Calculation/CalcMedian.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

void allocMedian(Parameter* para)
{
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		para->cudaAllocMedianOut(lev);
		for (unsigned int i = 0; i < para->getParH(lev)->size_Mat_SP; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vy_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->vz_SP_Med_Out[i]    = (real)0.0;
			para->getParH(lev)->rho_SP_Med_Out[i]   = (real)0.0;
			para->getParH(lev)->press_SP_Med_Out[i] = (real)0.0;
		}
	}
}





void calcMedian(Parameter* para, unsigned int tdiff)
{
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		for (unsigned int i = 0; i < para->getParH(lev)->size_Mat_SP; i++)
		{
			para->getParH(lev)->vx_SP_Med_Out[i]    = para->getParH(lev)->vx_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->vy_SP_Med_Out[i]    = para->getParH(lev)->vy_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->vz_SP_Med_Out[i]    = para->getParH(lev)->vz_SP_Med[i]   / (real)tdiff;
			para->getParH(lev)->rho_SP_Med_Out[i]   = para->getParH(lev)->rho_SP_Med[i]  / (real)tdiff;
			para->getParH(lev)->press_SP_Med_Out[i] = para->getParH(lev)->press_SP_Med[i]/ (real)tdiff;
		}
	}
}
