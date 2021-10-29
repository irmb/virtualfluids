//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "Calculation/CalcTurbulenceIntensity.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size)
{
    int lev = para->getFine();
	cudaManager->cudaAllocTurbulenceIntensity(lev, size);
    resetTurbulenceIntensity(para, cudaManager, size);
}


void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size)
{
    int lev = para->getFine();
    cudaManager->cudaCopyTurbulenceIntensityDH(lev, size);

	for (uint i = 0; i < size; i++)
	{
        // mean velocity
        para->getParH(lev)->vx_mean[i] = para->getParH(lev)->vx_mean[i]   / (real)tdiff;
        para->getParH(lev)->vy_mean[i] = para->getParH(lev)->vy_mean[i]   / (real)tdiff;
        para->getParH(lev)->vz_mean[i] = para->getParH(lev)->vz_mean[i]   / (real)tdiff;

        // fluctuations
        para->getParH(lev)->vxx[i]     = para->getParH(lev)->vxx[i]       / (real)tdiff;
        para->getParH(lev)->vyy[i]     = para->getParH(lev)->vyy[i]       / (real)tdiff;
        para->getParH(lev)->vzz[i]     = para->getParH(lev)->vzz[i]       / (real)tdiff;

        para->getParH(lev)->vxx[i] = para->getParH(lev)->vxx[i] - para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i];
        para->getParH(lev)->vyy[i] = para->getParH(lev)->vyy[i] - para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i];
        para->getParH(lev)->vzz[i] = para->getParH(lev)->vzz[i] - para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
	}
}


void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint size) {
    calcVelocityAndFluctuations(para, cudaManager, tdiff, size);
    int lev = para->getFine();

    para->getParH(lev)->turbulenceIntensity.resize(size);

    real fluc_squared;
    real v_mean_squared;
    for (uint i = 0; i < size; i++) {
        fluc_squared =
            (real)(1.0 / 3.0 * (para->getParH(lev)->vxx[i] + para->getParH(lev)->vyy[i] + para->getParH(lev)->vzz[i]));
        v_mean_squared = para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i] +
                         para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i] +
                         para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
        para->getParH(lev)->turbulenceIntensity[i] = (real)sqrt(fluc_squared / v_mean_squared);
    }
}


void resetTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size)
{
    int lev = para->getFine();

    for (unsigned int i = 0; i < size; i++) {
        para->getParH(lev)->vxx[i]     = (real)0.0;
        para->getParH(lev)->vyy[i]     = (real)0.0;
        para->getParH(lev)->vzz[i]     = (real)0.0;
        para->getParH(lev)->vx_mean[i] = (real)0.0;
        para->getParH(lev)->vy_mean[i] = (real)0.0;
        para->getParH(lev)->vz_mean[i] = (real)0.0;
    }

    cudaManager->cudaCopyTurbulenceIntensityHD(lev, size);
}

