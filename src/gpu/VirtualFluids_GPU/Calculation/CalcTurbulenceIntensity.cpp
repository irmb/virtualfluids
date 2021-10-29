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
#include <basics/Core/StringUtilities/StringUtil.h>

void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint size)
{
    int lev = para->getFine();
	cudaManager->cudaAllocTurbulenceIntensity(lev, size);
    resetVelocityFluctuationsAndMeans(para, cudaManager, size);
    para->getParH(lev)->turbulenceIntensity.resize(size);
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


void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaManager, uint size)
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

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices)
{
    //////////////////////////////////////////////////////////////////////////
    // set level
    int lev = para->getFine();
    //////////////////////////////////////////////////////////////////////////
    // set filename
    std::string ffname = para->getFName() + StringUtil::toString<int>(para->getMyID()) + "_" +
                         StringUtil::toString<int>(timestep) + "_TurbulenIntensity.txt";
    const char *fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    // set ofstream
    std::ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    // open file
    ostr.open(fname);
    //////////////////////////////////////////////////////////////////////////
    // add header
    ostr << "index_sp" << "\t" << "ti" << std::endl;
    //////////////////////////////////////////////////////////////////////////
    // fill file with data
    for (size_t i = 0; i < para->getParH(lev)->turbulenceIntensity.size(); i++) {
        ostr << vectorOfSparseIndices[i] << "\t" << para->getParH(lev)->turbulenceIntensity[i]
             << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////
    // close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
}
