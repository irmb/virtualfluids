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

void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint sizeOfTiArray)
{
    int lev = para->getFine();
	cudaManager->cudaAllocTurbulenceIntensity(lev, sizeOfTiArray);
    resetVelocityFluctuationsAndMeans(para, cudaManager, sizeOfTiArray);
    para->getParH(lev)->turbulenceIntensity.resize(sizeOfTiArray);
}


void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint sizeOfTiArray)
{
    int lev = para->getFine();
    cudaManager->cudaCopyTurbulenceIntensityDH(lev, sizeOfTiArray);

	for (uint i = 0; i < sizeOfTiArray; i++)
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


void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff, uint sizeOfTiArray) {
    
    calcVelocityAndFluctuations(para, cudaManager, tdiff, sizeOfTiArray);
    int lev = para->getFine();

    real fluc_squared;
    real v_mean_squared;
    for (uint i = 0; i < sizeOfTiArray; i++) {
        fluc_squared =
            (real)(1.0 / 3.0 * (para->getParH(lev)->vxx[i] + para->getParH(lev)->vyy[i] + para->getParH(lev)->vzz[i]));
        v_mean_squared = para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i] +
                         para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i] +
                         para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
        para->getParH(lev)->turbulenceIntensity[i] = (real)sqrt(fluc_squared / v_mean_squared);
    }
}


void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaManager, uint sizeOfTiArray)
{
    int lev = para->getFine();

    for (unsigned int i = 0; i < sizeOfTiArray; i++) {
        para->getParH(lev)->vxx[i]     = (real)0.0;
        para->getParH(lev)->vyy[i]     = (real)0.0;
        para->getParH(lev)->vzz[i]     = (real)0.0;
        para->getParH(lev)->vx_mean[i] = (real)0.0;
        para->getParH(lev)->vy_mean[i] = (real)0.0;
        para->getParH(lev)->vz_mean[i] = (real)0.0;
    }

    cudaManager->cudaCopyTurbulenceIntensityHD(lev, sizeOfTiArray);
}

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray)
{
    int lev = para->getFine();
    std::vector<real *> data = { para->getParH(lev)->turbulenceIntensity.data() };
    std::vector<std::string> datanames = { "ti" };
    writeTiStuffToFile(para, timestep, vectorOfSparseIndices, sizeOfTiArray, data, datanames);
}

void writeVeloFluctuationToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray) 
{
    int lev                            = para->getFine();
    std::vector<real *> data           = { para->getParH(lev)->vxx, para->getParH(lev)->vyy, para->getParH(lev)->vzz };
    std::vector<std::string> datanames = { "vxx", "vyy", "vzz" };
    writeTiStuffToFile(para, timestep, vectorOfSparseIndices, sizeOfTiArray, data, datanames);
}

void writeVeloMeansToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray) {
    int lev                            = para->getFine();
    std::vector<real *> data           = { para->getParH(lev)->vx_mean, para->getParH(lev)->vy_mean, para->getParH(lev)->vz_mean };
    std::vector<std::string> datanames = { "vx_mean", "vy_mean", "vz_mean" };
    writeTiStuffToFile(para, timestep, vectorOfSparseIndices, sizeOfTiArray, data, datanames);
}

void writeAllTiDatafToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray)
{
    int lev                            = para->getFine();
    std::vector<real *> data           = { para->getParH(lev)->vxx,
                                           para->getParH(lev)->vyy,
                                           para->getParH(lev)->vzz,
                                           para->getParH(lev)->vx_mean,
                                           para->getParH(lev)->vy_mean,
                                           para->getParH(lev)->vz_mean,
                                           para->getParH(lev)->turbulenceIntensity.data() };
    std::vector<std::string> datanames = { "vxx", "vyy", "vzz", "vx_mean", "vy_mean", "vz_mean", "ti" };
    writeTiStuffToFile(para, timestep, vectorOfSparseIndices, sizeOfTiArray, data, datanames);
}

void writeTiStuffToFile(Parameter *para, uint timestep, int *vectorOfSparseIndices, int sizeOfTiArray, std::vector<real *> &data,
                        std::vector<std::string> &datanames)
{
    //////////////////////////////////////////////////////////////////////////
    // set filename
    std::string names;
    std::for_each(datanames.begin(), datanames.end(), [&names](const std::string &s) { return names += "_" + s; });
    std::string ffname = para->getFName() + StringUtil::toString<int>(para->getMyID()) + "_" +
                         StringUtil::toString<int>(timestep) + names + "_ti.txt";
    const char *fname = ffname.c_str();
    ////////////////////////////////////////////////////////////////////////
    // set ofstream
    std::ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    // open file
    ostr.open(fname);
    //////////////////////////////////////////////////////////////////////////
    // add header
    ostr << "index_sp";
        for (auto name : datanames) ostr << "\t" << name;
    ostr << std::endl;
    ////////////////////////////////////////////////////////////////////////
    // fill file with data
    for (int i = 0; i < sizeOfTiArray; i++) {
        ostr << vectorOfSparseIndices[i];
        for (auto dataset : data)
            ostr << "\t" << dataset[i];
        ostr << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////
    // close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
}