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

void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaManager->cudaAllocTurbulenceIntensity(lev, para->getParH(lev)->size_Mat_SP);
        para->getParH(lev)->turbulenceIntensity.resize(para->getParH(lev)->size_Mat_SP);    
    }
        resetVelocityFluctuationsAndMeans(para, cudaManager);
}


void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaManager->cudaCopyTurbulenceIntensityDH(lev, para->getParH(lev)->size_Mat_SP);

        for (uint i = 0; i < para->getParH(lev)->size_Mat_SP; i++) {
            // mean velocity
            para->getParH(lev)->vx_mean[i] = para->getParH(lev)->vx_mean[i] / (real)tdiff;
            para->getParH(lev)->vy_mean[i] = para->getParH(lev)->vy_mean[i] / (real)tdiff;
            para->getParH(lev)->vz_mean[i] = para->getParH(lev)->vz_mean[i] / (real)tdiff;

            // fluctuations
            para->getParH(lev)->vxx[i] = para->getParH(lev)->vxx[i] / (real)tdiff;
            para->getParH(lev)->vyy[i] = para->getParH(lev)->vyy[i] / (real)tdiff;
            para->getParH(lev)->vzz[i] = para->getParH(lev)->vzz[i] / (real)tdiff;
            para->getParH(lev)->vxy[i] = para->getParH(lev)->vxy[i] / (real)tdiff;
            para->getParH(lev)->vxz[i] = para->getParH(lev)->vxz[i] / (real)tdiff;
            para->getParH(lev)->vyz[i] = para->getParH(lev)->vyz[i] / (real)tdiff;

            para->getParH(lev)->vxx[i] =
                para->getParH(lev)->vxx[i] - para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i];
            para->getParH(lev)->vyy[i] =
                para->getParH(lev)->vyy[i] - para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i];
            para->getParH(lev)->vzz[i] =
                para->getParH(lev)->vzz[i] - para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
            para->getParH(lev)->vxy[i] =
                para->getParH(lev)->vxy[i] - para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vy_mean[i];
            para->getParH(lev)->vxz[i] =
                para->getParH(lev)->vxz[i] - para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vz_mean[i];
            para->getParH(lev)->vyz[i] =
                para->getParH(lev)->vyz[i] - para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vz_mean[i];
        }
    }
}


void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaManager, uint tdiff) {
    

    real fluc_squared;
    real v_mean_squared;

    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
    calcVelocityAndFluctuations(para, cudaManager, tdiff);

        for (uint i = 0; i < para->getParH(lev)->size_Mat_SP; i++) {
            fluc_squared = (real)(
                1.0 / 3.0 * (para->getParH(lev)->vxx[i] + para->getParH(lev)->vyy[i] + para->getParH(lev)->vzz[i]));
            v_mean_squared = para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i] +
                             para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i] +
                             para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
            para->getParH(lev)->turbulenceIntensity[i] = (real)sqrt(fluc_squared / v_mean_squared);
        }
    }
}


void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (unsigned int i = 0; i < para->getParH(lev)->size_Mat_SP; i++) {
            para->getParH(lev)->vxx[i]     = (real)0.0;
            para->getParH(lev)->vyy[i]     = (real)0.0;
            para->getParH(lev)->vzz[i]     = (real)0.0;
            para->getParH(lev)->vxy[i]     = (real)0.0;
            para->getParH(lev)->vxz[i]     = (real)0.0;
            para->getParH(lev)->vyz[i]     = (real)0.0;
            para->getParH(lev)->vx_mean[i] = (real)0.0;
            para->getParH(lev)->vy_mean[i] = (real)0.0;
            para->getParH(lev)->vz_mean[i] = (real)0.0;
        }

        cudaManager->cudaCopyTurbulenceIntensityHD(lev, para->getParH(lev)->size_Mat_SP);
    }
}

void cudaFreeTurbulenceIntensityArrays(Parameter *para, CudaMemoryManager *cudaManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaManager->cudaFreeTurbulenceIntensity(lev);
    }
}

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data           = { para->getParH(lev)->turbulenceIntensity.data() };
        std::vector<std::string> datanames = { "ti" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->size_Mat_SP, data, datanames);
    }
}

void writeVeloFluctuationToFile(Parameter *para, uint timestep) 
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data = { para->getParH(lev)->vxx, para->getParH(lev)->vyy, para->getParH(lev)->vzz };
        std::vector<std::string> datanames = { "vxx", "vyy", "vzz" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->size_Mat_SP, data, datanames);
    }
}

void writeVeloMeansToFile(Parameter *para, uint timestep) {
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data           = { para->getParH(lev)->vx_mean, 
                                               para->getParH(lev)->vy_mean,
                                               para->getParH(lev)->vz_mean };
        std::vector<std::string> datanames = { "vx_mean", "vy_mean", "vz_mean" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->size_Mat_SP, data, datanames);
    }
}

void writeAllTiDatafToFile(Parameter *para, uint timestep)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data = { para->getParH(lev)->vxx,
                                     para->getParH(lev)->vyy,
                                     para->getParH(lev)->vzz,
                                     para->getParH(lev)->vx_mean,
                                     para->getParH(lev)->vy_mean,
                                     para->getParH(lev)->vz_mean,
                                     para->getParH(lev)->turbulenceIntensity.data() };
        std::vector<std::string> datanames = { "vxx", "vyy", "vzz", "vx_mean", "vy_mean", "vz_mean", "ti" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->size_Mat_SP, data, datanames);
    }
}

void writeTiStuffToFile(Parameter *para, uint timestep, int sizeOfTiArray, std::vector<real *> &data,
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
        ostr << i;
        for (auto dataset : data)
            ostr << "\t" << dataset[i];
        ostr << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////
    // close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
}