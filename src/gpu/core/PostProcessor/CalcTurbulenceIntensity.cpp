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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//=======================================================================================
#include "CalcTurbulenceIntensity.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <basics/StringUtilities/StringUtil.h>

#include "Cuda/CudaMemoryManager.h"
#include "Parameter/Parameter.h"

void allocTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaMemoryManager->cudaAllocTurbulenceIntensity(lev, para->getParH(lev)->numberOfNodes);
        para->getParH(lev)->turbulenceIntensity.resize(para->getParH(lev)->numberOfNodes);    
    }
        resetVelocityFluctuationsAndMeans(para, cudaMemoryManager);
}


void calcVelocityAndFluctuations(Parameter *para, CudaMemoryManager *cudaMemoryManager, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaMemoryManager->cudaCopyTurbulenceIntensityDH(lev, para->getParH(lev)->numberOfNodes);

        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            // mean velocity
            para->getParH(lev)->vx_mean[pos] = para->getParH(lev)->vx_mean[pos] / (real)tdiff;
            para->getParH(lev)->vy_mean[pos] = para->getParH(lev)->vy_mean[pos] / (real)tdiff;
            para->getParH(lev)->vz_mean[pos] = para->getParH(lev)->vz_mean[pos] / (real)tdiff;

            // fluctuations
            para->getParH(lev)->vxx[pos] = para->getParH(lev)->vxx[pos] / (real)tdiff;
            para->getParH(lev)->vyy[pos] = para->getParH(lev)->vyy[pos] / (real)tdiff;
            para->getParH(lev)->vzz[pos] = para->getParH(lev)->vzz[pos] / (real)tdiff;
            para->getParH(lev)->vxy[pos] = para->getParH(lev)->vxy[pos] / (real)tdiff;
            para->getParH(lev)->vxz[pos] = para->getParH(lev)->vxz[pos] / (real)tdiff;
            para->getParH(lev)->vyz[pos] = para->getParH(lev)->vyz[pos] / (real)tdiff;

            para->getParH(lev)->vxx[pos] =
                para->getParH(lev)->vxx[pos] - para->getParH(lev)->vx_mean[pos] * para->getParH(lev)->vx_mean[pos];
            para->getParH(lev)->vyy[pos] =
                para->getParH(lev)->vyy[pos] - para->getParH(lev)->vy_mean[pos] * para->getParH(lev)->vy_mean[pos];
            para->getParH(lev)->vzz[pos] =
                para->getParH(lev)->vzz[pos] - para->getParH(lev)->vz_mean[pos] * para->getParH(lev)->vz_mean[pos];
            para->getParH(lev)->vxy[pos] =
                para->getParH(lev)->vxy[pos] - para->getParH(lev)->vx_mean[pos] * para->getParH(lev)->vy_mean[pos];
            para->getParH(lev)->vxz[pos] =
                para->getParH(lev)->vxz[pos] - para->getParH(lev)->vx_mean[pos] * para->getParH(lev)->vz_mean[pos];
            para->getParH(lev)->vyz[pos] =
                para->getParH(lev)->vyz[pos] - para->getParH(lev)->vy_mean[pos] * para->getParH(lev)->vz_mean[pos];
        }
    }
}


void calcTurbulenceIntensity(Parameter *para, CudaMemoryManager *cudaMemoryManager, uint tdiff) {
    

    real fluc_squared;
    real v_mean_squared;

    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
    calcVelocityAndFluctuations(para, cudaMemoryManager, tdiff);

        for (uint i = 0; i < para->getParH(lev)->numberOfNodes; i++) {
            fluc_squared = (real)(
                1.0 / 3.0 * (para->getParH(lev)->vxx[i] + para->getParH(lev)->vyy[i] + para->getParH(lev)->vzz[i]));
            v_mean_squared = para->getParH(lev)->vx_mean[i] * para->getParH(lev)->vx_mean[i] +
                             para->getParH(lev)->vy_mean[i] * para->getParH(lev)->vy_mean[i] +
                             para->getParH(lev)->vz_mean[i] * para->getParH(lev)->vz_mean[i];
            para->getParH(lev)->turbulenceIntensity[i] = (real)sqrt(fluc_squared / v_mean_squared);
        }
    }
}


void resetVelocityFluctuationsAndMeans(Parameter *para, CudaMemoryManager *cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (unsigned int i = 0; i < para->getParH(lev)->numberOfNodes; i++) {
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

        cudaMemoryManager->cudaCopyTurbulenceIntensityHD(lev, para->getParH(lev)->numberOfNodes);
    }
}

void cudaFreeTurbulenceIntensityArrays(Parameter *para, CudaMemoryManager *cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaMemoryManager->cudaFreeTurbulenceIntensity(lev);
    }
}

void writeTurbulenceIntensityToFile(Parameter *para, uint timestep)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data           = { para->getParH(lev)->turbulenceIntensity.data() };
        std::vector<std::string> datanames = { "ti" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->numberOfNodes, data, datanames);
    }
}

void writeVeloFluctuationToFile(Parameter *para, uint timestep) 
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data = { para->getParH(lev)->vxx, para->getParH(lev)->vyy, para->getParH(lev)->vzz };
        std::vector<std::string> datanames = { "vxx", "vyy", "vzz" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->numberOfNodes, data, datanames);
    }
}

void writeVeloMeansToFile(Parameter *para, uint timestep) {
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        std::vector<real *> data           = { para->getParH(lev)->vx_mean, 
                                               para->getParH(lev)->vy_mean,
                                               para->getParH(lev)->vz_mean };
        std::vector<std::string> datanames = { "vx_mean", "vy_mean", "vz_mean" };
        writeTiStuffToFile(para, timestep, para->getParH(lev)->numberOfNodes, data, datanames);
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
        writeTiStuffToFile(para, timestep, para->getParH(lev)->numberOfNodes, data, datanames);
    }
}

void writeTiStuffToFile(Parameter *para, uint timestep, unsigned long long sizeOfTiArray, std::vector<real *> &data,
                        std::vector<std::string> &datanames)
{
    ////////////////////////////////////////////////////////////////////////
    // set filename
    std::string names;
    std::for_each(datanames.begin(), datanames.end(), [&names](const std::string &s) { return names += "_" + s; });
    std::string ffname = para->getFName() + StringUtil::toString<int>(para->getMyProcessID()) + "_" +
                         StringUtil::toString<int>(timestep) + names + "_ti.txt";
    const char *fname = ffname.c_str();
    ////////////////////////////////////////////////////////////////////////
    // set ofstream
    std::ofstream ostr;
    ////////////////////////////////////////////////////////////////////////
    // open file
    ostr.open(fname);
    ////////////////////////////////////////////////////////////////////////
    // add header
    ostr << "index_sp";
        for (auto name : datanames) ostr << "\t" << name;
    ostr << std::endl;
    ////////////////////////////////////////////////////////////////////////
    // fill file with data
    for (size_t pos = 0; pos < sizeOfTiArray; pos++) {
        ostr << pos;
        for (auto dataset : data)
            ostr << "\t" << dataset[pos];
        ostr << std::endl;
    }
    ////////////////////////////////////////////////////////////////////////
    // close file
    ostr.close();
    ////////////////////////////////////////////////////////////////////////
}

//! \}
