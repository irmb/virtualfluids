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
#include "ForceCalculations.h"

//////////////////////////////////////////////////////////////////////////

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "Cuda/CudaMemoryManager.h"

#include "PostProcessor/MacroscopicQuantities.cuh"

#include "StringUtilities/StringUtil.h"
//using namespace std;
//////////////////////////////////////////////////////////////////////////

ForceCalculations::ForceCalculations(Parameter* para)
{
    Ta = 10.0;        // number of time steps between adjusting
    vx1Targed = para->getVelocity(); // objected LB velocity

    Kpcrit = 3.0 / Ta;// 0.3;
    Tcrit = 3.0 * Ta; // 30.0;
    Tn = 1.0 * Tcrit; //0.5 * Tcrit;
    Tv = 0.24 * Tcrit; //0.12 * Tcrit;

    Kp = 0.6 * Kpcrit;
    Ki = Kp / Tn;
    Kd = Kp * Tv;

    y = 0.0;
    e = 0.0;
    esum = 0.0;
    eold = 0.0;

    isPID = true;
}

void ForceCalculations::calcPIDControllerForForce(Parameter* para, CudaMemoryManager* cudaMemoryManager)
 {
     //////////////////////////////////////////////////////////////////////////
     double tempVeloX = 0.0;
     double veloAverageX = 0.0;
     double levelVeloAverageX = 0.0;
     int counter = 0;
     //////////////////////////////////////////////////////////////////////////
     for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
     {
         //////////////////////////////////////////////////////////////////////
         //measure the velocity
         unsigned long long numberOfElements = para->getParH(lev)->numberOfNodes;
         if (numberOfElements > 0)
         {
             calculateMacroscopicQuantitiesCompressible(para->getParD(lev)->velocityX,
                             para->getParD(lev)->velocityY,
                             para->getParD(lev)->velocityZ,
                             para->getParD(lev)->rho,
                             para->getParD(lev)->pressure,
                             para->getParD(lev)->typeOfGridNode,
                             para->getParD(lev)->neighborX,
                             para->getParD(lev)->neighborY,
                             para->getParD(lev)->neighborZ,
                             para->getParD(lev)->numberOfNodes,
                             para->getParD(lev)->numberofthreads,
                             para->getParD(lev)->distributions.f[0],
                             para->getParD(lev)->isEvenTimestep);
             getLastCudaError("calculateMacroscopicQuantities execution failed");
             //////////////////////////////////////////////////////////////////
             cudaMemoryManager->cudaCopyPrint(lev);
             //////////////////////////////////////////////////////////////////
             for (size_t pos = 0; pos < numberOfElements; pos++)
             {
                 tempVeloX += (double)para->getParH(lev)->velocityX[pos];
             }
             tempVeloX /= (double)numberOfElements;
             //////////////////////////////////////////////////////////////////
             levelVeloAverageX += tempVeloX;
             //////////////////////////////////////////////////////////////////
             counter++;
             //////////////////////////////////////////////////////////////////
         }
     }
     //////////////////////////////////////////////////////////////////////////
     veloAverageX = levelVeloAverageX / (double)counter;
     //////////////////////////////////////////////////////////////////////////
     if (isPID)
     {
         //PID-Controller
         e = vx1Targed - veloAverageX;
         esum = esum + e;
         y = Kp * e + Ki * Ta * esum + Kd * (e - eold) / Ta;
         eold = e;

         y = y / 2.0;
     }
     //////////////////////////////////////////////////////////////////////////
     para->getForcesDouble()[0] = (para->getForcesDouble()[0] + y);
     para->getForcesDouble()[1] = (para->getForcesDouble()[1] + y) * 0.0;
     para->getForcesDouble()[2] = (para->getForcesDouble()[2] + y) * 0.0;
     //////////////////////////////////////////////////////////////////////////
     para->getForcesHost()[0] = (real)(para->getForcesHost()[0] + y);
     para->getForcesHost()[1] = (real)(para->getForcesHost()[1] + y) * (real)0.0;
     para->getForcesHost()[2] = (real)(para->getForcesHost()[2] + y) * (real)0.0;
     //////////////////////////////////////////////////////////////////////////
     cudaMemoryManager->cudaCopyForcingToDevice();
     //////////////////////////////////////////////////////////////////////////
 }


void ForceCalculations::printForcing(Parameter* para)
{
    //////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffname = para->getFName() + StringUtil::toString<int>(para->getMyProcessID()) + "_forcing.txt";
    const char* fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostr.open(fname, std::fstream::app);
    ostr << para->getForcesHost()[0] << " " << para->getForcesHost()[1] << " " << para->getForcesHost()[2];
    ostr << std::endl;
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
}

//! \}
