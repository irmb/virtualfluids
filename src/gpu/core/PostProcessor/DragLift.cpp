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
//! \author Martin Schoenherr
//=======================================================================================

#include "DragLift.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cstdio>
#include <fstream>
#include <sstream>

#include "Cuda/CudaMemoryManager.h"
#include "Parameter/Parameter.h"
#include "StringUtilities/StringUtil.h"

using namespace std;

//! \brief Calculate drag and lift for a geometry
//! \details note, that the drag/lift calculations are build for being used with geometry boundary nodes and the related area has to be defined here
void calcDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev)
{
    //////////////////////////////////////////////////////////////////////////
    //copy to host
    //finest Grid ... with the geometry nodes
    cudaMemoryManager->cudaCopyDragLift(lev, para->getParH(lev)->geometryBC.numberOfBCnodes);
    //////////////////////////////////////////////////////////////////////////
    double dragX = 0., dragY = 0., dragZ = 0.;
    double CDX   = 0., CDY   = 0., CDZ   = 0.;
    //////////////////////////////////////////////////////////////////////////
    // the calculation of the related area (A) will change, please take care 
    double delta_x_F = 0.00625; //[m]
    double A  = 2.19/(delta_x_F*delta_x_F);
    //////////////////////////////////////////////////////////////////////////

    for (unsigned int it = 0; it < para->getParH(lev)->geometryBC.numberOfBCnodes; it++)
    {
        dragX += (double) (para->getParH(lev)->DragLiftPreProcessingInXdirection[it] - para->getParH(lev)->DragLiftPostProcessingInXdirection[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
        dragY += (double) (para->getParH(lev)->DragLiftPreProcessingInYdirection[it] - para->getParH(lev)->DragLiftPostProcessingInYdirection[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
        dragZ += (double) (para->getParH(lev)->DragLiftPreProcessingInZdirection[it] - para->getParH(lev)->DragLiftPostProcessingInZdirection[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
    }
    //////////////////////////////////////////////////////////////////////////
    //calc CD
    CDX = 2.0 * dragX / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    CDY = 2.0 * dragY / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    CDZ = 2.0 * dragZ / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    //////////////////////////////////////////////////////////////////////////
    //Copy to vector x,y,z
    para->getParH(lev)->DragLiftVectorInXdirection.push_back(CDX);
    para->getParH(lev)->DragLiftVectorInYdirection.push_back(CDY);
    para->getParH(lev)->DragLiftVectorInZdirection.push_back(CDZ);
    //////////////////////////////////////////////////////////////////////////
}



void allocDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    //////////////////////////////////////////////////////////////////////////
    //set level
    int lev = para->getMaxLevel();
    //////////////////////////////////////////////////////////////////////////
    //allocation
    //finest Grid ... with the geometry nodes
    cudaMemoryManager->cudaAllocDragLift(lev, para->getParH(lev)->geometryBC.numberOfBCnodes);
    //////////////////////////////////////////////////////////////////////////
    printf("\n number of elements for drag and lift = %d \n", para->getParH(lev)->geometryBC.numberOfBCnodes);
}



void printDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int timestep)
{
    //////////////////////////////////////////////////////////////////////////
    //set level
    int lev = para->getMaxLevel();
    //////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffname = para->getFName()+StringUtil::toString<int>(para->getMyProcessID())+"_"+StringUtil::toString<int>(timestep)+"_DragLift.txt";
    const char* fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostr.open(fname);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    for (size_t i = 0; i < para->getParH(lev)->DragLiftVectorInXdirection.size(); i++)
    {
        ostr << para->getParH(lev)->DragLiftVectorInXdirection[i] << "\t" << para->getParH(lev)->DragLiftVectorInYdirection[i] << "\t" << para->getParH(lev)->DragLiftVectorInZdirection[i] << endl ;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
    if (timestep == (int)para->getTimestepEnd())
    {
        cudaMemoryManager->cudaFreeDragLift(lev);
    }
    //////////////////////////////////////////////////////////////////////////
}
//! \}
