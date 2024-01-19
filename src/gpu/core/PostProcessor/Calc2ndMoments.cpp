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
//======================================================================================
#include "Calc2ndMoments.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Calc2ndMoments27.cuh"

#include "Cuda/CudaMemoryManager.h"
#include "Parameter/Parameter.h"

void alloc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //allocation (device-memory + host-memory)
        cudaMemoryManager->cudaAlloc2ndMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}



void init2ndMoments(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //init host arrays
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->kxyFromfcNEQ[pos]    = 0.0;
            para->getParH(lev)->kyzFromfcNEQ[pos]    = 0.0;
            para->getParH(lev)->kxzFromfcNEQ[pos]    = 0.0;
            para->getParH(lev)->kxxMyyFromfcNEQ[pos] = 0.0;
            para->getParH(lev)->kxxMzzFromfcNEQ[pos] = 0.0;
        }
        //////////////////////////////////////////////////////////////////////////
    }
}



void calc2ndMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //calc 2nd Moments on device
        Calc2ndMomentsCompSP27( para->getParD(lev)->kxyFromfcNEQ,
                                para->getParD(lev)->kyzFromfcNEQ,
                                para->getParD(lev)->kxzFromfcNEQ,
                                para->getParD(lev)->kxxMyyFromfcNEQ,
                                para->getParD(lev)->kxxMzzFromfcNEQ,
                                para->getParD(lev)->typeOfGridNode,       
                                para->getParD(lev)->neighborX, 
                                para->getParD(lev)->neighborY, 
                                para->getParD(lev)->neighborZ,
                                para->getParD(lev)->numberOfNodes, 
                                para->getParD(lev)->numberofthreads, 
                                para->getParD(lev)->distributions.f[0],    
                                para->getParD(lev)->isEvenTimestep);
        //////////////////////////////////////////////////////////////////////////
        //copy results to host
        cudaMemoryManager->cudaCopy2ndMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}





































void alloc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //allocation (device-memory + host-memory)
        cudaMemoryManager->cudaAlloc3rdMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}



void init3rdMoments(Parameter* para)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //init host arrays
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->CUMbbb[pos] = 0.0;
            para->getParH(lev)->CUMabc[pos] = 0.0;
            para->getParH(lev)->CUMbac[pos] = 0.0;
            para->getParH(lev)->CUMbca[pos] = 0.0;
            para->getParH(lev)->CUMcba[pos] = 0.0;
            para->getParH(lev)->CUMacb[pos] = 0.0;
            para->getParH(lev)->CUMcab[pos] = 0.0;
        }
        //////////////////////////////////////////////////////////////////////////
    }
}



void calc3rdMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //calc 2nd Moments on device
        Calc3rdMomentsCompSP27( para->getParD(lev)->CUMbbb,
                                para->getParD(lev)->CUMabc,
                                para->getParD(lev)->CUMbac,
                                para->getParD(lev)->CUMbca,
                                para->getParD(lev)->CUMcba,
                                para->getParD(lev)->CUMacb,
                                para->getParD(lev)->CUMcab,
                                para->getParD(lev)->typeOfGridNode,       
                                para->getParD(lev)->neighborX, 
                                para->getParD(lev)->neighborY, 
                                para->getParD(lev)->neighborZ,
                                para->getParD(lev)->numberOfNodes, 
                                para->getParD(lev)->numberofthreads, 
                                para->getParD(lev)->distributions.f[0],    
                                para->getParD(lev)->isEvenTimestep);
        //////////////////////////////////////////////////////////////////////////
        //copy results to host
        cudaMemoryManager->cudaCopy3rdMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}





































void allocHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //allocation (device-memory + host-memory)
        cudaMemoryManager->cudaAllocHigherMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}



void initHigherOrderMoments(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //init host arrays
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++)
        {
            para->getParH(lev)->CUMcbb[pos] = 0.0;
            para->getParH(lev)->CUMbcb[pos] = 0.0;
            para->getParH(lev)->CUMbbc[pos] = 0.0;
            para->getParH(lev)->CUMcca[pos] = 0.0;
            para->getParH(lev)->CUMcac[pos] = 0.0;
            para->getParH(lev)->CUMacc[pos] = 0.0;
            para->getParH(lev)->CUMbcc[pos] = 0.0;
            para->getParH(lev)->CUMcbc[pos] = 0.0;
            para->getParH(lev)->CUMccb[pos] = 0.0;
            para->getParH(lev)->CUMccc[pos] = 0.0;
        }
        //////////////////////////////////////////////////////////////////////////
    }
}



void calcHigherOrderMoments(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
        //////////////////////////////////////////////////////////////////////////
        //calc 2nd Moments on device
        CalcHigherMomentsCompSP27(  para->getParD(lev)->CUMcbb,
                                    para->getParD(lev)->CUMbcb,
                                    para->getParD(lev)->CUMbbc,
                                    para->getParD(lev)->CUMcca,
                                    para->getParD(lev)->CUMcac,
                                    para->getParD(lev)->CUMacc,
                                    para->getParD(lev)->CUMbcc,
                                    para->getParD(lev)->CUMcbc,
                                    para->getParD(lev)->CUMccb,
                                    para->getParD(lev)->CUMccc,
                                    para->getParD(lev)->typeOfGridNode,       
                                    para->getParD(lev)->neighborX, 
                                    para->getParD(lev)->neighborY, 
                                    para->getParD(lev)->neighborZ,
                                    para->getParD(lev)->numberOfNodes, 
                                    para->getParD(lev)->numberofthreads, 
                                    para->getParD(lev)->distributions.f[0],    
                                    para->getParD(lev)->isEvenTimestep);
        //////////////////////////////////////////////////////////////////////////
        //copy results to host
        cudaMemoryManager->cudaCopyHigherMoments(lev, para->getParH(lev)->numberOfNodes);
        //////////////////////////////////////////////////////////////////////////
    }
}

//! \}
