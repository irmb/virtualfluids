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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LBKernelManager.cpp
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>
#include <stdexcept>
#include <string>

#include "BCKernelManager.h"
#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Calculation/Cp.h"
#include "Calculation/DragLift.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

BCKernelManager::BCKernelManager(SPtr<Parameter> parameter, BoundaryConditionFactory *bcFactory) : para(parameter)
{
    this->velocityBoundaryConditionPost = bcFactory->getVelocityBoundaryConditionPost();
    this->noSlipBoundaryConditionPost   = bcFactory->getNoSlipBoundaryConditionPost();
    this->slipBoundaryConditionPost     = bcFactory->getSlipBoundaryConditionPost();
    this->pressureBoundaryConditionPre  = bcFactory->getPressureBoundaryConditionPre();
    this->geometryBoundaryConditionPost = bcFactory->getGeometryBoundaryConditionPost();
    this->stressBoundaryConditionPost   = bcFactory->getStressBoundaryConditionPost();

    try {
        checkBoundaryCondition(this->velocityBoundaryConditionPost, this->para->getParD(0)->velocityBC,
                               "velocityBoundaryConditionPost");
        checkBoundaryCondition(this->noSlipBoundaryConditionPost, this->para->getParD(0)->noSlipBC,
                               "noSlipBoundaryConditionPost");
        checkBoundaryCondition(this->slipBoundaryConditionPost, this->para->getParD(0)->slipBC,
                               "slipBoundaryConditionPost");
        checkBoundaryCondition(this->pressureBoundaryConditionPre, this->para->getParD(0)->pressureBC,
                               "pressureBoundaryConditionPre");
        checkBoundaryCondition(this->geometryBoundaryConditionPost, this->para->getParD(0)->geometryBC,
                               "geometryBoundaryConditionPost");
        checkBoundaryCondition(this->stressBoundaryConditionPost, this->para->getParD(0)->stressBC,
                               "stressBoundaryConditionPost");
    } catch (const std::runtime_error &e) {
        std::cout << e.what() << std::endl;
        throw;
    } catch (const std::exception &e) {
        std::cout << "Unknown exception in BCKernelManager: " << e.what() << std::endl;
        throw;
    }
}

void BCKernelManager::checkBoundaryCondition(const boundaryCondition &bc, const QforBoundaryConditions &bcStruct,
                                             const std::string &bcName)
{
    if (!bc && bcStruct.numberOfBCnodes > 0)
        throw std::runtime_error("The boundary condition " + bcName + " was not set!");
}

void BCKernelManager::checkBoundaryCondition(const boundaryConditionPara &bc, const QforBoundaryConditions &bcStruct,
                                             const std::string &bcName)
{
    if (!bc && bcStruct.numberOfBCnodes > 0)
        throw std::runtime_error("The boundary condition " + bcName + " was not set!");
}

void BCKernelManager::runVelocityBCKernelPre(const int level) const
{
    if (para->getParD(level)->velocityBC.numberOfBCnodes > 0)
    {
        // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
        // if ( myid == 0)
        // {
        //    VelSchlaffer27(para->getParD(level)->numberofthreads, t,
        //                   para->getParD(level)->distributions.f[0],       para->getParD(level)->velocityBC.Vz,
        //                   para->getParD(level)->velocityBC.deltaVz, para->getParD(level)->velocityBC.k,
        //                   para->getParD(level)->velocityBC.kN,      para->getParD(level)->velocityBC.numberOfBCnodes,
        //                   para->getParD(level)->omega,           para->getParD(level)->neighborX,
        //                   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
        //                   para->getParD(level)->numberOfNodes,     para->getParD(level)->isEvenTimestep);
        //    getLastCudaError("VelSchlaffer27 execution failed");
        // }
        ////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        // QVelDevIncompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->velocityBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        ////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        // QVelDevCompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->velocityBC.Vx,
        //     para->getParD(level)->velocityBC.Vy,
        //     para->getParD(level)->velocityBC.Vz,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->velocityBC.k,
        //     para->getParD(level)->velocityBC.q27[0],
        //     para->getParD(level)->velocityBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void BCKernelManager::runVelocityBCKernelPost(const int level) const
{
     if (para->getParD(level)->velocityBC.numberOfBCnodes > 0)
     {
        velocityBoundaryConditionPost(para->getParD(level).get(), &(para->getParD(level)->velocityBC));

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////

        // QVelDevice1h27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                para->getParD(level)->velocityBC.Vx,      para->getParD(level)->velocityBC.Vy,   para->getParD(level)->velocityBC.Vz,
        //                para->getParD(level)->distributions.f[0],       para->getParD(level)->velocityBC.k,    para->getParD(level)->velocityBC.q27[0],
        //                para->getParD(level)->velocityBC.numberOfBCnodes,      para->getParD(level)->omega,
        //                para->getPhi(),                        para->getAngularVelocity(),
        //                para->getParD(level)->neighborX,    para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
        //                para->getParD(level)->coordinateX,       para->getParD(level)->coordinateY,    para->getParD(level)->coordinateZ,
        //                para->getParD(level)->numberOfNodes,     para->getParD(level)->isEvenTimestep);
        // getLastCudaError("QVelDev27 execution failed");
     }
}

void BCKernelManager::runGeoBCKernelPre(const int level, unsigned int t, CudaMemoryManager* cudaMemoryManager) const{
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0){
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part II
            DragLiftPreD27(
                para->getParD(level)->distributions.f[0],
                para->getParD(level)->geometryBC.k,
                para->getParD(level)->geometryBC.q27[0],
                para->getParD(level)->geometryBC.numberOfBCnodes,
                para->getParD(level)->DragPreX,
                para->getParD(level)->DragPreY,
                para->getParD(level)->DragPreZ,
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level)->isEvenTimestep,
                para->getParD(level)->numberofthreads);
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of Drag and Lift
            ////////////////////////////////////////////////////////////////////////////////
            calcDragLift(para.get(), cudaMemoryManager, level);
            ////////////////////////////////////////////////////////////////////////////////
        }

        if (para->getCalcCp())
        {
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of cp
            ////////////////////////////////////////////////////////////////////////////////

            if(t > para->getTStartOut())
            {
                ////////////////////////////////////////////////////////////////////////////////
                CalcCPtop27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpTopIndex,
                    para->getParD(level)->numberOfPointsCpTop,
                    para->getParD(level)->cpPressTop,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpBottomIndex,
                    para->getParD(level)->numberOfPointsCpBottom,
                    para->getParD(level)->cpPressBottom,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->cpBottom2Index,
                    para->getParD(level)->numberOfPointsCpBottom2,
                    para->getParD(level)->cpPressBottom2,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->isEvenTimestep,
                    para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                calcCp(para.get(), cudaMemoryManager, level);
            }            
        }

        ////////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        // QDevIncompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        // QDevCompHighNu27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->geometryBC.k,
        //     para->getParD(level)->geometryBC.q27[0],
        //     para->getParD(level)->geometryBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);

    }
}

void BCKernelManager::runGeoBCKernelPost(const int level) const
{
    if (para->getParD(level)->geometryBC.numberOfBCnodes > 0)
    {
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part I
            DragLiftPostD27(para->getParD(level)->distributions.f[0],
                            para->getParD(level)->geometryBC.k,
                            para->getParD(level)->geometryBC.q27[0],
                            para->getParD(level)->geometryBC.numberOfBCnodes,
                            para->getParD(level)->DragPostX,
                            para->getParD(level)->DragPostY,
                            para->getParD(level)->DragPostZ,
                            para->getParD(level)->neighborX,
                            para->getParD(level)->neighborY,
                            para->getParD(level)->neighborZ,
                            para->getParD(level)->numberOfNodes,
                            para->getParD(level)->isEvenTimestep,
                            para->getParD(level)->numberofthreads);
            getLastCudaError("DragLift27 execution failed");
        }

        geometryBoundaryConditionPost(para->getParD(level).get(), &(para->getParD(level)->geometryBC));

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////
        // the GridGenerator does currently not provide normals!

        //     QSlipGeomDevComp27(
        //         para->getParD(level)->numberofthreads,
        //         para->getParD(level)->distributions.f[0],
        //         para->getParD(level)->geometryBC.k,
        //         para->getParD(level)->geometryBC.q27[0],
        //         para->getParD(level)->geometryBC.numberOfBCnodes,
        //         para->getParD(level)->omega,
        //         para->getParD(level)->geometryBCnormalX.q27[0],
        //         para->getParD(level)->geometryBCnormalY.q27[0],
        //         para->getParD(level)->geometryBCnormalZ.q27[0],
        //         para->getParD(level)->neighborX,
        //         para->getParD(level)->neighborY,
        //         para->getParD(level)->neighborZ,
        //         para->getParD(level)->numberOfNodes,
        //         para->getParD(level)->isEvenTimestep);

        //     QSlipNormDevComp27(
        //         para->getParD(level)->numberofthreads,
        //         para->getParD(level)->distributions.f[0],
        //         para->getParD(level)->geometryBC.k,
        //         para->getParD(level)->geometryBC.q27[0],
        //         para->getParD(level)->geometryBC.numberOfBCnodes,
        //         para->getParD(level)->omega,
        //         para->getParD(level)->geometryBCnormalX.q27[0],
        //         para->getParD(level)->geometryBCnormalY.q27[0],
        //         para->getParD(level)->geometryBCnormalZ.q27[0],
        //         para->getParD(level)->neighborX,
        //         para->getParD(level)->neighborY,
        //         para->getParD(level)->neighborZ,
        //         para->getParD(level)->numberOfNodes,
        //         para->getParD(level)->isEvenTimestep);
    }
}

void BCKernelManager::runOutflowBCKernelPre(const int level) const{
    if (para->getParD(level)->outflowBC.numberOfBCnodes > 0)
    {
        QPressNoRhoDev27(para->getParD(level).get(), &(para->getParD(level)->outflowBC));

        // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
        // if (  myid == numprocs - 1)
        // PressSchlaffer27(
        //     para->getParD(level)->numberofthreads,
        //     para->getParD(level)->outflowBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->outflowBC.Vx,
        //     para->getParD(level)->outflowBC.Vy,
        //     para->getParD(level)->outflowBC.Vz,
        //     para->getParD(level)->outflowBC.deltaVz,
        //     para->getParD(level)->outflowBC.k,
        //     para->getParD(level)->outflowBC.kN,
        //     para->getParD(level)->outflowBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void BCKernelManager::runPressureBCKernelPre(const int level) const{
    if (para->getParD(level)->pressureBC.numberOfBCnodes > 0)
    {
        this->pressureBoundaryConditionPre(para->getParD(level).get(), &(para->getParD(level)->pressureBC));
    }
}

void BCKernelManager::runPressureBCKernelPost(const int level) const{
    if (para->getParD(level)->pressureBC.numberOfBCnodes > 0)
    {
        // QPressDev27_IntBB(
        //     para->getParD(level)->numberofthreads, 
        //     para->getParD(level)->pressureBC.RhoBC,
        //     para->getParD(level)->distributions.f[0],
        //     para->getParD(level)->pressureBC.k,
        //     para->getParD(level)->pressureBC.q27[0],
        //     para->getParD(level)->pressureBC.numberOfBCnodes,
        //     para->getParD(level)->omega,
        //     para->getParD(level)->neighborX,
        //     para->getParD(level)->neighborY,
        //     para->getParD(level)->neighborZ,
        //     para->getParD(level)->numberOfNodes,
        //     para->getParD(level)->isEvenTimestep);
    }
}

void BCKernelManager::runStressWallModelKernelPost(const int level) const{
    if (para->getParD(level)->stressBC.numberOfBCnodes > 0)
    {
        stressBoundaryConditionPost(para.get(), &(para->getParD(level)->stressBC), level);
    }
}

void BCKernelManager::runSlipBCKernelPost(const int level) const{
    if (para->getParD(level)->slipBC.numberOfBCnodes > 0)
    {
        slipBoundaryConditionPost(para->getParD(level).get(), &(para->getParD(level)->slipBC));
    }
}

void BCKernelManager::runNoSlipBCKernelPost(const int level) const{
    if (para->getParD(level)->noSlipBC.numberOfBCnodes > 0)
    {
        noSlipBoundaryConditionPost(para->getParD(level).get(), &(para->getParD(level)->noSlipBC));
    }
}
