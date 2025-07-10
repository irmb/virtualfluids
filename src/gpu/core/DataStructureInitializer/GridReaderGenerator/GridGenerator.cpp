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
//! \addtogroup gpu_DataStructureInitializer DataStructureInitializer
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "GridGenerator.h"

#include <algorithm>

#include <GridGenerator/TransientBCSetter/TransientBCSetter.h>
#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GridGenerator/utilities/math/Math.h>

#include <parallel/Communicator.h>

#include <logger/Logger.h>

#include "Cuda/CudaMemoryManager.h"
#include "IndexRearrangementForStreams.h"
#include "InterpolationCellGrouper.h"
#include "Calculation/Calculation.h"
#include "Output/QDebugWriter.hpp"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "utilities/communication.h"
#include "BoundaryConditions/BoundaryConditionFactory.h"

using namespace vf::lbm::dir;

GridGenerator::GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para,
                             std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::parallel::Communicator& communicator)
    : mpiProcessID(communicator.getProcessID()), builder(builder)
{
    this->para = para;
    this->cudaMemoryManager = cudaMemoryManager;
    this->indexRearrangement = std::make_unique<IndexRearrangementForStreams>(para, builder, communicator);
    this->interpolationGrouper =
        std::make_unique<InterpolationCellGrouper>(para->getParHallLevels(), para->getParDallLevels(), builder);
}

GridGenerator::~GridGenerator() = default;

void GridGenerator::setIndexRearrangementForStreams(std::unique_ptr<IndexRearrangementForStreams> &&indexRearrangement)
{
    this->indexRearrangement = std::move(indexRearrangement);
}

void GridGenerator::initalGridInformations()
{
    if (para->getKernelNeedsFluidNodeIndicesToRun())
        builder->findFluidNodes(para->getUseStreams());
    std::vector<int> gridX, gridY, gridZ;
    std::vector<int> distX, distY, distZ;
    const int numberOfGridLevels = builder->getNumberOfGridLevels();
    builder->getGridInformations(gridX, gridY, gridZ, distX, distY, distZ);
    para->setMaxLevel(numberOfGridLevels);
    para->setGridX(gridX);
    para->setGridY(gridY);
    para->setGridZ(gridZ);
}

void GridGenerator::allocArrays_CoordNeighborGeo()
{
    const uint numberOfLevels = builder->getNumberOfGridLevels();
    VF_LOG_INFO("Number of Level: {}", numberOfLevels);
    int numberOfNodesGlobal = 0;
    VF_LOG_INFO("Number of Nodes: ");

    for (uint level = 0; level < numberOfLevels; level++)
    {
        const uint numberOfNodesPerLevel = builder->getNumberOfNodes(level) + 1;
        numberOfNodesGlobal += numberOfNodesPerLevel;
        VF_LOG_INFO("Level {} = {} Nodes", level, numberOfNodesPerLevel);

        setNumberOfNodes(numberOfNodesPerLevel, level);

        cudaMemoryManager->cudaAllocCoord(level);
        cudaMemoryManager->cudaAllocSP(level);
        //cudaMemoryManager->cudaAllocF3SP(level);
        cudaMemoryManager->cudaAllocNeighborWSB(level);

        if(para->getUseTurbulentViscosity())
            cudaMemoryManager->cudaAllocTurbulentViscosity(level);

        if(para->getIsBodyForce())
            cudaMemoryManager->cudaAllocBodyForce(level);

        if(para->getDiffOn())
        {
            cudaMemoryManager->cudaAllocConcentration(level);
            cudaMemoryManager->cudaAllocConcentrationFs(level);
            if(para->getUseTurbulentDiffusivity())
                cudaMemoryManager->cudaAllocTurbulentDiffusivity(level);
            if(para->getBuoyancyEnabled())
                cudaMemoryManager->cudaAllocLocalReferenceTemperature(level);
        }
        builder->getNodeValues(
            para->getParH(level)->coordinateX,
            para->getParH(level)->coordinateY,
            para->getParH(level)->coordinateZ,
            para->getParH(level)->neighborX,
            para->getParH(level)->neighborY,
            para->getParH(level)->neighborZ,
            para->getParH(level)->neighborInverse,
            para->getParH(level)->typeOfGridNode,
            level);

        setInitialNodeValues(numberOfNodesPerLevel, level);

        if(para->getDiffOn())
            setInitialNodeValuesAD(numberOfNodesPerLevel, level);

        cudaMemoryManager->cudaCopyNeighborWSB(level);
        cudaMemoryManager->cudaCopySP(level);
        cudaMemoryManager->cudaCopyCoord(level);
        if(para->getIsBodyForce())
            cudaMemoryManager->cudaCopyBodyForce(level);
        if(para->getDiffOn())
        {
            cudaMemoryManager->cudaCopyConcentrationHostToDevice(level);
            if(para->getUseTurbulentDiffusivity())
                cudaMemoryManager->cudaCopyTurbulentDiffusivityHostToDevice(level);
            if(para->getBuoyancyEnabled())
                cudaMemoryManager->cudaCopyLocalReferenceTemperatureHostToDevice(level);
        }
    }

    for (int i = 0; i <= para->getMaxLevel(); i++) {
        para->getParH(i)->gridSpacing = builder->getGrid(i)->getDelta();
        para->getParD(i)->gridSpacing = builder->getGrid(i)->getDelta();
    }

    VF_LOG_INFO("Number of Nodes: {}", numberOfNodesGlobal);
    VF_LOG_TRACE("-----finish Coord, Neighbor, Geo------");
}

void GridGenerator::allocArrays_taggedFluidNodes() {

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
    {
        for ( CollisionTemplate tag: all_CollisionTemplate )
        {   //TODO: Need to add CollisionTemplate to GridBuilder to allow as argument and get rid of individual get functions for fluid node indices... and clean up this mess
            switch(tag)
            {
                case CollisionTemplate::Default:
                    this->setNumberOfTaggedFluidNodes(builder->getNumberOfFluidNodes(level), CollisionTemplate::Default, level);
                    cudaMemoryManager->cudaAllocTaggedFluidNodeIndices(CollisionTemplate::Default, level);
                    builder->getFluidNodeIndices(para->getParH(level)->taggedFluidNodeIndices[CollisionTemplate::Default], level);
                    cudaMemoryManager->cudaCopyTaggedFluidNodeIndices(CollisionTemplate::Default, level);
                    if(para->getParH(level)->numberOfTaggedFluidNodes[tag]>0)
                        para->getParH(level)->allocatedBulkFluidNodeTags.push_back(tag);
                    break;
                case CollisionTemplate::SubDomainBorder:
                    this->setNumberOfTaggedFluidNodes(builder->getNumberOfFluidNodesBorder(level), CollisionTemplate::SubDomainBorder, level);
                    cudaMemoryManager->cudaAllocTaggedFluidNodeIndices(CollisionTemplate::SubDomainBorder, level);
                    builder->getFluidNodeIndicesBorder(para->getParH(level)->taggedFluidNodeIndices[CollisionTemplate::SubDomainBorder], level);
                    cudaMemoryManager->cudaCopyTaggedFluidNodeIndices(CollisionTemplate::SubDomainBorder, level);
                    break;
                case CollisionTemplate::WriteMacroVars:
                    this->setNumberOfTaggedFluidNodes(builder->getNumberOfFluidNodesMacroVars(level), CollisionTemplate::WriteMacroVars, level);
                    cudaMemoryManager->cudaAllocTaggedFluidNodeIndices(CollisionTemplate::WriteMacroVars, level);
                    builder->getFluidNodeIndicesMacroVars(para->getParH(level)->taggedFluidNodeIndices[CollisionTemplate::WriteMacroVars], level);
                    cudaMemoryManager->cudaCopyTaggedFluidNodeIndices(CollisionTemplate::WriteMacroVars, level);
                    if(para->getParH(level)->numberOfTaggedFluidNodes[tag]>0)
                        para->getParH(level)->allocatedBulkFluidNodeTags.push_back(tag);
                    break;
                case CollisionTemplate::ApplyBodyForce:
                    this->setNumberOfTaggedFluidNodes(builder->getNumberOfFluidNodesApplyBodyForce(level), CollisionTemplate::ApplyBodyForce, level);
                    cudaMemoryManager->cudaAllocTaggedFluidNodeIndices(CollisionTemplate::ApplyBodyForce, level);
                    builder->getFluidNodeIndicesApplyBodyForce(para->getParH(level)->taggedFluidNodeIndices[CollisionTemplate::ApplyBodyForce], level);
                    cudaMemoryManager->cudaCopyTaggedFluidNodeIndices(CollisionTemplate::ApplyBodyForce, level);
                    if(para->getParH(level)->numberOfTaggedFluidNodes[tag]>0)
                        para->getParH(level)->allocatedBulkFluidNodeTags.push_back(tag);
                    break;
                case CollisionTemplate::AllFeatures:
                    this->setNumberOfTaggedFluidNodes(builder->getNumberOfFluidNodesAllFeatures(level), CollisionTemplate::AllFeatures, level);
                    cudaMemoryManager->cudaAllocTaggedFluidNodeIndices(CollisionTemplate::AllFeatures, level);
                    builder->getFluidNodeIndicesAllFeatures(para->getParH(level)->taggedFluidNodeIndices[CollisionTemplate::AllFeatures], level);
                    cudaMemoryManager->cudaCopyTaggedFluidNodeIndices(CollisionTemplate::AllFeatures, level);
                    if(para->getParH(level)->numberOfTaggedFluidNodes[tag]>0)
                        para->getParH(level)->allocatedBulkFluidNodeTags.push_back(tag);
                    break;
                default:
                    break;
            }
        }
        VF_LOG_INFO("Number of tagged nodes on level {}:", level);
        VF_LOG_INFO("Default: {}, Border: {}, WriteMacroVars: {}, ApplyBodyForce: {}, AllFeatures: {}",
                    para->getParH(level)->numberOfTaggedFluidNodes[CollisionTemplate::Default],
                    para->getParH(level)->numberOfTaggedFluidNodes[CollisionTemplate::SubDomainBorder],
                    para->getParH(level)->numberOfTaggedFluidNodes[CollisionTemplate::WriteMacroVars],
                    para->getParH(level)->numberOfTaggedFluidNodes[CollisionTemplate::ApplyBodyForce],
                    para->getParH(level)->numberOfTaggedFluidNodes[CollisionTemplate::AllFeatures]    );
    }
}

void GridGenerator::tagFluidNodeIndices(const std::vector<uint>& taggedFluidNodeIndices, CollisionTemplate tag, uint level) {
    switch(tag)
    {
        case CollisionTemplate::WriteMacroVars:
            builder->addFluidNodeIndicesMacroVars( taggedFluidNodeIndices, level );
            break;
        case CollisionTemplate::ApplyBodyForce:
            builder->addFluidNodeIndicesApplyBodyForce( taggedFluidNodeIndices, level );
            break;
        case CollisionTemplate::AllFeatures:
            builder->addFluidNodeIndicesAllFeatures( taggedFluidNodeIndices, level );
            break;
        case CollisionTemplate::Default:
        case CollisionTemplate::SubDomainBorder:
            throw std::runtime_error("Cannot tag fluid nodes as Default or SubDomainBorder!");
        default:
            throw std::runtime_error("Tagging fluid nodes with invald tag!");
            break;

    }

}

void GridGenerator::sortFluidNodeTags() {
    VF_LOG_INFO("Start sorting tagged fluid nodes...");
    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
    {
        if(para->getAllNodesAllFeatures() || para->getDiffOn())
            builder->addAllFluidNodeIndicesToAllFeatures(level);
        builder->sortFluidNodeIndicesAllFeatures(level); //has to be called first!
        builder->sortFluidNodeIndicesMacroVars(level);
        builder->sortFluidNodeIndicesApplyBodyForce(level);
    }
    VF_LOG_INFO("done.");
}

void GridGenerator::initPressureBoundaryCondition()
{
    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfPressureValues = int(builder->getPressureSize(level));
        VF_LOG_INFO("size pressure level {}: {}", level, numberOfPressureValues);

        para->getParH(level)->pressureBC.numberOfBCnodes = 0;
        para->getParD(level)->outflowPressureCorrectionFactor = para->getOutflowPressureCorrectionFactor();
        if (numberOfPressureValues > 1) {
            para->getParH(level)->pressureBC.numberOfBCnodes = numberOfPressureValues;
            cudaMemoryManager->cudaAllocPress(level);
            builder->getPressureValues(para->getParH(level)->pressureBC.RhoBC, para->getParH(level)->pressureBC.k,
                                       para->getParH(level)->pressureBC.kN, level);
            cudaMemoryManager->cudaCopyPress(level);
        }
        para->getParD(level)->pressureBC.numberOfBCnodes = para->getParH(level)->pressureBC.numberOfBCnodes;
    }
}

void GridGenerator::initDirectionalPressureBoundaryConditions()
{
    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfBoundaryConditions = builder->getNumberOfPressureBoundaryConditions(level);
        VF_LOG_INFO("Number of pressure boundary conditions on level {}: {}", level, numberOfBoundaryConditions);

        auto& parH = para->getParHostAsReference(level);
        auto& parD = para->getParDeviceAsReference(level);

        parH.pressureBCDirectional.resize(numberOfBoundaryConditions);
        parD.pressureBCDirectional.resize(numberOfBoundaryConditions);

        parD.outflowPressureCorrectionFactor = para->getOutflowPressureCorrectionFactor();

        for (uint bcIndex = 0; bcIndex < numberOfBoundaryConditions; bcIndex++) {
            QforDirectionalBoundaryCondition& bcHost = parH.pressureBCDirectional[bcIndex];
            QforDirectionalBoundaryCondition& bcDevice = parD.pressureBCDirectional[bcIndex];

            const auto numberOfPressureBoundaryNodes = builder->getSizeOfPressureBoundaryCondition(level, bcIndex);
            VF_LOG_INFO("Size of pressure boundary condition {} on level {}: {}", bcIndex, level,
                        static_cast<uint>(numberOfPressureBoundaryNodes));
            bcHost.numberOfBCnodes = static_cast<uint>(numberOfPressureBoundaryNodes);
            bcDevice.numberOfBCnodes = bcHost.numberOfBCnodes;

            bcHost.direction = builder->getPressureBoundaryConditionDirection(level, bcIndex);
            bcDevice.direction = bcHost.direction;

            cudaMemoryManager->cudaAllocDirectionalBoundaryCondition(bcHost, bcDevice);
            builder->getPressureValues(bcHost.RhoBC, bcHost.k, bcHost.kN, level, bcIndex);
            cudaMemoryManager->cudaCopyDirectionalBoundaryCondition(bcHost, bcDevice);
        }
    }
}

void GridGenerator::allocArrays_BoundaryValues(const BoundaryConditionFactory* bcFactory)
{
    VF_LOG_TRACE("-----alloc BoundaryValues------");

    if (bcFactory->hasDirectionalPressureBoundaryCondition())
        initDirectionalPressureBoundaryConditions();
    else
        initPressureBoundaryCondition();

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfSlipValues = int(builder->getSlipSize(level));
        VF_LOG_INFO("size slip level {}: {}", level, numberOfSlipValues);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->slipBC.numberOfBCnodes = 0;
        if (numberOfSlipValues > 1)
        {
            para->getParH(level)->slipBC.numberOfBCnodes = numberOfSlipValues;
            cudaMemoryManager->cudaAllocSlipBC(level);
            builder->getSlipValues(para->getParH(level)->slipBC.normalX, para->getParH(level)->slipBC.normalY,
                                   para->getParH(level)->slipBC.normalZ, para->getParH(level)->slipBC.k, level);
            cudaMemoryManager->cudaCopySlipBC(level);
        }
        para->getParD(level)->slipBC.numberOfBCnodes = para->getParH(level)->slipBC.numberOfBCnodes;
    }

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfStressValues = int(builder->getStressSize(level));
        VF_LOG_INFO("size stress level {}: {}", level, numberOfStressValues);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->stressBC.numberOfBCnodes = 0;
        if (numberOfStressValues > 1)
        {
            para->getParH(level)->stressBC.numberOfBCnodes = numberOfStressValues;
            cudaMemoryManager->cudaAllocStressBC(level);
            cudaMemoryManager->cudaAllocWallModel(level, para->getHasWallModelMonitor());
            builder->getStressValues(   para->getParH(level)->stressBC.normalX,  para->getParH(level)->stressBC.normalY,  para->getParH(level)->stressBC.normalZ,
                                        para->getParH(level)->stressBC.Vx,       para->getParH(level)->stressBC.Vy,       para->getParH(level)->stressBC.Vz,
                                        para->getParH(level)->stressBC.Vx1,      para->getParH(level)->stressBC.Vy1,      para->getParH(level)->stressBC.Vz1,
                                        para->getParH(level)->stressBC.k,        para->getParH(level)->stressBC.kN,
                                        para->getParH(level)->wallModel.samplingOffset, para->getParH(level)->wallModel.z0,
                                        level);

            cudaMemoryManager->cudaCopyStressBC(level);
            cudaMemoryManager->cudaCopyWallModel(level, para->getHasWallModelMonitor());
        }
        para->getParD(level)->stressBC.numberOfBCnodes = para->getParH(level)->stressBC.numberOfBCnodes;
    }


    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfVelocityValues = int(builder->getVelocitySize(level));
        VF_LOG_INFO("size velocity level {}: {}", level, numberOfVelocityValues);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->velocityBC.numberOfBCnodes = numberOfVelocityValues;
        para->getParD(level)->velocityBC.numberOfBCnodes = para->getParH(level)->velocityBC.numberOfBCnodes;
        if (numberOfVelocityValues > 1)
        {
            cudaMemoryManager->cudaAllocVeloBC(level);
            builder->getVelocityValues(para->getParH(level)->velocityBC.Vx, para->getParH(level)->velocityBC.Vy,
                                       para->getParH(level)->velocityBC.Vz, para->getParH(level)->velocityBC.k, level);
            cudaMemoryManager->cudaCopyVeloBC(level);
        }
    }

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfPrecursorValues = int(builder->getPrecursorSize(level));
        VF_LOG_INFO("size precursor level {}: {}", level, numberOfPrecursorValues);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        auto blocks = (numberOfPrecursorValues / para->getParH(level)->numberofthreads) + 1;
        para->getParH(level)->precursorBC.sizeQ = blocks * para->getParH(level)->numberofthreads;
        para->getParD(level)->precursorBC.sizeQ = para->getParH(level)->precursorBC.sizeQ;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->precursorBC.numberOfBCnodes = numberOfPrecursorValues;
        para->getParD(level)->precursorBC.numberOfBCnodes = numberOfPrecursorValues;

        if (numberOfPrecursorValues > 1)
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaAllocPrecursorBC(level);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            builder->getPrecursorValues(
                    para->getParH(level)->precursorBC.planeNeighbor0PP, para->getParH(level)->precursorBC.planeNeighbor0PM,
                    para->getParH(level)->precursorBC.planeNeighbor0MP, para->getParH(level)->precursorBC.planeNeighbor0MM,
                    para->getParH(level)->precursorBC.weights0PP, para->getParH(level)->precursorBC.weights0PM,
                    para->getParH(level)->precursorBC.weights0MP, para->getParH(level)->precursorBC.weights0MM,
                    para->getParH(level)->precursorBC.k, para->getParH(level)->transientBCInputFileReader, para->getParH(level)->precursorBC.numberOfPrecursorNodes,
                    para->getParH(level)->precursorBC.numberOfQuantities, para->getParH(level)->precursorBC.timeStepsBetweenReads,
                    para->getParH(level)->precursorBC.velocityX, para->getParH(level)->precursorBC.velocityY, para->getParH(level)->precursorBC.velocityZ,
                    level);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->getParD(level)->precursorBC.numberOfPrecursorNodes = para->getParH(level)->precursorBC.numberOfPrecursorNodes;
            para->getParD(level)->precursorBC.numberOfQuantities = para->getParH(level)->precursorBC.numberOfQuantities;
            para->getParD(level)->precursorBC.timeStepsBetweenReads = para->getParH(level)->precursorBC.timeStepsBetweenReads;
            para->getParD(level)->precursorBC.velocityX = para->getParH(level)->precursorBC.velocityX;
            para->getParD(level)->precursorBC.velocityY = para->getParH(level)->precursorBC.velocityY;
            para->getParD(level)->precursorBC.velocityZ = para->getParH(level)->precursorBC.velocityZ;

            for(auto& reader : para->getParH(level)->transientBCInputFileReader)
            {
                if(reader->getNumberOfQuantities() != para->getParD(level)->precursorBC.numberOfQuantities)
                    throw std::runtime_error(
                        "Number of quantities in reader and number of quantities needed for precursor don't match!");
            }

            cudaMemoryManager->cudaCopyPrecursorBC(level);
            cudaMemoryManager->cudaAllocPrecursorData(level);
            para->getParD(level)->precursorBC.streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::Precursor);
            

            // read first timestep of precursor into next and copy to next on device
            for(auto& reader : para->getParH(level)->transientBCInputFileReader)
            {
                reader->getNextData(para->getParH(level)->precursorBC.next, para->getParH(level)->precursorBC.numberOfPrecursorNodes, 0);
            }

            cudaMemoryManager->cudaCopyPrecursorData(level);

            //switch next with last pointers
            real* tmp = para->getParD(level)->precursorBC.last;
            para->getParD(level)->precursorBC.last = para->getParD(level)->precursorBC.next;
            para->getParD(level)->precursorBC.next = tmp;

            //read second timestep of precursor into next and copy next to device
            real nextTime = para->getParD(level)->precursorBC.timeStepsBetweenReads*pow(2,-((real)level))*para->getTimeRatio();
            for(auto& reader : para->getParH(level)->transientBCInputFileReader)
            {
                reader->getNextData(para->getParH(level)->precursorBC.next, para->getParH(level)->precursorBC.numberOfPrecursorNodes, nextTime);
            }

            cudaMemoryManager->cudaCopyPrecursorData(level);

            para->getParD(level)->precursorBC.nPrecursorReads = 1;


            //switch next with current pointers
            tmp = para->getParD(level)->precursorBC.current;
            para->getParD(level)->precursorBC.current = para->getParD(level)->precursorBC.next;
            para->getParD(level)->precursorBC.next = tmp;

            //start usual cycle of loading, i.e. read velocities of timestep after current and copy asynchronously to device
            for(auto& reader : para->getParH(level)->transientBCInputFileReader)
            {
                reader->getNextData(para->getParH(level)->precursorBC.next, para->getParH(level)->precursorBC.numberOfPrecursorNodes, 2*nextTime);
            }

            cudaMemoryManager->cudaCopyPrecursorData(level);

            para->getParD(level)->precursorBC.nPrecursorReads = 2;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn())
                throw std::runtime_error(" Advection Diffusion not implemented for Precursor!");
            
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }



    if (builder->hasGeometryValues()) {
        para->setUseGeometryValues(true);
        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
            int numberOfGeometryValues = builder->getGeometrySize(level);
            VF_LOG_INFO("size geometry values, Level {} : {}", level, numberOfGeometryValues);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            para->getParH(level)->geometryBC.numberOfBCnodes = 0;
            if (numberOfGeometryValues > 0)
            {;
                para->getParH(level)->geometryBC.numberOfBCnodes = numberOfGeometryValues;
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocGeomValuesBC(level);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                builder->getGeometryValues(para->getParH(level)->geometryBC.Vx, para->getParH(level)->geometryBC.Vy, para->getParH(level)->geometryBC.Vz, level);

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for (uint m = 0; m < para->getParH(level)->geometryBC.numberOfBCnodes; m++)
                {
                    para->getParH(level)->geometryBC.Vx[m] = para->getParH(level)->geometryBC.Vx[m] / para->getVelocityRatio();
                    para->getParH(level)->geometryBC.Vy[m] = para->getParH(level)->geometryBC.Vy[m] / para->getVelocityRatio();
                    para->getParH(level)->geometryBC.Vz[m] = para->getParH(level)->geometryBC.Vz[m] / para->getVelocityRatio();
                }
                cudaMemoryManager->cudaCopyGeomValuesBC(level);
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// advection - diffusion stuff
                //if (para->getDiffOn()==true){
                //    //////////////////////////////////////////////////////////////////////////
                //    para->getParH(i)->Temp.kTemp = temp4;
                //    cout << "Groesse kTemp = " << para->getParH(i)->Temp.kTemp << "\n";
                //    //////////////////////////////////////////////////////////////////////////
                //    para->cudaAllocTempNoSlipBC(i);
                //    //////////////////////////////////////////////////////////////////////////
                //    for (int m = 0; m < temp4; m++)
                //    {
                //        para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
                //        para->getParH(i)->Temp.k[m]    = para->getParH(i)->geometryBC.k[m];
                //    }
                //    //////////////////////////////////////////////////////////////////////////
                //    para->cudaCopyTempNoSlipBCHD(i);
                //    //////////////////////////////////////////////////////////////////////////
                //}
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
            para->getParD(level)->geometryBC.numberOfBCnodes = para->getParH(level)->geometryBC.numberOfBCnodes;

        }
    }//ende geo

    if (para->getDiffOn()) {
        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
        {
            const int numberOfADNoFluxValues = int(builder->getADNoFluxSize(level));
            para->getParH(level)->AdvectionDiffusionNoFluxBC.numberOfBCnodes = numberOfADNoFluxValues;
            para->getParD(level)->AdvectionDiffusionNoFluxBC.numberOfBCnodes = numberOfADNoFluxValues;
            if(numberOfADNoFluxValues < 1) continue;
            VF_LOG_INFO("size No Flux AD level {}: {}", level, numberOfADNoFluxValues);

            cudaMemoryManager->cudaAllocConcentrationNoFluxBC(level);
            builder->getADNoFluxValues(para->getParH(level)->AdvectionDiffusionNoFluxBC.BCNodeIndices, level);
            cudaMemoryManager->cudaCopyConcentrationNoFluxBCHostToDevice(level);
        }

        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
        {
            const int numberOfADFluxValues = int(builder->getADFluxSize(level));
            para->getParH(level)->AdvectionDiffusionFluxBC.numberOfBCnodes = numberOfADFluxValues;
            para->getParD(level)->AdvectionDiffusionFluxBC.numberOfBCnodes = numberOfADFluxValues;
            
            if(numberOfADFluxValues < 1) continue;
            VF_LOG_INFO("size Flux AD level {}: {}", level, numberOfADFluxValues);

            cudaMemoryManager->cudaAllocConcentrationFluxBC(level);
            builder->getADFluxValues(para->getParH(level)->AdvectionDiffusionFluxBC.normalX,
                                             para->getParH(level)->AdvectionDiffusionFluxBC.normalY,
                                             para->getParH(level)->AdvectionDiffusionFluxBC.normalZ,
                                             para->getParH(level)->AdvectionDiffusionFluxBC.gradient,
                                             para->getParH(level)->AdvectionDiffusionFluxBC.BCNodeIndices, level);
            cudaMemoryManager->cudaCopyConcentrationFluxBCHostToDevice(level);
        }

        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
        {
            const int numberOfADDirichletValues = int(builder->getADDirichletSize(level));
            para->getParH(level)->AdvectionDiffusionDirichletBC.numberOfBCnodes = numberOfADDirichletValues;
            para->getParD(level)->AdvectionDiffusionDirichletBC.numberOfBCnodes = numberOfADDirichletValues;
            if (numberOfADDirichletValues < 1)
                continue;
            VF_LOG_INFO("size Dirichlet AD level {}: {}", level, numberOfADDirichletValues);
            cudaMemoryManager->cudaAllocConcentrationDirichletBC(level);
            builder->getADDirichletValues(para->getParH(level)->AdvectionDiffusionDirichletBC.concentration,
                                          para->getParH(level)->AdvectionDiffusionDirichletBC.vx,
                                          para->getParH(level)->AdvectionDiffusionDirichletBC.vy,
                                          para->getParH(level)->AdvectionDiffusionDirichletBC.vz,
                                          para->getParH(level)->AdvectionDiffusionDirichletBC.BCNodeIndices, level);
            cudaMemoryManager->cudaCopyConcentrationDirichletBCHostToDevice(level);
        }

        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
        {
            const int numberOfNeumannADValues = int(builder->getADNeumannSize(level));
            para->getParH(level)->AdvectionDiffusionNeumannBC.numberOfBCnodes = numberOfNeumannADValues;
            para->getParD(level)->AdvectionDiffusionNeumannBC.numberOfBCnodes = numberOfNeumannADValues;
            if (numberOfNeumannADValues < 1)
                continue;
            VF_LOG_INFO("size Neumann AD level {}: {}", level, numberOfNeumannADValues);

            cudaMemoryManager->cudaAllocConcentrationNeumannBC(level);
            builder->getADNeumannValues(para->getParH(level)->AdvectionDiffusionNeumannBC.gradient,
                                        para->getParH(level)->AdvectionDiffusionNeumannBC.vx,
                                        para->getParH(level)->AdvectionDiffusionNeumannBC.vy,
                                        para->getParH(level)->AdvectionDiffusionNeumannBC.vz,
                                        para->getParH(level)->AdvectionDiffusionNeumannBC.BCNodeIndices, level);
            cudaMemoryManager->cudaCopyConcentrationNeumannBCHostToDevice(level);
        }
    }

    initalValuesDomainDecompostion();
}

void GridGenerator::initalValuesDomainDecompostion()
{
    if (para->getNumprocs() < 2)
        return;
    if ((para->getNumprocs() > 1) /*&& (procNeighborsSendX.size() == procNeighborsRecvX.size())*/) {

        // direction has to be changed in case of periodic BCs and multiple sub domains
        std::vector<int> fillOrder = { 0, 1, 2, 3, 4, 5 };

        for (int direction = 0; direction < 6; direction++) {
            if (direction % 2 > 0 && mpiProcessID % 2 > 0 && (builder->getCommunicationProcess(direction) == builder->getCommunicationProcess(direction - 1)))
            {
                int temp = fillOrder[direction];
                fillOrder[direction] = fillOrder[direction-1];
                fillOrder[direction-1] = temp;
            }
        }

        for (int direction : fillOrder) {
            if (builder->getCommunicationProcess(direction) == INVALID_INDEX)
                continue;

            for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
                auto& parD = para->getParDeviceAsReference(level);
                auto& parH = para->getParHostAsReference(level);
                const uint nSendIndices = builder->getNumberOfSendIndices(direction, level);
                const uint nRecvIndices = builder->getNumberOfReceiveIndices(direction, level);

                if (nSendIndices == 0)
                    continue;

                const uint rankNeighbor = builder->getCommunicationProcess(direction);
                ProcessNeighbor27 sendNeighborHost(nSendIndices, rankNeighbor);
                ProcessNeighbor27 sendNeighborDevice(nSendIndices, rankNeighbor);
                ProcessNeighbor27 recvNeighborHost(nRecvIndices, rankNeighbor);
                ProcessNeighbor27 recvNeighborDevice(nRecvIndices, rankNeighbor);
                cudaMemoryManager->cudaAllocProcessNeighbor(sendNeighborHost, sendNeighborDevice, recvNeighborHost,
                                                            recvNeighborDevice);
                builder->getSendIndices(sendNeighborHost.index, direction, level);
                builder->getReceiveIndices(recvNeighborHost.index, direction, level);

                if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC) {
                    ProcessNeighbor27 sendNeighborAfterFtoCHost;
                    ProcessNeighbor27 recvNeighborAfterFtoCHost;
                    ProcessNeighbor27 sendNeighborAfterFtoCDevice;
                    ProcessNeighbor27 recvNeighborAfterFtoCDevice;
                    indexRearrangement->initCommunicationArraysForCommAfterFinetoCoarse(
                        sendNeighborHost, sendNeighborDevice, sendNeighborAfterFtoCHost, sendNeighborAfterFtoCDevice,
                        recvNeighborHost, recvNeighborDevice, recvNeighborAfterFtoCHost, recvNeighborAfterFtoCDevice, level,
                        direction);
                    switch (direction) {
                        case communication_directions::MX:
                        case communication_directions::PX: {
                            parH.sendProcessNeighborsAfterFtoCX.push_back(sendNeighborAfterFtoCHost);
                            parH.recvProcessNeighborsAfterFtoCX.push_back(recvNeighborAfterFtoCHost);
                            parD.sendProcessNeighborsAfterFtoCX.push_back(sendNeighborAfterFtoCDevice);
                            parD.recvProcessNeighborsAfterFtoCX.push_back(recvNeighborAfterFtoCDevice);
                        } break;
                        case communication_directions::MY:
                        case communication_directions::PY: {
                            parH.sendProcessNeighborsAfterFtoCY.push_back(sendNeighborAfterFtoCHost);
                            parH.recvProcessNeighborsAfterFtoCY.push_back(recvNeighborAfterFtoCHost);
                            parD.sendProcessNeighborsAfterFtoCY.push_back(sendNeighborAfterFtoCDevice);
                            parD.recvProcessNeighborsAfterFtoCY.push_back(recvNeighborAfterFtoCDevice);
                        } break;
                        case communication_directions::MZ:
                        case communication_directions::PZ: {
                            parH.sendProcessNeighborsAfterFtoCZ.push_back(sendNeighborAfterFtoCHost);
                            parH.recvProcessNeighborsAfterFtoCZ.push_back(recvNeighborAfterFtoCHost);
                            parD.sendProcessNeighborsAfterFtoCZ.push_back(sendNeighborAfterFtoCDevice);
                            parD.recvProcessNeighborsAfterFtoCZ.push_back(recvNeighborAfterFtoCDevice);
                        } break;
                    }
                }
                cudaMemoryManager->cudaCopyProcessNeighborIndex(sendNeighborHost, sendNeighborDevice, recvNeighborHost,
                                                                recvNeighborDevice);
                switch (direction) {
                    case communication_directions::MX:
                    case communication_directions::PX: {
                        parH.sendProcessNeighborsX.push_back(sendNeighborHost);
                        parD.sendProcessNeighborsX.push_back(sendNeighborDevice);

                        parH.recvProcessNeighborsX.push_back(recvNeighborHost);
                        parD.recvProcessNeighborsX.push_back(recvNeighborDevice);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        VF_LOG_INFO("size of Data for X send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nSendIndices, rankNeighbor);
                        VF_LOG_INFO("size of Data for X receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nRecvIndices, rankNeighbor);
                    } break;
                    case communication_directions::MY:
                    case communication_directions::PY: {
                        parH.sendProcessNeighborsY.push_back(sendNeighborHost);
                        parD.sendProcessNeighborsY.push_back(sendNeighborDevice);

                        parH.recvProcessNeighborsY.push_back(recvNeighborHost);
                        parD.recvProcessNeighborsY.push_back(recvNeighborDevice);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        VF_LOG_INFO("size of Data for Y send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nSendIndices, rankNeighbor);
                        VF_LOG_INFO("size of Data for Y receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nRecvIndices, rankNeighbor);
                    } break;
                    case communication_directions::MZ:
                    case communication_directions::PZ: {
                        parH.sendProcessNeighborsZ.push_back(sendNeighborHost);
                        parD.sendProcessNeighborsZ.push_back(sendNeighborDevice);

                        parH.recvProcessNeighborsZ.push_back(recvNeighborHost);
                        parD.recvProcessNeighborsZ.push_back(recvNeighborDevice);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        VF_LOG_INFO("size of Data for Z send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nSendIndices, rankNeighbor);
                        VF_LOG_INFO("size of Data for Z receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level,
                                    nRecvIndices, rankNeighbor);
                    } break;
                }
            }
        }
    }
}

void GridGenerator::initSubgridDistancesOfPressureBoundaryCondition(uint level)
{
    VF_LOG_INFO("size Pressure: {}: {}", level, builder->getPressureSize(level));
    QforBoundaryConditions& Q = para->getParH(level)->pressureBC;
    initPointersToSubgridDistances(Q);
    builder->getPressureQs(Q.q27, level);
    cudaMemoryManager->cudaCopyPress(level);
}

void GridGenerator::initSubgridDistancesOfDirectionalPressureBoundaryCondition(uint level)
{
    for (size_t indexInBoundaryConditionVector = 0;
         indexInBoundaryConditionVector < para->getParH(level)->pressureBCDirectional.size();
         indexInBoundaryConditionVector++) {
        QforDirectionalBoundaryCondition& pressureBCHost =
            para->getParH(level)->pressureBCDirectional[indexInBoundaryConditionVector];
        QforDirectionalBoundaryCondition& pressureBCDevice =
            para->getParD(level)->pressureBCDirectional[indexInBoundaryConditionVector];

        VF_LOG_INFO("Size of pressure boundary condition {} on level {}: {}", indexInBoundaryConditionVector, level,
                    static_cast<uint>(pressureBCHost.numberOfBCnodes));

        initPointersToSubgridDistances(pressureBCHost);
        builder->getPressureQs(pressureBCHost.q27, level, static_cast<uint>(indexInBoundaryConditionVector));
        cudaMemoryManager->cudaCopyDirectionalBoundaryCondition(pressureBCHost, pressureBCDevice);
    }
}

void GridGenerator::allocArrays_BoundaryQs()
{
    VF_LOG_TRACE("allocArrays_BoundaryQs()");

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        auto& parH = para->getParHostAsReference(level);
        if (parH.pressureBCDirectional.size() > 0)
            initDirectionalPressureBoundaryConditions();
        else if (builder->getPressureSize(level) > 0)
            initPressureBoundaryCondition();
    }

    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        int numberOfSlipValues = (int)builder->getSlipSize(i);
        if (numberOfSlipValues > 0)
        {
            VF_LOG_INFO("size Slip:  {}: {}", i, numberOfSlipValues);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            QforBoundaryConditions &Q = para->getParH(i)->slipBC;
            initPointersToSubgridDistances(Q);

            builder->getSlipQs(Q.q27, i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaCopySlipBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }//ende if
    }//ende oberste for schleife

    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        int numberOfStressValues = (int)builder->getStressSize(i);
        if (numberOfStressValues > 0)
        {
            VF_LOG_INFO("size Stress:  {}: {}", i, numberOfStressValues);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            QforBoundaryConditions &Q = para->getParH(i)->stressBC;
            initPointersToSubgridDistances(Q);
            
            builder->getStressQs(Q.q27, i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaCopyStressBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }//ende if
    }//ende oberste for schleife

    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const auto numberOfVelocityNodes = int(builder->getVelocitySize(i));
        if (numberOfVelocityNodes > 0)
        {
            VF_LOG_INFO("size velocity level {}: {}", i, numberOfVelocityNodes);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            QforBoundaryConditions &Q = para->getParH(i)->velocityBC;
            initPointersToSubgridDistances(Q);
            builder->getVelocityQs(Q.q27, i);
            cudaMemoryManager->cudaCopyVeloBC(i);
        }
    }

    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const auto numberOfPrecursorNodes = int(builder->getPrecursorSize(i));
        if (numberOfPrecursorNodes > 0)
        {
            VF_LOG_INFO("size velocity level {}: {}", i, numberOfPrecursorNodes);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->precursorBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->precursorBC.numberOfBCnodes;
            QforBoundaryConditions Q;
            getPointersToBoundaryConditions(Q.q27, QQ, sizeQ);

            builder->getPrecursorQs(Q.q27, i);

            if (para->getDiffOn()) {
                throw std::runtime_error("Advection diffusion not implemented for Precursor!");
            }
            cudaMemoryManager->cudaCopyPrecursorBC(i);
        }
    }



    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const int numberOfGeometryNodes = builder->getGeometrySize(i);
        VF_LOG_INFO("size of GeomBoundaryQs, Level {}: {}", i, numberOfGeometryNodes);

        para->getParH(i)->geometryBC.numberOfBCnodes = numberOfGeometryNodes;
        para->getParD(i)->geometryBC.numberOfBCnodes = para->getParH(i)->geometryBC.numberOfBCnodes;
        if (numberOfGeometryNodes > 0)
        {
            //cout << "Groesse der Daten GeomBoundaryQs, Level:  " << i << " : " << numberOfGeometryNodes << "MyID: " << para->getMyID() << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //para->getParH(i)->geometryBC.numberOfBCnodes = temp4;
            //para->getParD(i)->geometryBC.numberOfBCnodes = para->getParH(i)->geometryBC.numberOfBCnodes;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaAllocGeomBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////////////////////
            //Indexarray
            builder->getGeometryIndices(para->getParH(i)->geometryBC.k, i);
            //////////////////////////////////////////////////////////////////////////
            //preprocessing
            QforBoundaryConditions &Q = para->getParH(i)->geometryBC;
            initPointersToSubgridDistances(Q);
            //////////////////////////////////////////////////////////////////

            builder->getGeometryQs(Q.q27, i);
            //QDebugWriter::writeQValues(Q, para->getParH(i)->geometryBC.k, para->getParH(i)->geometryBC.numberOfBCnodes, "M:/TestGridGeneration/results/GeomGPU.dat");
            //////////////////////////////////////////////////////////////////
            for (int node_i = 0; node_i < numberOfGeometryNodes; node_i++)
            {
                Q.q27[d000][node_i] = 0.0f;
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaCopyGeomBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }

    if(para->getDiffOn())
    {
        for(uint level=0; level<builder->getNumberOfGridLevels(); level++)
        {
            const uint numberOfBoundaryNodes = builder->getADNoFluxSize(level);
            if(numberOfBoundaryNodes == 0) continue;
            VF_LOG_INFO("size NoSlipQs, Level {}: {}", level, numberOfBoundaryNodes);

            real* QQ = para->getParH(level)->AdvectionDiffusionNoFluxBC.q27[0];
            getPointersToBoundaryConditions(para->getParH(level)->AdvectionDiffusionNoFluxBC.q27, QQ, numberOfBoundaryNodes);
            builder->getADNoFluxQs(para->getParH(level)->AdvectionDiffusionNoFluxBC.q27, level);
            cudaMemoryManager->cudaCopyConcentrationNoFluxBCHostToDevice(level); 
        }

        for(uint level=0; level<builder->getNumberOfGridLevels(); level++)
        {
            const uint numberOfBoundaryNodes = builder->getADFluxSize(level);
            if(numberOfBoundaryNodes == 0) continue;
            VF_LOG_INFO("size FluxQs, Level {}: {}", level, numberOfBoundaryNodes);
            
            real* QQ = para->getParH(level)->AdvectionDiffusionFluxBC.q27[0];
            getPointersToBoundaryConditions(para->getParH(level)->AdvectionDiffusionFluxBC.q27, QQ, numberOfBoundaryNodes);
            builder->getADFluxQs(para->getParH(level)->AdvectionDiffusionFluxBC.q27, level);
            cudaMemoryManager->cudaCopyConcentrationFluxBCHostToDevice(level); 
        }

        for(uint level=0; level<builder->getNumberOfGridLevels(); level++)
        {
            const uint numberOfBoundaryNodes = builder->getADDirichletSize(level);
            if(numberOfBoundaryNodes == 0) continue;
            VF_LOG_INFO("size DirichletQs, Level {}: {}", level, numberOfBoundaryNodes);
            real* QQ = para->getParH(level)->AdvectionDiffusionDirichletBC.q27[0];
            getPointersToBoundaryConditions(para->getParH(level)->AdvectionDiffusionDirichletBC.q27, QQ, numberOfBoundaryNodes);
            builder->getADDirichletQs(para->getParH(level)->AdvectionDiffusionDirichletBC.q27, level);
            cudaMemoryManager->cudaCopyConcentrationDirichletBCHostToDevice(level); 
        }

        for(uint level=0; level<builder->getNumberOfGridLevels(); level++)
        {
            const uint numberOfBoundaryNodes = builder->getADNeumannSize(level);
            if(numberOfBoundaryNodes == 0) continue;
            VF_LOG_INFO("size NeumannQs, Level {}: {}", level, numberOfBoundaryNodes);
            
            real* QQ = para->getParH(level)->AdvectionDiffusionNeumannBC.q27[0];
            getPointersToBoundaryConditions(para->getParH(level)->AdvectionDiffusionNeumannBC.q27, QQ, numberOfBoundaryNodes);
            builder->getADNeumannQs(para->getParH(level)->AdvectionDiffusionNeumannBC.q27, level);
            cudaMemoryManager->cudaCopyConcentrationNeumannBCHostToDevice(level); 
        }
    }

    VF_LOG_TRACE("-----finish BoundaryQs------");
}

void GridGenerator::allocArrays_OffsetScale()
{
    for (uint level = 0; level < builder->getNumberOfGridLevels() - 1; level++)
    {
        const uint numberOfNodesPerLevelCF = builder->getNumberOfNodesCF(level);
        const uint numberOfNodesPerLevelFC = builder->getNumberOfNodesFC(level);

        VF_LOG_INFO("number of nodes CF Level {}: {}", level, numberOfNodesPerLevelCF);
        VF_LOG_INFO("number of nodes FC Level {}: {}", level, numberOfNodesPerLevelFC);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size CF
        para->getParH(level)->coarseToFine.numberOfCells = numberOfNodesPerLevelCF;
        para->getParD(level)->coarseToFine.numberOfCells = para->getParH(level)->coarseToFine.numberOfCells;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size FC
        para->getParH(level)->fineToCoarse.numberOfCells = numberOfNodesPerLevelFC;
        para->getParD(level)->fineToCoarse.numberOfCells = para->getParH(level)->fineToCoarse.numberOfCells;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //alloc
        cudaMemoryManager->cudaAllocInterfaceCF(level);
        cudaMemoryManager->cudaAllocInterfaceFC(level);
        cudaMemoryManager->cudaAllocInterfaceOffCF(level);
        cudaMemoryManager->cudaAllocInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //init
        builder->getOffsetCF(para->getParH(level)->neighborCoarseToFine.x, para->getParH(level)->neighborCoarseToFine.y, para->getParH(level)->neighborCoarseToFine.z, level);
        builder->getOffsetFC(para->getParH(level)->neighborFineToCoarse.x, para->getParH(level)->neighborFineToCoarse.y, para->getParH(level)->neighborFineToCoarse.z, level);
        builder->getGridInterfaceIndices(para->getParH(level)->coarseToFine.coarseCellIndices, para->getParH(level)->coarseToFine.fineCellIndices, para->getParH(level)->fineToCoarse.coarseCellIndices, para->getParH(level)->fineToCoarse.fineCellIndices, level);

        if (para->getUseStreams() || para->getNumprocs() > 1) {
            // split fine-to-coarse indices into border and bulk
            interpolationGrouper->splitFineToCoarseIntoBorderAndBulk(level);
            // split coarse-to-fine indices into border and bulk
            interpolationGrouper->splitCoarseToFineIntoBorderAndBulk(level);
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //copy
        cudaMemoryManager->cudaCopyInterfaceCF(level);
        cudaMemoryManager->cudaCopyInterfaceFC(level);
        cudaMemoryManager->cudaCopyInterfaceOffCF(level);
        cudaMemoryManager->cudaCopyInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}

void GridGenerator::setDimensions()
{
    //std::vector<int> localGridNX(1);
    //std::vector<int> localGridNY(1);
    //std::vector<int> localGridNZ(1);

    //builder->getDimensions(localGridNX[0], localGridNY[0], localGridNZ[0], 0);

    //para->setGridX(localGridNX);
    //para->setGridY(localGridNY);
    //para->setGridZ(localGridNZ);
}

void GridGenerator::setBoundingBox()
{
    std::vector<int> localGridNX(1);
    std::vector<int> localGridNY(1);
    std::vector<int> localGridNZ(1);
    builder->getDimensions(localGridNX[0], localGridNY[0], localGridNZ[0], 0);

    std::vector<real> minX, maxX, minY, maxY, minZ, maxZ;
    minX.push_back(0);
    minY.push_back(0);
    minZ.push_back(0);

    maxX.push_back((real)localGridNX[0]);
    maxY.push_back((real)localGridNY[0]);
    maxZ.push_back((real)localGridNZ[0]);

    para->setMinCoordX(minX);
    para->setMinCoordY(minY);
    para->setMinCoordZ(minZ);
    para->setMaxCoordX(maxX);
    para->setMaxCoordY(maxY);
    para->setMaxCoordZ(maxZ);
}

void GridGenerator::initPeriodicNeigh(std::vector<std::vector<std::vector<uint> > > periodV, std::vector<std::vector<uint> > periodIndex, std::string way)
{

}





std::string GridGenerator::verifyNeighborIndices(int level) const
{
    std::ostringstream oss;
    oss << "---------report start---------\n";
    oss << "Checking neighbor indices in grid \n";

    int invalidNodes = 0;
    int wrongNeighbors = 0;
    int stopperNodes = 0;

    for (size_t index = 0; index < para->getParH(level)->numberOfNodes; index++)
        oss << verifyNeighborIndex(level, (int)index, invalidNodes, stopperNodes, wrongNeighbors);


    oss << "invalid nodes found: " << invalidNodes << "\n";
    oss << "wrong neighbors found: " << wrongNeighbors << "\n";
    oss << "stopper nodes found : " << stopperNodes << "\n";
    oss << "---------report end---------\n";
    return oss.str();
}

std::string GridGenerator::verifyNeighborIndex(int level, int index , int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const
{
    std::ostringstream oss;

    const int geo = para->getParH(level)->typeOfGridNode[index];
    if (geo == 16)
    {
        stopperNodes++;
        return "";
    }

    real x = para->getParH(level)->coordinateX[index];
    real y = para->getParH(level)->coordinateY[index];
    real z = para->getParH(level)->coordinateZ[index];

    real delta = para->getParH(level)->coordinateX[2] - para->getParH(level)->coordinateX[1];

    //std::cout << para->getParH(level)->coordinateX[1] << ", " << para->getParH(level)->coordinateY[1] << ", " << para->getParH(level)->coordinateZ[1] << std::endl;
    //std::cout << para->getParH(level)->coordinateX[para->getParH(level)->numberOfNodes - 1] << ", " << para->getParH(level)->coordinateY[para->getParH(level)->numberOfNodes - 1] << ", " << para->getParH(level)->coordinateZ[para->getParH(level)->numberOfNodes - 1] << std::endl;

    real maxX = para->getParH(level)->coordinateX[para->getParH(level)->numberOfNodes - 1] - delta;
    real maxY = para->getParH(level)->coordinateY[para->getParH(level)->numberOfNodes - 1] - delta;
    real maxZ = para->getParH(level)->coordinateZ[para->getParH(level)->numberOfNodes - 1] - delta;
    real realNeighborX = vf::Math::lessEqual(x + delta, maxX) ? x + delta : para->getParH(level)->coordinateX[1];
    real realNeighborY = vf::Math::lessEqual(y + delta, maxY) ? y + delta : para->getParH(level)->coordinateY[1];
    real realNeighborZ = vf::Math::lessEqual(z + delta, maxZ) ? z + delta : para->getParH(level)->coordinateZ[1];

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborX[index], realNeighborX, y, z, "X");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborY[index], x, realNeighborY, z, "Y");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ[index], x, y, realNeighborZ, "Z");

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborY[this->para->getParH(level)->neighborX[index]], realNeighborX, realNeighborY, z, "XY");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ[this->para->getParH(level)->neighborX[index]], realNeighborX, y, realNeighborZ, "XZ");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ[this->para->getParH(level)->neighborY[index]], x, realNeighborY, realNeighborZ, "YZ");

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ[this->para->getParH(level)->neighborY[this->para->getParH(level)->neighborX[index]]], realNeighborX, realNeighborY, realNeighborZ, "XYZ");

    return oss.str();
}

std::string GridGenerator::checkNeighbor(int level, real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const
{
    std::ostringstream oss("");
    //if (neighborIndex == -1 || neighborIndex >= size)
    //{
    //    oss << "index broken... \n";
    //    oss << "NeighborX invalid from: (" << x << ", " << y << ", " << z << "), new index: " << newIndex << ", "
    //        << direction << " neighborIndex: " << neighborIndex << "\n";
    //    numberOfWrongNeihgbors++;
    //    return oss.str();
    //}

    real neighborCoordX = para->getParH(level)->coordinateX[neighborIndex];
    real neighborCoordY = para->getParH(level)->coordinateY[neighborIndex];
    real neighborCoordZ = para->getParH(level)->coordinateZ[neighborIndex];

    const bool neighborValid = vf::Math::equal(neighborX, neighborCoordX) && vf::Math::equal(neighborY, neighborCoordY) && vf::Math::equal(neighborZ, neighborCoordZ);

    if (!neighborValid) {
        oss << "NeighborX invalid from: (" << x << ", " << y << ", " << z << "), index: " << index << ", "
            << direction << " neighborIndex: " << neighborIndex <<
            ", actual neighborCoords : (" << neighborCoordX << ", " << neighborCoordY << ", " << neighborCoordZ <<
            "), expected neighborCoords : (" << neighborX << ", " << neighborY << ", " << neighborZ << ")\n";
        numberOfWrongNeihgbors++;
    }
    return oss.str();
}

void GridGenerator::getPointersToBoundaryConditions(real** subgridDistancesInDirections, real* subgridDistances, const unsigned int numberOfBCnodes){
    subgridDistancesInDirections[dP00] = &subgridDistances[dP00 * numberOfBCnodes];
    subgridDistancesInDirections[dM00] = &subgridDistances[dM00 * numberOfBCnodes];
    subgridDistancesInDirections[d0P0] = &subgridDistances[d0P0 * numberOfBCnodes];
    subgridDistancesInDirections[d0M0] = &subgridDistances[d0M0 * numberOfBCnodes];
    subgridDistancesInDirections[d00P] = &subgridDistances[d00P * numberOfBCnodes];
    subgridDistancesInDirections[d00M] = &subgridDistances[d00M * numberOfBCnodes];
    subgridDistancesInDirections[dPP0] = &subgridDistances[dPP0 * numberOfBCnodes];
    subgridDistancesInDirections[dMM0] = &subgridDistances[dMM0 * numberOfBCnodes];
    subgridDistancesInDirections[dPM0] = &subgridDistances[dPM0 * numberOfBCnodes];
    subgridDistancesInDirections[dMP0] = &subgridDistances[dMP0 * numberOfBCnodes];
    subgridDistancesInDirections[dP0P] = &subgridDistances[dP0P * numberOfBCnodes];
    subgridDistancesInDirections[dM0M] = &subgridDistances[dM0M * numberOfBCnodes];
    subgridDistancesInDirections[dP0M] = &subgridDistances[dP0M * numberOfBCnodes];
    subgridDistancesInDirections[dM0P] = &subgridDistances[dM0P * numberOfBCnodes];
    subgridDistancesInDirections[d0PP] = &subgridDistances[d0PP * numberOfBCnodes];
    subgridDistancesInDirections[d0MM] = &subgridDistances[d0MM * numberOfBCnodes];
    subgridDistancesInDirections[d0PM] = &subgridDistances[d0PM * numberOfBCnodes];
    subgridDistancesInDirections[d0MP] = &subgridDistances[d0MP * numberOfBCnodes];
    subgridDistancesInDirections[d000] = &subgridDistances[d000 * numberOfBCnodes];
    subgridDistancesInDirections[dPPP] = &subgridDistances[dPPP * numberOfBCnodes];
    subgridDistancesInDirections[dMMP] = &subgridDistances[dMMP * numberOfBCnodes];
    subgridDistancesInDirections[dPMP] = &subgridDistances[dPMP * numberOfBCnodes];
    subgridDistancesInDirections[dMPP] = &subgridDistances[dMPP * numberOfBCnodes];
    subgridDistancesInDirections[dPPM] = &subgridDistances[dPPM * numberOfBCnodes];
    subgridDistancesInDirections[dMMM] = &subgridDistances[dMMM * numberOfBCnodes];
    subgridDistancesInDirections[dPMM] = &subgridDistances[dPMM * numberOfBCnodes];
    subgridDistancesInDirections[dMPM] = &subgridDistances[dMPM * numberOfBCnodes];
}

void GridGenerator::initPointersToSubgridDistances(QforBoundaryConditions& boundaryCondition)
{
    GridGenerator::getPointersToBoundaryConditions(boundaryCondition.q27, boundaryCondition.q27[0],
                                                   boundaryCondition.numberOfBCnodes);
}

void GridGenerator::initPointersToSubgridDistances(QforDirectionalBoundaryCondition& boundaryCondition)
{
    GridGenerator::getPointersToBoundaryConditions(boundaryCondition.q27, boundaryCondition.q27[0],
                                                   boundaryCondition.numberOfBCnodes);
}
//! \}
