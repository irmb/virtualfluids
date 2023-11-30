#include "GridGenerator.h"

#include <algorithm>

#include <GridGenerator/TransientBCSetter/TransientBCSetter.h>
#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GridGenerator/utilities/math/Math.h>

#include <parallel/Communicator.h>

#include <logger/Logger.h>

#include "GPU/CudaMemoryManager.h"
#include "IndexRearrangementForStreams.h"
#include "InterpolationCellGrouper.h"
#include "LBM/LB.h"
#include "Output/QDebugWriter.hpp"
#include "Parameter/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "utilities/communication.h"

using namespace vf::lbm::dir;

GridGenerator::GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para,
                             std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::parallel::Communicator &communicator)
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
    para->setDistX(distX);
    para->setDistY(distY);
    para->setDistZ(distZ);
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

        cudaMemoryManager->cudaCopyNeighborWSB(level);
        cudaMemoryManager->cudaCopySP(level);
        cudaMemoryManager->cudaCopyCoord(level);
        if(para->getIsBodyForce())
            cudaMemoryManager->cudaCopyBodyForce(level);
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
        {   //TODO: Need to add CollisionTemplate to GridBuilder to allow as argument and get rid of indivual get funtions for fluid node indices... and clean up this mess
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
        builder->sortFluidNodeIndicesAllFeatures(level); //has to be called first!
        builder->sortFluidNodeIndicesMacroVars(level);
        builder->sortFluidNodeIndicesApplyBodyForce(level);
    }
    VF_LOG_INFO("done.");
}

void GridGenerator::allocArrays_BoundaryValues()
{
    VF_LOG_TRACE("-----alloc BoundaryValues------");

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfPressureValues = int(builder->getPressureSize(level));
        VF_LOG_INFO("size pressure level {}: {}", level, numberOfPressureValues);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->pressureBC.numberOfBCnodes = 0;
        para->getParD(level)->outflowPressureCorrectionFactor = para->getOutflowPressureCorrectionFactor();
        if (numberOfPressureValues > 1)
        {
            para->getParH(level)->pressureBC.numberOfBCnodes = numberOfPressureValues;
            cudaMemoryManager->cudaAllocPress(level);
            builder->getPressureValues(para->getParH(level)->pressureBC.RhoBC, para->getParH(level)->pressureBC.k, para->getParH(level)->pressureBC.kN, level);
            cudaMemoryManager->cudaCopyPress(level);
        }
        para->getParD(level)->pressureBC.numberOfBCnodes = para->getParH(level)->pressureBC.numberOfBCnodes;
        para->getParH(level)->numberOfPressureBCnodesRead = para->getParH(level)->pressureBC.numberOfBCnodes * para->getD3Qxx();
        para->getParD(level)->numberOfPressureBCnodesRead = para->getParH(level)->pressureBC.numberOfBCnodes * para->getD3Qxx();
    }

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
        para->getParH(level)->numberOfSlipBCnodesRead = para->getParH(level)->slipBC.numberOfBCnodes * para->getD3Qxx();
        para->getParD(level)->numberOfSlipBCnodesRead = para->getParH(level)->slipBC.numberOfBCnodes * para->getD3Qxx();
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
        para->getParH(level)->numberOfStressBCnodesRead = para->getParH(level)->stressBC.numberOfBCnodes * para->getD3Qxx();
        para->getParD(level)->numberOfStressBCnodesRead = para->getParH(level)->stressBC.numberOfBCnodes * para->getD3Qxx();
    }


    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfVelocityValues = int(builder->getVelocitySize(level));
        VF_LOG_INFO("size velocity level {}: {}", level, numberOfVelocityValues);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        para->getParH(level)->velocityBC.numberOfBCnodes = 0;

        if (numberOfVelocityValues > 1)
        {
            para->getParH(level)->velocityBC.numberOfBCnodes = numberOfVelocityValues;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaAllocVeloBC(level);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            builder->getVelocityValues(para->getParH(level)->velocityBC.Vx, para->getParH(level)->velocityBC.Vy,
                                       para->getParH(level)->velocityBC.Vz, para->getParH(level)->velocityBC.k, level);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            cudaMemoryManager->cudaCopyVeloBC(level);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn()==true){
                //////////////////////////////////////////////////////////////////////////
                para->getParH(level)->TempVel.kTemp = para->getParH(level)->velocityBC.numberOfBCnodes;
                //cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
                VF_LOG_INFO("getTemperatureInit = {}", para->getTemperatureInit());
                VF_LOG_INFO("getTemperatureBC = {}", para->getTemperatureBC());
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempVeloBC(level);
                //cout << "nach alloc " << endl;
                //////////////////////////////////////////////////////////////////////////
                for (uint m = 0; m < para->getParH(level)->velocityBC.numberOfBCnodes; m++)
                {
                    para->getParH(level)->TempVel.temp[m]      = para->getTemperatureInit();
                    para->getParH(level)->TempVel.tempPulse[m] = para->getTemperatureBC();
                    para->getParH(level)->TempVel.velo[m]      = para->getVelocity();
                    para->getParH(level)->TempVel.k[m]         = para->getParH(level)->velocityBC.k[m];
                }
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor copy " << endl;
                cudaMemoryManager->cudaCopyTempVeloBCHD(level);
                //cout << "nach copy " << endl;
                //////////////////////////////////////////////////////////////////////////
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        para->getParD(level)->velocityBC.numberOfBCnodes = para->getParH(level)->velocityBC.numberOfBCnodes;
        para->getParH(level)->numberOfVeloBCnodesRead = para->getParH(level)->velocityBC.numberOfBCnodes * para->getD3Qxx();
        para->getParD(level)->numberOfVeloBCnodesRead = para->getParH(level)->velocityBC.numberOfBCnodes * para->getD3Qxx();
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
        para->getParH(level)->numberOfPrecursorBCnodesRead = numberOfPrecursorValues * para->getD3Qxx();
        para->getParD(level)->numberOfPrecursorBCnodesRead = numberOfPrecursorValues * para->getD3Qxx();

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

            for(auto reader : para->getParH(level)->transientBCInputFileReader)
            {
                if(reader->getNumberOfQuantities() != para->getParD(level)->precursorBC.numberOfQuantities)
                    throw std::runtime_error(
                        "Number of quantities in reader and number of quantities needed for precursor don't match!");
            }

            cudaMemoryManager->cudaCopyPrecursorBC(level);
            cudaMemoryManager->cudaAllocPrecursorData(level);
            para->getParD(level)->precursorBC.streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::Precursor);
            

            // read first timestep of precursor into next and copy to next on device
            for(auto reader : para->getParH(level)->transientBCInputFileReader)
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
            for(auto reader : para->getParH(level)->transientBCInputFileReader)
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
            for(auto reader : para->getParH(level)->transientBCInputFileReader)
            {
                reader->getNextData(para->getParH(level)->precursorBC.next, para->getParH(level)->precursorBC.numberOfPrecursorNodes, 2*nextTime);
            }

            cudaMemoryManager->cudaCopyPrecursorData(level);

            para->getParD(level)->precursorBC.nPrecursorReads = 2;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // advection - diffusion stuff
        if (para->getDiffOn()==true){
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
                if (direction == CommunicationDirections::MX || direction == CommunicationDirections::PX) {
                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);

                    if (tempSend > 0) {
                        int indexProcessNeighbor = (int)para->getParH(level)->sendProcessNeighborX.size();

                        para->getParH(level)->sendProcessNeighborX.emplace_back();
                        para->getParD(level)->sendProcessNeighborX.emplace_back();
                        para->getParH(level)->recvProcessNeighborX.emplace_back();
                        para->getParD(level)->recvProcessNeighborX.emplace_back();
                        if (para->getDiffOn() == true) {
                            para->getParH(level)->sendProcessNeighborADX.emplace_back();
                            para->getParD(level)->sendProcessNeighborADX.emplace_back();
                            para->getParH(level)->recvProcessNeighborADX.emplace_back();
                            para->getParD(level)->recvProcessNeighborADX.emplace_back();
                        }
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        para->getParH(level)->sendProcessNeighborX.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        VF_LOG_INFO("size of Data for X send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborX.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborX.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborX.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParD(level)->sendProcessNeighborX.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParH(level)->sendProcessNeighborX.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborX.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborX.back().memsizeFs = sizeof(real) * tempSend;
                        para->getParD(level)->sendProcessNeighborX.back().memsizeFs = sizeof(real) * tempSend;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for X receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborX.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborX.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborX.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborX.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParD(level)->recvProcessNeighborX.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParH(level)->recvProcessNeighborX.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborX.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborX.back().memsizeFs = sizeof(real) * tempRecv;
                        para->getParD(level)->recvProcessNeighborX.back().memsizeFs = sizeof(real) * tempRecv;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborX(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborX[indexProcessNeighbor].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborX[indexProcessNeighbor].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            indexRearrangement->initCommunicationArraysForCommAfterFinetoCoarseX(level, indexProcessNeighbor, direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborXIndex(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MY || direction == CommunicationDirections::PY) {
                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);

                    if (tempSend > 0) {
                        int indexProcessNeighbor = (int)para->getParH(level)->sendProcessNeighborY.size();

                        para->getParH(level)->sendProcessNeighborY.emplace_back();
                        para->getParD(level)->sendProcessNeighborY.emplace_back();
                        para->getParH(level)->recvProcessNeighborY.emplace_back();
                        para->getParD(level)->recvProcessNeighborY.emplace_back();
                        if (para->getDiffOn() == true) {
                            para->getParH(level)->sendProcessNeighborADY.emplace_back();
                            para->getParD(level)->sendProcessNeighborADY.emplace_back();
                            para->getParH(level)->recvProcessNeighborADY.emplace_back();
                            para->getParD(level)->recvProcessNeighborADY.emplace_back();
                        }
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        VF_LOG_INFO("size of Data for Y send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborY.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborY.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborY.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborY.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParD(level)->sendProcessNeighborY.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParH(level)->sendProcessNeighborY.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborY.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborY.back().memsizeFs = sizeof(real) * tempSend;
                        para->getParD(level)->sendProcessNeighborY.back().memsizeFs = sizeof(real) * tempSend;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for Y receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborY.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborY.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborY.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborY.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParD(level)->recvProcessNeighborY.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParH(level)->recvProcessNeighborY.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborY.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborY.back().memsizeFs = sizeof(real) * tempRecv;
                        para->getParD(level)->recvProcessNeighborY.back().memsizeFs = sizeof(real) * tempRecv;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborY(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborY[indexProcessNeighbor].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborY[indexProcessNeighbor].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            indexRearrangement->initCommunicationArraysForCommAfterFinetoCoarseY(level, indexProcessNeighbor, direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborYIndex(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MZ || direction == CommunicationDirections::PZ) {
                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);

                    if (tempSend > 0) {
                        int indexProcessNeighbor = (int)para->getParH(level)->sendProcessNeighborZ.size();

                        para->getParH(level)->sendProcessNeighborZ.emplace_back();
                        para->getParD(level)->sendProcessNeighborZ.emplace_back();
                        para->getParH(level)->recvProcessNeighborZ.emplace_back();
                        para->getParD(level)->recvProcessNeighborZ.emplace_back();
                        if (para->getDiffOn() == true) {
                            para->getParH(level)->sendProcessNeighborADZ.emplace_back();
                            para->getParD(level)->sendProcessNeighborADZ.emplace_back();
                            para->getParH(level)->recvProcessNeighborADZ.emplace_back();
                            para->getParD(level)->recvProcessNeighborADZ.emplace_back();
                        }
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        VF_LOG_INFO("size of Data for Z send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborZ.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborZ.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborZ.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborZ.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParD(level)->sendProcessNeighborZ.back().numberOfFs    = para->getD3Qxx() * tempSend;
                        para->getParH(level)->sendProcessNeighborZ.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborZ.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborZ.back().memsizeFs = sizeof(real) * tempSend;
                        para->getParD(level)->sendProcessNeighborZ.back().memsizeFs = sizeof(real) * tempSend;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for Z receive buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborZ.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborZ.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborZ.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborZ.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParD(level)->recvProcessNeighborZ.back().numberOfFs    = para->getD3Qxx() * tempRecv;
                        para->getParH(level)->recvProcessNeighborZ.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborZ.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborZ.back().memsizeFs = sizeof(real) * tempRecv;
                        para->getParD(level)->recvProcessNeighborZ.back().memsizeFs = sizeof(real) * tempRecv;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborZ(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborZ[indexProcessNeighbor].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborZ[indexProcessNeighbor].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            indexRearrangement->initCommunicationArraysForCommAfterFinetoCoarseZ(level, indexProcessNeighbor, direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborZIndex(level, indexProcessNeighbor);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }
            }
        }
    }

    // data exchange for F3 / G6
    if ((para->getNumprocs() > 1) && (para->getIsF3())) {
        for (int direction = 0; direction < 6; direction++) {
            if (builder->getCommunicationProcess(direction) == INVALID_INDEX)
                continue;

            for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
                if (direction == CommunicationDirections::MX || direction == CommunicationDirections::PX) {
                    int j = (int)para->getParH(level)->sendProcessNeighborF3X.size();

                    para->getParH(level)->sendProcessNeighborF3X.emplace_back();
                    para->getParD(level)->sendProcessNeighborF3X.emplace_back();
                    para->getParH(level)->recvProcessNeighborF3X.emplace_back();
                    para->getParD(level)->recvProcessNeighborF3X.emplace_back();

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        VF_LOG_INFO("size of Data for X send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3X.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3X.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborF3X.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs    = 6 * tempSend;
                        para->getParD(level)->sendProcessNeighborF3X.back().numberOfGs    = 6 * tempSend;
                        para->getParH(level)->sendProcessNeighborF3X.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborF3X.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborF3X.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs;
                        para->getParD(level)->sendProcessNeighborF3X.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for X recv buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3X.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3X.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborF3X.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs    = 6 * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3X.back().numberOfGs    = 6 * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3X.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3X.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3X.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs;
                        para->getParD(level)->recvProcessNeighborF3X.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborF3X(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3X[j].index, direction,
                                                level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3X[j].index, direction,
                                                   level);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborF3XIndex(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MY || direction == CommunicationDirections::PY) {
                    int j = (int)para->getParH(level)->sendProcessNeighborF3Y.size();

                    para->getParH(level)->sendProcessNeighborF3Y.emplace_back();
                    para->getParD(level)->sendProcessNeighborF3Y.emplace_back();
                    para->getParH(level)->recvProcessNeighborF3Y.emplace_back();
                    para->getParD(level)->recvProcessNeighborF3Y.emplace_back();

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        VF_LOG_INFO("size of Data for Y send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3Y.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3Y.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborF3Y.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs    = 6 * tempSend;
                        para->getParD(level)->sendProcessNeighborF3Y.back().numberOfGs    = 6 * tempSend;
                        para->getParH(level)->sendProcessNeighborF3Y.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborF3Y.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborF3Y.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs;
                        para->getParD(level)->sendProcessNeighborF3Y.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for Y recv buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3Y.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3Y.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Y.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs    = 6 * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Y.back().numberOfGs    = 6 * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Y.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Y.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Y.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs;
                        para->getParD(level)->recvProcessNeighborF3Y.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborF3Y(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3Y[j].index, direction,
                                                level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3Y[j].index, direction,
                                                   level);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborF3YIndex(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MZ || direction == CommunicationDirections::PZ) {
                    int j = (int)para->getParH(level)->sendProcessNeighborF3Z.size();

                    para->getParH(level)->sendProcessNeighborF3Z.emplace_back();
                    para->getParD(level)->sendProcessNeighborF3Z.emplace_back();
                    para->getParH(level)->recvProcessNeighborF3Z.emplace_back();
                    para->getParD(level)->recvProcessNeighborF3Z.emplace_back();

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        VF_LOG_INFO("size of Data for Z send buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempSend, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3Z.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborF3Z.back().numberOfNodes = tempSend;
                        para->getParD(level)->sendProcessNeighborF3Z.back().numberOfNodes = tempSend;
                        para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs    = 6 * tempSend;
                        para->getParD(level)->sendProcessNeighborF3Z.back().numberOfGs    = 6 * tempSend;
                        para->getParH(level)->sendProcessNeighborF3Z.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParD(level)->sendProcessNeighborF3Z.back().memsizeIndex =
                            sizeof(unsigned int) * tempSend;
                        para->getParH(level)->sendProcessNeighborF3Z.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs;
                        para->getParD(level)->sendProcessNeighborF3Z.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // recv
                        VF_LOG_INFO("size of Data for Z recv buffer, \t\tLevel {}: {} \t(neighbor rank: {})", level, tempRecv, builder->getCommunicationProcess(direction));
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3Z.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->recvProcessNeighborF3Z.back().numberOfNodes = tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Z.back().numberOfNodes = tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs    = 6 * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Z.back().numberOfGs    = 6 * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Z.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParD(level)->recvProcessNeighborF3Z.back().memsizeIndex =
                            sizeof(unsigned int) * tempRecv;
                        para->getParH(level)->recvProcessNeighborF3Z.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs;
                        para->getParD(level)->recvProcessNeighborF3Z.back().memsizeGs =
                            sizeof(real) * para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborF3Z(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3Z[j].index, direction,
                                                level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3Z[j].index, direction,
                                                   level);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborF3ZIndex(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }
            }
        }
    }
}

void GridGenerator::allocArrays_BoundaryQs()
{
    VF_LOG_TRACE("allocArrays_BoundaryQs()");


    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const auto numberOfPressureValues = (int)builder->getPressureSize(i);
        if (numberOfPressureValues > 0)
        {
            VF_LOG_INFO("size Pressure: {}: {}", i, numberOfPressureValues);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->pressureBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->pressureBC.numberOfBCnodes;
            QforBoundaryConditions &Q = para->getParH(i)->pressureBC;
            getPointersToBoundaryConditions(Q, QQ, sizeQ);

            builder->getPressureQs(Q.q27, i);


            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion
            if (para->getDiffOn()) {
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor setzen von kTemp" << endl;
                para->getParH(i)->TempPress.kTemp = numberOfPressureValues;
                para->getParD(i)->TempPress.kTemp = numberOfPressureValues;
                VF_LOG_INFO("size TempPress.kTemp: {}: {}", i, para->getParH(i)->TempPress.kTemp);
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempPressBC(i);
                //////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfPressureValues; m++)
                {
                    para->getParH(i)->TempPress.temp[m] = para->getTemperatureInit();
                    para->getParH(i)->TempPress.velo[m] = (real)0.0;
                    para->getParH(i)->TempPress.k[m] = para->getParH(i)->pressureBC.k[m];
                }
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaCopyTempPressBCHD(i);
                //////////////////////////////////////////////////////////////////////////
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaCopyPress(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }//ende if
    }//ende oberste for schleife

    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        int numberOfSlipValues = (int)builder->getSlipSize(i);
        if (numberOfSlipValues > 0)
        {
            VF_LOG_INFO("size Slip:  {}: {}", i, numberOfSlipValues);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->slipBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->slipBC.numberOfBCnodes;
            QforBoundaryConditions &Q = para->getParH(i)->slipBC;
            getPointersToBoundaryConditions(Q, QQ, sizeQ);

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
            real* QQ = para->getParH(i)->stressBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->stressBC.numberOfBCnodes;
            QforBoundaryConditions &Q = para->getParH(i)->stressBC;
            getPointersToBoundaryConditions(Q, QQ, sizeQ);
            
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
            real* QQ = para->getParH(i)->velocityBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->velocityBC.numberOfBCnodes;
            QforBoundaryConditions &Q = para->getParH(i)->velocityBC;
            getPointersToBoundaryConditions(Q, QQ, sizeQ);
            builder->getVelocityQs(Q.q27, i);

            if (para->getDiffOn()) {
                //////////////////////////////////////////////////////////////////////////
                para->getParH(i)->TempVel.kTemp = numberOfVelocityNodes;
                para->getParD(i)->TempVel.kTemp = numberOfVelocityNodes;
                VF_LOG_INFO("size TempVel.kTemp: {}",  para->getParH(i)->TempVel.kTemp);
                VF_LOG_INFO("getTemperatureInit: {}",  para->getTemperatureInit());
                VF_LOG_INFO("getTemperatureBC: {}",  para->getTemperatureBC());
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempVeloBC(i);
                //////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfVelocityNodes; m++)
                {
                    para->getParH(i)->TempVel.temp[m] = para->getTemperatureInit();
                    para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
                    para->getParH(i)->TempVel.velo[m] = para->getVelocity();
                    para->getParH(i)->TempVel.k[m] = para->getParH(i)->velocityBC.k[m];
                }
                cudaMemoryManager->cudaCopyTempVeloBCHD(i);
            }
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
            getPointersToBoundaryConditions(Q, QQ, sizeQ);

            builder->getPrecursorQs(Q.q27, i);

            if (para->getDiffOn()) {
                throw std::runtime_error("Advection diffusion not implemented for Precursor!");
                //////////////////////////////////////////////////////////////////////////
                // para->getParH(i)->TempVel.kTemp = numberOfVelocityNodes;
                // para->getParD(i)->TempVel.kTemp = numberOfVelocityNodes;
                // std::cout << "Groesse TempVel.kTemp = " << para->getParH(i)->TempPress.kTemp << std::endl;
                // std::cout << "getTemperatureInit = " << para->getTemperatureInit() << std::endl;
                // std::cout << "getTemperatureBC = " << para->getTemperatureBC() << std::endl;
                // //////////////////////////////////////////////////////////////////////////
                // cudaMemoryManager->cudaAllocTempVeloBC(i);
                // //cout << "nach alloc " << std::endl;
                // //////////////////////////////////////////////////////////////////////////
                // for (int m = 0; m < numberOfVelocityNodes; m++)
                // {
                //     para->getParH(i)->TempVel.temp[m] = para->getTemperatureInit();
                //     para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
                //     para->getParH(i)->TempVel.velo[m] = para->getVelocity();
                //     para->getParH(i)->TempVel.k[m] = para->getParH(i)->Qinflow.k[m];
                // }
                // //////////////////////////////////////////////////////////////////////////
                // //cout << "vor copy " << std::endl;
                // cudaMemoryManager->cudaCopyTempVeloBCHD(i);
                // //cout << "nach copy " << std::endl;
                //////////////////////////////////////////////////////////////////////////
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
            real* QQ = para->getParH(i)->geometryBC.q27[0];
            unsigned int sizeQ = para->getParH(i)->geometryBC.numberOfBCnodes;
            QforBoundaryConditions &Q = para->getParH(i)->geometryBC;
            getPointersToBoundaryConditions(Q, QQ, sizeQ);
            //////////////////////////////////////////////////////////////////

            builder->getGeometryQs(Q.q27, i);
            //QDebugWriter::writeQValues(Q, para->getParH(i)->geometryBC.k, para->getParH(i)->geometryBC.numberOfBCnodes, "M:/TestGridGeneration/results/GeomGPU.dat");
            //////////////////////////////////////////////////////////////////
            for (int node_i = 0; node_i < numberOfGeometryNodes; node_i++)
            {
                Q.q27[d000][node_i] = 0.0f;
            }
            //for(int test = 0; test < 3; test++)
            //{
            //    for (int tmp = 0; tmp < 27; tmp++)
            //    {
            //        cout <<"Kuhs: " << Q.q27[tmp][test]  << "\n";
            //    }
            //}

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn() == true) {
                    //////////////////////////////////////////////////////////////////////////
                    para->getParH(i)->Temp.kTemp = numberOfGeometryNodes;
                    para->getParD(i)->Temp.kTemp = numberOfGeometryNodes;
                    std::cout << "Groesse Temp.kTemp = " << para->getParH(i)->Temp.kTemp << std::endl;
                    //////////////////////////////////////////////////////////////////////////
                    cudaMemoryManager->cudaAllocTempNoSlipBC(i);
                    //////////////////////////////////////////////////////////////////////////
                    for (int m = 0; m < numberOfGeometryNodes; m++)
                    {
                        para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
                        para->getParH(i)->Temp.k[m] = para->getParH(i)->geometryBC.k[m];
                    }
                    //////////////////////////////////////////////////////////////////////////
                    cudaMemoryManager->cudaCopyTempNoSlipBCHD(i);
                    //////////////////////////////////////////////////////////////////////////
                }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaCopyGeomBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

void GridGenerator::getPointersToBoundaryConditions(QforBoundaryConditions& boundaryConditionStruct, real* subgridDistances, const unsigned int numberOfBCnodes){
    boundaryConditionStruct.q27[dP00] = &subgridDistances[dP00 * numberOfBCnodes];
    boundaryConditionStruct.q27[dM00] = &subgridDistances[dM00 * numberOfBCnodes];
    boundaryConditionStruct.q27[d0P0] = &subgridDistances[d0P0 * numberOfBCnodes];
    boundaryConditionStruct.q27[d0M0] = &subgridDistances[d0M0 * numberOfBCnodes];
    boundaryConditionStruct.q27[d00P] = &subgridDistances[d00P * numberOfBCnodes];
    boundaryConditionStruct.q27[d00M] = &subgridDistances[d00M * numberOfBCnodes];
    boundaryConditionStruct.q27[dPP0] = &subgridDistances[dPP0 * numberOfBCnodes];
    boundaryConditionStruct.q27[dMM0] = &subgridDistances[dMM0 * numberOfBCnodes];
    boundaryConditionStruct.q27[dPM0] = &subgridDistances[dPM0 * numberOfBCnodes];
    boundaryConditionStruct.q27[dMP0] = &subgridDistances[dMP0 * numberOfBCnodes];
    boundaryConditionStruct.q27[dP0P] = &subgridDistances[dP0P * numberOfBCnodes];
    boundaryConditionStruct.q27[dM0M] = &subgridDistances[dM0M * numberOfBCnodes];
    boundaryConditionStruct.q27[dP0M] = &subgridDistances[dP0M * numberOfBCnodes];
    boundaryConditionStruct.q27[dM0P] = &subgridDistances[dM0P * numberOfBCnodes];
    boundaryConditionStruct.q27[d0PP] = &subgridDistances[d0PP * numberOfBCnodes];
    boundaryConditionStruct.q27[d0MM] = &subgridDistances[d0MM * numberOfBCnodes];
    boundaryConditionStruct.q27[d0PM] = &subgridDistances[d0PM * numberOfBCnodes];
    boundaryConditionStruct.q27[d0MP] = &subgridDistances[d0MP * numberOfBCnodes];
    boundaryConditionStruct.q27[d000] = &subgridDistances[d000 * numberOfBCnodes];
    boundaryConditionStruct.q27[dPPP] = &subgridDistances[dPPP * numberOfBCnodes];
    boundaryConditionStruct.q27[dMMP] = &subgridDistances[dMMP * numberOfBCnodes];
    boundaryConditionStruct.q27[dPMP] = &subgridDistances[dPMP * numberOfBCnodes];
    boundaryConditionStruct.q27[dMPP] = &subgridDistances[dMPP * numberOfBCnodes];
    boundaryConditionStruct.q27[dPPM] = &subgridDistances[dPPM * numberOfBCnodes];
    boundaryConditionStruct.q27[dMMM] = &subgridDistances[dMMM * numberOfBCnodes];
    boundaryConditionStruct.q27[dPMM] = &subgridDistances[dPMM * numberOfBCnodes];
    boundaryConditionStruct.q27[dMPM] = &subgridDistances[dMPM * numberOfBCnodes];
}
