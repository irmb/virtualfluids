#include "GridGenerator.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GPU/CudaMemoryManager.h>

#include <sstream>
#include <iostream>
#include <algorithm>
#include "utilities/math/Math.h"
#include "LBM/LB.h"
#include "Output/QDebugWriter.hpp"

#include "utilities/communication.h"

#include "Communication/Communicator.h"


GridGenerator::GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
	this->builder = builder;
    this->para = para;
    this->cudaMemoryManager = cudaManager;
}

GridGenerator::~GridGenerator()
{

}

void GridGenerator::initalGridInformations()
{
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
	std::cout << "Number of Level: " << numberOfLevels << std::endl;
	int numberOfNodesGlobal = 0;
	std::cout << "Number of Nodes: " << std::endl;
	
	for (uint level = 0; level < numberOfLevels; level++) 
	{
		const int numberOfNodesPerLevel = builder->getNumberOfNodes(level) + 1;
		numberOfNodesGlobal += numberOfNodesPerLevel;
		std::cout << "Level " << level << " = " << numberOfNodesPerLevel << " Nodes" << std::endl;
	
		setNumberOfNodes(numberOfNodesPerLevel, level);
	
		cudaMemoryManager->cudaAllocCoord(level);
        cudaMemoryManager->cudaAllocSP(level);
        //cudaMemoryManager->cudaAllocF3SP(level);
		cudaMemoryManager->cudaAllocNeighborWSB(level);

        if(para->getUseWale())
            cudaMemoryManager->cudaAllocTurbulentViscosity(level);

		builder->getNodeValues(
			para->getParH(level)->coordX_SP,
			para->getParH(level)->coordY_SP,
			para->getParH(level)->coordZ_SP,
			para->getParH(level)->neighborX_SP,
			para->getParH(level)->neighborY_SP,
			para->getParH(level)->neighborZ_SP,
			para->getParH(level)->neighborWSB_SP,
			para->getParH(level)->geoSP,
			level);

		setInitalNodeValues(numberOfNodesPerLevel, level);

        cudaMemoryManager->cudaCopyNeighborWSB(level);
        cudaMemoryManager->cudaCopySP(level);
        cudaMemoryManager->cudaCopyCoord(level);

        //std::cout << verifyNeighborIndices(level);
	}
	std::cout << "Number of Nodes: " << numberOfNodesGlobal << std::endl;
	std::cout << "-----finish Coord, Neighbor, Geo------" << std::endl;
}

void GridGenerator::allocArrays_fluidNodeIndices() {
    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        setNumberOfFluidNodes(builder->getNumberOfFluidNodes(level), level);
        cudaMemoryManager->cudaAllocFluidNodeIndices(level);
        builder->getFluidNodeIndices(para->getParH(level)->fluidNodeIndices, level);
        cudaMemoryManager->cudaCopyFluidNodeIndices(level);
    }    
}

void GridGenerator::allocArrays_fluidNodeIndicesBorder() {
    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        setNumberOfFluidNodesBorder(builder->getNumberOfFluidNodesBorder(level), level);
        cudaMemoryManager->cudaAllocFluidNodeIndicesBorder(level);
        builder->getFluidNodeIndicesBorder(para->getParH(level)->fluidNodeIndicesBorder, level);
        cudaMemoryManager->cudaCopyFluidNodeIndicesBorder(level);
    }
}

void GridGenerator::allocArrays_BoundaryValues()
{
	std::cout << "------read BoundaryValues------" << std::endl;

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfPressureValues = int(builder->getPressureSize(level));

        std::cout << "size pressure level " << level << " : " << numberOfPressureValues << std::endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->QPress.kQ = numberOfPressureValues;
        para->getParD(level)->QPress.kQ = numberOfPressureValues;
        para->getParH(level)->kPressQread = numberOfPressureValues * para->getD3Qxx();
        para->getParD(level)->kPressQread = numberOfPressureValues * para->getD3Qxx();
        if (numberOfPressureValues > 1)
        {
            cudaMemoryManager->cudaAllocPress(level);
            builder->getPressureValues(para->getParH(level)->QPress.RhoBC, para->getParH(level)->QPress.k, para->getParH(level)->QPress.kN, level);
            cudaMemoryManager->cudaCopyPress(level);
        }
    }
    

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfVelocityValues = int(builder->getVelocitySize(level));
        std::cout << "size velocity level " << level << " : " << numberOfVelocityValues << std::endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int blocks = (numberOfVelocityValues / para->getParH(level)->numberofthreads) + 1;
        para->getParH(level)->Qinflow.kArray = blocks * para->getParH(level)->numberofthreads;
        para->getParD(level)->Qinflow.kArray = para->getParH(level)->Qinflow.kArray;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->Qinflow.kQ = numberOfVelocityValues;
        para->getParD(level)->Qinflow.kQ = numberOfVelocityValues;
        para->getParH(level)->kInflowQ = numberOfVelocityValues;
        para->getParD(level)->kInflowQ = numberOfVelocityValues;
        para->getParH(level)->kInflowQread = numberOfVelocityValues * para->getD3Qxx();
        para->getParD(level)->kInflowQread = numberOfVelocityValues * para->getD3Qxx();

        if (numberOfVelocityValues > 1)
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaAllocVeloBC(level);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            builder->getVelocityValues(para->getParH(level)->Qinflow.Vx, para->getParH(level)->Qinflow.Vy, para->getParH(level)->Qinflow.Vz, para->getParH(level)->Qinflow.k, level);


            //for (int i = 0; i < numberOfVelocityValues; i++)
            //{
            //    std::cout << "index: " << para->getParH(level)->Qinflow.k[i];
            //    std::cout << " (x,y,z)" << para->getParH(level)->coordX_SP[para->getParH(level)->Qinflow.k[i]];
            //    std::cout << ", " << para->getParH(level)->coordY_SP[para->getParH(level)->Qinflow.k[i]];
            //    std::cout << ", " << para->getParH(level)->coordZ_SP[para->getParH(level)->Qinflow.k[i]];
            //    std::cout << " geo: " << para->getParH(level)->geoSP[para->getParH(level)->Qinflow.k[i]];
            //    std::cout << std::endl;
            //}


            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            cudaMemoryManager->cudaCopyVeloBC(level);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn()==true){
            	//////////////////////////////////////////////////////////////////////////
            	para->getParH(level)->TempVel.kTemp = numberOfVelocityValues;
            	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
            	std::cout << "getTemperatureInit = " << para->getTemperatureInit() << std::endl;
            	std::cout << "getTemperatureBC = " << para->getTemperatureBC() << std::endl;
            	//////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempVeloBC(level);
            	//cout << "nach alloc " << endl;
            	//////////////////////////////////////////////////////////////////////////
            	for (int m = 0; m < numberOfVelocityValues; m++)
            	{
            		para->getParH(level)->TempVel.temp[m]      = para->getTemperatureInit();
            		para->getParH(level)->TempVel.tempPulse[m] = para->getTemperatureBC();
            		para->getParH(level)->TempVel.velo[m]      = para->getVelocity();
            		para->getParH(level)->TempVel.k[m]         = para->getParH(level)->Qinflow.k[m];
            	}
            	//////////////////////////////////////////////////////////////////////////
            	//cout << "vor copy " << endl;
                cudaMemoryManager->cudaCopyTempVeloBCHD(level);
            	//cout << "nach copy " << endl;
            	//////////////////////////////////////////////////////////////////////////
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }



    if (builder->hasGeometryValues()) {
        para->setGeometryValues(true);
        for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
            int numberOfGeometryValues = builder->getGeometrySize(i);
            std::cout << "size geometry values, Level " << i << " : " << numberOfGeometryValues << std::endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->getParH(i)->QGeom.kQ = numberOfGeometryValues;
            para->getParD(i)->QGeom.kQ = numberOfGeometryValues;
            if (numberOfGeometryValues > 0)
            {

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocGeomValuesBC(i);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Indexarray

                builder->getGeometryValues(para->getParH(i)->QGeom.Vx, para->getParH(i)->QGeom.Vy, para->getParH(i)->QGeom.Vz, i);

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfGeometryValues; m++)
                {
                    para->getParH(i)->QGeom.Vx[m] = para->getParH(i)->QGeom.Vx[m] / para->getVelocityRatio();
                    para->getParH(i)->QGeom.Vy[m] = para->getParH(i)->QGeom.Vy[m] / para->getVelocityRatio();
                    para->getParH(i)->QGeom.Vz[m] = para->getParH(i)->QGeom.Vz[m] / para->getVelocityRatio();
                    //para->getParH(i)->QGeom.Vx[m] = para->getParH(i)->QGeom.Vx[m] / 100.0f;
                    //para->getParH(i)->QGeom.Vy[m] = para->getParH(i)->QGeom.Vy[m] / 100.0f;
                    //para->getParH(i)->QGeom.Vz[m] = para->getParH(i)->QGeom.Vz[m] / 100.0f;
                    //para->getParH(i)->QGeom.Vx[m] = 0.0f;
                    //para->getParH(i)->QGeom.Vy[m] = 0.0f;
                    //para->getParH(i)->QGeom.Vz[m] = 0.0f;
                }
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ////T�st
                //for (int m = 0; m < temp4; m++)
                //{
                //	para->getParH(i)->QGeom.Vx[m] = para->getVelocity();//0.035f;
                //	para->getParH(i)->QGeom.Vy[m] = 0.0f;//para->getVelocity();//0.0f;
                //	para->getParH(i)->QGeom.Vz[m] = 0.0f;
                //}
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaCopyGeomValuesBC(i);
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// advection - diffusion stuff
                //if (para->getDiffOn()==true){
                //	//////////////////////////////////////////////////////////////////////////
                //	para->getParH(i)->Temp.kTemp = temp4;
                //	cout << "Groesse kTemp = " << para->getParH(i)->Temp.kTemp << std::endl;
                //	//////////////////////////////////////////////////////////////////////////
                //	para->cudaAllocTempNoSlipBC(i);
                //	//////////////////////////////////////////////////////////////////////////
                //	for (int m = 0; m < temp4; m++)
                //	{
                //		para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
                //		para->getParH(i)->Temp.k[m]    = para->getParH(i)->QGeom.k[m];
                //	}
                //	//////////////////////////////////////////////////////////////////////////
                //	para->cudaCopyTempNoSlipBCHD(i);
                //	//////////////////////////////////////////////////////////////////////////
                //}
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
    }//ende geo

    initalValuesDomainDecompostion();
}

void GridGenerator::initalValuesDomainDecompostion()
{
    if (para->getNumprocs() < 2)
        return;
    if ((para->getNumprocs() > 1) /*&& (procNeighborsSendX.size() == procNeighborsRecvX.size())*/) {
        for (int direction = 0; direction < 6; direction++) {
            if (builder->getCommunicationProcess(direction) == INVALID_INDEX)
                continue;

            for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
                if (direction == CommunicationDirections::MX || direction == CommunicationDirections::PX) {
                    int j = (int)para->getParH(level)->sendProcessNeighborX.size();

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

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
                        ////////////////////////////////////////////////////////////////////////////////////////
                        para->getParH(level)->sendProcessNeighborX.back().rankNeighbor =
                            builder->getCommunicationProcess(direction);
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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
                        cudaMemoryManager->cudaAllocProcessNeighborX(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborX[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborX[j].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            initCommunicationArraysForCommAfterFinetoCoarseX(level, j, direction);                        
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborXIndex(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MY || direction == CommunicationDirections::PY) {
                    int j = (int)para->getParH(level)->sendProcessNeighborY.size();

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

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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
                        cudaMemoryManager->cudaAllocProcessNeighborY(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////                        
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborY[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborY[j].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            initCommunicationArraysForCommAfterFinetoCoarseY(level, j, direction);                       
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborYIndex(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                    }
                }

                if (direction == CommunicationDirections::MZ || direction == CommunicationDirections::PZ) {
                    int j = (int)para->getParH(level)->sendProcessNeighborZ.size();

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

                    int tempSend = builder->getNumberOfSendIndices(direction, level);
                    int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
                    if (tempSend > 0) {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // send
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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
                        cudaMemoryManager->cudaAllocProcessNeighborZ(level, j);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // init index arrays
                        builder->getSendIndices(para->getParH(level)->sendProcessNeighborZ[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborZ[j].index, direction,
                                                   level);
                        if (level != builder->getNumberOfGridLevels() - 1 && para->useReducedCommunicationAfterFtoC)
                            initCommunicationArraysForCommAfterFinetoCoarseZ(level, j, direction);
                        ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborZIndex(level, j);
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
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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
                        std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend
                                  << std::endl;
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
                        std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv
                                  << std::endl;
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

void GridGenerator::initCommunicationArraysForCommAfterFinetoCoarseX(const uint &level, int j, int direction)
{
    // init send indices for communication after coarse to fine
    std::cout << "communication: reorder send indices X ";
    para->initNumberOfProcessNeighborsAfterFtoCX(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCX(direction, level, j, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCX(para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].numberOfNodes,
                                            level, j);

    // send sendIndicesForCommAfterFtoCPositions to receiving process and receive recvIndicesForCommAfterFtoCPositions from sending process
    std::cout << "mpi send and receive ";
    std::vector<uint> recvIndicesForCommAfterFtoCPositions;
    recvIndicesForCommAfterFtoCPositions.resize(
        (size_t)para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].numberOfNodes *
        2); // give vector an arbitraty size (larger than needed) // TODO: This is stupid! Find a better way
    auto comm = vf::gpu::Communicator::getInstanz();
    comm->exchangeIndices(recvIndicesForCommAfterFtoCPositions.data(), recvIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->recvProcessNeighborX[j].rankNeighbor,
                          sendIndicesForCommAfterFtoCPositions.data(), sendIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->sendProcessNeighborX[j].rankNeighbor);
    // resize receiving vector to correct size
    auto it = std::unique(recvIndicesForCommAfterFtoCPositions.begin(), recvIndicesForCommAfterFtoCPositions.end());
    recvIndicesForCommAfterFtoCPositions.erase(std::prev(it, 1), recvIndicesForCommAfterFtoCPositions.end());

    // init receive indices for communication after coarse to fine
    std::cout << "reorder receive indices ";
    reorderRecvIndicesForCommAfterFtoCX(direction, level, j, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCX(para->getParH(level)->recvProcessNeighborsAfterFtoCX[j].numberOfNodes,
                                            level, j);
    copyProcessNeighborToAfterFtoCX(level, j);

    std::cout << "done." << std::endl;
}

void GridGenerator::initCommunicationArraysForCommAfterFinetoCoarseY(const uint &level, int j, int direction)
{
    // init send indices for communication after coarse to fine
    std::cout << "communication: reorder send indices Y ";
    para->initNumberOfProcessNeighborsAfterFtoCY(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCY(direction, level, j, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCY(para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].numberOfNodes,
                                            level, j);

    // send sendIndicesForCommAfterFtoCPositions to receiving process and receive recvIndicesForCommAfterFtoCPositions from sending process
    std::cout << "mpi send and receive ";
    std::vector<uint> recvIndicesForCommAfterFtoCPositions; 
    recvIndicesForCommAfterFtoCPositions.resize((size_t) para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].numberOfNodes *
                                                2); // give vector an arbitraty size (larger than needed) // TODO: This is stupid! Find a better way
    auto comm = vf::gpu::Communicator::getInstanz();
    comm->exchangeIndices(recvIndicesForCommAfterFtoCPositions.data(), recvIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->recvProcessNeighborY[j].rankNeighbor,
                          sendIndicesForCommAfterFtoCPositions.data(), sendIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->sendProcessNeighborY[j].rankNeighbor);
    // resize receiving vector to correct size
    auto it = std::unique(recvIndicesForCommAfterFtoCPositions.begin(), recvIndicesForCommAfterFtoCPositions.end());
    recvIndicesForCommAfterFtoCPositions.erase(std::prev(it, 1), recvIndicesForCommAfterFtoCPositions.end());

    // init receive indices for communication after coarse to fine
    std::cout << "reorder receive indices ";
    reorderRecvIndicesForCommAfterFtoCY(direction, level, j, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCY(para->getParH(level)->recvProcessNeighborsAfterFtoCY[j].numberOfNodes,
                                            level, j);

    copyProcessNeighborToAfterFtoCY(level, j);

    std::cout << "done." << std::endl;
}

void GridGenerator::initCommunicationArraysForCommAfterFinetoCoarseZ(const uint &level, int j, int direction)
{
    // init send indices for communication after coarse to fine
    std::cout << "communication: reorder send indices Z ";
    para->initNumberOfProcessNeighborsAfterFtoCZ(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCZ(direction, level, j, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCZ(para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].numberOfNodes,
                                            level, j);

    // send sendIndicesForCommAfterFtoCPositions to receiving process and receive recvIndicesForCommAfterFtoCPositions from sending process
    std::cout << "mpi send and receive ";
    std::vector<uint> recvIndicesForCommAfterFtoCPositions; 
    recvIndicesForCommAfterFtoCPositions.resize((size_t) para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].numberOfNodes *
                                                2); // give vector an arbitraty size (larger than needed) // TODO: This is stupid! Find a better way
    auto comm = vf::gpu::Communicator::getInstanz();
    comm->exchangeIndices(recvIndicesForCommAfterFtoCPositions.data(), recvIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->recvProcessNeighborZ[j].rankNeighbor,
                          sendIndicesForCommAfterFtoCPositions.data(), sendIndicesForCommAfterFtoCPositions.size(),
                          para->getParH(level)->sendProcessNeighborZ[j].rankNeighbor);
    // resize receiving vector to correct size
    auto it = std::unique(recvIndicesForCommAfterFtoCPositions.begin(), recvIndicesForCommAfterFtoCPositions.end());
    recvIndicesForCommAfterFtoCPositions.erase(std::prev(it, 1), recvIndicesForCommAfterFtoCPositions.end());

    // init receive indices for communication after coarse to fine
    std::cout << "reorder receive indices ";
    reorderRecvIndicesForCommAfterFtoCZ(direction, level, j, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCZ(para->getParH(level)->recvProcessNeighborsAfterFtoCZ[j].numberOfNodes,
                                            level, j);

    copyProcessNeighborToAfterFtoCZ(level, j);

    std::cout << "done." << std::endl;
}

void GridGenerator::copyProcessNeighborToAfterFtoCX(const uint &level, int j)
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCX[j].f[0] = para->getParD(level)->sendProcessNeighborX[j].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].f[0] = para->getParH(level)->sendProcessNeighborX[j].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCX[j].f[0] = para->getParD(level)->recvProcessNeighborX[j].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[j].f[0] = para->getParH(level)->recvProcessNeighborX[j].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCX[j].index = para->getParD(level)->sendProcessNeighborX[j].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].index = para->getParH(level)->sendProcessNeighborX[j].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCX[j].index = para->getParD(level)->recvProcessNeighborX[j].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[j].index = para->getParH(level)->recvProcessNeighborX[j].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].rankNeighbor = para->getParH(level)->sendProcessNeighborX[j].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[j].rankNeighbor = para->getParH(level)->recvProcessNeighborX[j].rankNeighbor;
}

void GridGenerator::copyProcessNeighborToAfterFtoCY(const uint &level, int j)
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCY[j].f[0] = para->getParD(level)->sendProcessNeighborY[j].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].f[0] = para->getParH(level)->sendProcessNeighborY[j].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCY[j].f[0] = para->getParD(level)->recvProcessNeighborY[j].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[j].f[0] = para->getParH(level)->recvProcessNeighborY[j].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCY[j].index = para->getParD(level)->sendProcessNeighborY[j].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].index = para->getParH(level)->sendProcessNeighborY[j].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCY[j].index = para->getParD(level)->recvProcessNeighborY[j].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[j].index = para->getParH(level)->recvProcessNeighborY[j].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].rankNeighbor = para->getParH(level)->sendProcessNeighborY[j].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[j].rankNeighbor = para->getParH(level)->recvProcessNeighborY[j].rankNeighbor;
}

void GridGenerator::copyProcessNeighborToAfterFtoCZ(const uint &level, int j)
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCZ[j].f[0] = para->getParD(level)->sendProcessNeighborZ[j].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].f[0] = para->getParH(level)->sendProcessNeighborZ[j].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCZ[j].f[0] = para->getParD(level)->recvProcessNeighborZ[j].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[j].f[0] = para->getParH(level)->recvProcessNeighborZ[j].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCZ[j].index = para->getParD(level)->sendProcessNeighborZ[j].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].index = para->getParH(level)->sendProcessNeighborZ[j].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCZ[j].index = para->getParD(level)->recvProcessNeighborZ[j].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[j].index = para->getParH(level)->recvProcessNeighborZ[j].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].rankNeighbor = para->getParH(level)->sendProcessNeighborZ[j].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[j].rankNeighbor = para->getParH(level)->recvProcessNeighborZ[j].rankNeighbor;
}

void GridGenerator::reorderSendIndicesForCommAfterFtoCX(int direction, int level, int j,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    int *sendIndices                    = para->getParH(level)->sendProcessNeighborX[j].index;
    int &numberOfSendNeighborsAfterFtoC = para->getParH(level)->sendProcessNeighborsAfterFtoCX[j].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderSendIndicesForCommAfterFtoCY(int direction, int level, int j,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    int *sendIndices                    = para->getParH(level)->sendProcessNeighborY[j].index;
    int &numberOfSendNeighborsAfterFtoC = para->getParH(level)->sendProcessNeighborsAfterFtoCY[j].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderSendIndicesForCommAfterFtoCZ(int direction, int level, int j,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions) 
{
    int *sendIndices                    = para->getParH(level)->sendProcessNeighborZ[j].index;
    int &numberOfSendNeighborsAfterFtoC = para->getParH(level)->sendProcessNeighborsAfterFtoCZ[j].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderSendIndicesForCommAfterFtoC(int *sendIndices, int &numberOfSendNeighborsAfterFtoC,
                                                       int direction, int level, int j,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "reorder send indices for communication after fine to coarse: level: " << level
                  << " direction: " << direction;
    if (para->getParH(level)->K_CF == 0 || para->getParH(level)->K_FC == 0)
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderSendIndicesForCommAfterFtoC(): iCellFCC needs to be inititalized before calling "
                         "this function "
                      << "\n";

    int sparseIndexSend;
    std::vector<int> sendIndicesAfterFtoC;
    std::vector<int> sendIndicesOther;
    std::array<int, 7> neighbors;
    uint numberOfSendIndices = builder->getNumberOfSendIndices(direction, level);

    //iCellFCC
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        neighbors.fill(-1);
        sparseIndexSend = sendIndices[posInSendIndices];  
        if (isSparseIndexInICellFCC(para->getParH(level)->K_CF, sparseIndexSend, level))
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndexSend,
                                                 sendIndicesForCommAfterFtoCPositions, posInSendIndices);
    }

    // iCellCFC
    std::vector<uint> nodesCFC;
    aggregateNodesInICellCFC(level, nodesCFC);
    for (auto sparseIndex : nodesCFC)
        findIfSparseIndexIsInSendIndicesAndAddToCommVectors(sparseIndex, sendIndices, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);

    numberOfSendNeighborsAfterFtoC = (int)sendIndicesAfterFtoC.size();

    findIndicesNotInCommAfterFtoC(numberOfSendIndices, sendIndices, sendIndicesAfterFtoC, sendIndicesOther);

    // copy new vectors back to sendIndices array
    for (int i = 0; i < numberOfSendNeighborsAfterFtoC; i++)
        sendIndices[i] = sendIndicesAfterFtoC[i];
    for (uint i = 0; i < (uint)sendIndicesOther.size(); i++)
        sendIndices[i + numberOfSendNeighborsAfterFtoC] = sendIndicesOther[i];

    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "... numberOfSendNeighborsAfterFtoC: " << numberOfSendNeighborsAfterFtoC << "\n";

    if (numberOfSendNeighborsAfterFtoC + sendIndicesOther.size() != numberOfSendIndices) {
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderSendIndicesForCommAfterFtoC(): incorrect number of nodes"
                      << "\n";
        std::cout << "numberOfSendNeighborsAfterFtoC = " << numberOfSendNeighborsAfterFtoC
                  << ", sendOrIndicesOther.size() = " << sendIndicesOther.size()
                  << ", numberOfSendOrRecvIndices = " << numberOfSendIndices << std::endl;
    }
}

bool GridGenerator::isSparseIndexInICellFCC(uint sizeOfICellFCC, int sparseIndex, int level)
{
    for (uint j = 0; j < sizeOfICellFCC; j++) {
        if (sparseIndex < 0)
            return false;
        if (para->getParH(level)->intFC.ICellFCC[j] == (uint)sparseIndex) {
            return true;
        }
    }
    return false;
}

void GridGenerator::aggregateNodesInICellCFC(int level, std::vector<uint> &nodesCFC)
{
    uint sparseIndex;
    uint *neighborX = para->getParH(level)->neighborX_SP;
    uint *neighborY = para->getParH(level)->neighborY_SP;
    uint *neighborZ = para->getParH(level)->neighborZ_SP;

    for (uint x = 0; x < para->getParH(level)->K_FC; x++) {
        sparseIndex = para->getParH(level)->intCF.ICellCFC[x];
        nodesCFC.push_back(sparseIndex);
        nodesCFC.push_back(neighborX[sparseIndex]);
        nodesCFC.push_back(neighborY[sparseIndex]);
        nodesCFC.push_back(neighborZ[sparseIndex]);
        nodesCFC.push_back(neighborY[neighborX[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborX[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborY[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborY[neighborX[sparseIndex]]]);           
    }
    std::sort(nodesCFC.begin(), nodesCFC.end());
    auto iterator = std::unique(nodesCFC.begin(), nodesCFC.end());
    nodesCFC.erase(iterator, nodesCFC.end());
}

void GridGenerator::addUniqueIndexToCommunicationVectors(
    std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
    std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions, uint &posInSendIndices) const
{
    // add index to corresponding vectors but omit indices which are already in sendIndicesAfterFtoC
    if (std::find(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), sparseIndexSend) == sendIndicesAfterFtoC.end()) {
        sendIndicesAfterFtoC.push_back(sparseIndexSend);
        sendIndicesForCommAfterFtoCPositions.push_back(posInSendIndices);
    }
}

void GridGenerator::findIfSparseIndexIsInSendIndicesAndAddToCommVectors(
    int sparseIndex, int *sendIndices, uint numberOfSendIndices, std::vector<int> &sendIndicesAfterFtoC,
    std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int sparseIndexSend;
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        sparseIndexSend = sendIndices[posInSendIndices];
        if (sparseIndex == sparseIndexSend) {
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndex,
                                                 sendIndicesForCommAfterFtoCPositions, posInSendIndices);
            break;
        }
    }
}

void GridGenerator::findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices, 
                                                      int *sendOrReceiveIndices, std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                                      std::vector<int> &sendOrIndicesOther)
{
    int sparseIndexSend;
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendOrRecvIndices; posInSendIndices++) {
        sparseIndexSend = sendOrReceiveIndices[posInSendIndices];
        if (std::find(sendOrReceiveIndicesAfterFtoC.begin(), sendOrReceiveIndicesAfterFtoC.end(), sparseIndexSend) ==
            sendOrReceiveIndicesAfterFtoC.end())
            sendOrIndicesOther.push_back(sparseIndexSend);
    }
}

void GridGenerator::reorderRecvIndicesForCommAfterFtoCX(int direction, int level, int j,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    int *recvIndices                    = para->getParH(level)->recvProcessNeighborX[j].index;
    int &numberOfRecvNeighborsAfterFtoC = para->getParH(level)->recvProcessNeighborsAfterFtoCX[j].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderRecvIndicesForCommAfterFtoCY(int direction, int level, int j,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    int *recvIndices                    = para->getParH(level)->recvProcessNeighborY[j].index;
    int &numberOfRecvNeighborsAfterFtoC = para->getParH(level)->recvProcessNeighborsAfterFtoCY[j].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderRecvIndicesForCommAfterFtoCZ(int direction, int level, int j,
                                                        std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    int *recvIndices                    = para->getParH(level)->recvProcessNeighborZ[j].index;
    int &numberOfRecvNeighborsAfterFtoC = para->getParH(level)->recvProcessNeighborsAfterFtoCZ[j].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNeighborsAfterFtoC, direction, level, j,
                                       sendIndicesForCommAfterFtoCPositions);
}

void GridGenerator::reorderRecvIndicesForCommAfterFtoC(int *recvIndices,
                                                       int &numberOfRecvNeighborsAfterFtoC, int direction, int level,
                                                       int j,
                                                       std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "reorder receive indices for communication after fine to coarse: level: " << level
                  << " direction: " << direction;
    if (sendIndicesForCommAfterFtoCPositions.size() == 0)
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderRecvIndicesForCommAfterFtoC(): sendIndicesForCommAfterFtoCPositions is empty."
                      << "\n";

    uint numberOfRecvIndices = builder->getNumberOfReceiveIndices(direction, level);
    std::vector<int> recvIndicesAfterFtoC;
    std::vector<int> recvIndicesOther;

    // find recvIndices for Communication after fine to coarse
    for (uint vectorPos : sendIndicesForCommAfterFtoCPositions)
        recvIndicesAfterFtoC.push_back(recvIndices[vectorPos]);

    findIndicesNotInCommAfterFtoC(numberOfRecvIndices, recvIndices, recvIndicesAfterFtoC, recvIndicesOther);

    numberOfRecvNeighborsAfterFtoC = (int)recvIndicesAfterFtoC.size();

    // copy new vectors back to sendIndices array
    for (int i = 0; i < numberOfRecvNeighborsAfterFtoC; i++)
        recvIndices[i] = recvIndicesAfterFtoC[i];
    for (uint i = 0; i < (uint)recvIndicesOther.size(); i++)
        recvIndices[i + numberOfRecvNeighborsAfterFtoC] = recvIndicesOther[i];

    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "... numberOfRecvNeighborsAfterFtoC: " << numberOfRecvNeighborsAfterFtoC << "\n";

    if (numberOfRecvNeighborsAfterFtoC + recvIndicesOther.size() != numberOfRecvIndices) {
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderRecvIndicesForCommAfterFtoC(): incorrect number of nodes"
                      << "\n";
        std::cout << "numberOfRecvNeighborsAfterFtoC = " << numberOfRecvNeighborsAfterFtoC
                  << ", recvIndicesOther.size() = " << recvIndicesOther.size()
                  << ", numberOfRecvIndices = " << numberOfRecvIndices << std::endl;
    }
}

void GridGenerator::allocArrays_BoundaryQs()
{
	std::cout << "------read BoundaryQs-------" << std::endl;


    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        int numberOfPressureValues = (int)builder->getPressureSize(i);
        if (numberOfPressureValues > 0)
        {
            std::cout << "size Pressure:  " << i << " : " << numberOfPressureValues << std::endl;
            //cout << "Groesse Pressure:  " << i << " : " << temp1 << "MyID: " << para->getMyID() << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->QPress.q27[0];
            unsigned int sizeQ = para->getParH(i)->QPress.kQ;
            QforBoundaryConditions Q;
            Q.q27[dirE] = &QQ[dirE   *sizeQ];
            Q.q27[dirW] = &QQ[dirW   *sizeQ];
            Q.q27[dirN] = &QQ[dirN   *sizeQ];
            Q.q27[dirS] = &QQ[dirS   *sizeQ];
            Q.q27[dirT] = &QQ[dirT   *sizeQ];
            Q.q27[dirB] = &QQ[dirB   *sizeQ];
            Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
            Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
            Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
            Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
            Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
            Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
            Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
            Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
            Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
            Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
            Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
            Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
            Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
            Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
            Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
            Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
            Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
            Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
            Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
            Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
            Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
            
            builder->getPressureQs(Q.q27, i);


            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            //cout << "vor advec diff" << endl;
            if (para->getDiffOn() == true) {
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor setzen von kTemp" << endl;
                para->getParH(i)->TempPress.kTemp = numberOfPressureValues;
                para->getParD(i)->TempPress.kTemp = numberOfPressureValues;
                std::cout << "Groesse TempPress.kTemp = " << para->getParH(i)->TempPress.kTemp << std::endl;
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempPressBC(i);
                //cout << "nach alloc" << endl;
                //////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfPressureValues; m++)
                {
                    para->getParH(i)->TempPress.temp[m] = para->getTemperatureInit();
                    para->getParH(i)->TempPress.velo[m] = (real)0.0;
                    para->getParH(i)->TempPress.k[m] = para->getParH(i)->QPress.k[m];
                }
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor copy" << endl;
                cudaMemoryManager->cudaCopyTempPressBCHD(i);
                //cout << "nach copy" << endl;
                //////////////////////////////////////////////////////////////////////////
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaCopyPress(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }//ende if
    }//ende oberste for schleife



    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const auto numberOfVelocityNodes = int(builder->getVelocitySize(i));
        if (numberOfVelocityNodes > 0)
        {
            std::cout << "size velocity level " << i << " : " << numberOfVelocityNodes << std::endl;
            //cout << "Groesse velocity level:  " << i << " : " << temp3 << "MyID: " << para->getMyID() << std::endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->Qinflow.q27[0];
            unsigned int sizeQ = para->getParH(i)->Qinflow.kQ;
            QforBoundaryConditions Q;
            Q.q27[dirE] = &QQ[dirE   *sizeQ];
            Q.q27[dirW] = &QQ[dirW   *sizeQ];
            Q.q27[dirN] = &QQ[dirN   *sizeQ];
            Q.q27[dirS] = &QQ[dirS   *sizeQ];
            Q.q27[dirT] = &QQ[dirT   *sizeQ];
            Q.q27[dirB] = &QQ[dirB   *sizeQ];
            Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
            Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
            Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
            Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
            Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
            Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
            Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
            Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
            Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
            Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
            Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
            Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
            Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
            Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
            Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
            Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
            Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
            Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
            Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
            Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
            Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];

            builder->getVelocityQs(Q.q27, i);

            if (para->getDiffOn()) {
                //////////////////////////////////////////////////////////////////////////
                para->getParH(i)->TempVel.kTemp = numberOfVelocityNodes;
                para->getParD(i)->TempVel.kTemp = numberOfVelocityNodes;
                std::cout << "Groesse TempVel.kTemp = " << para->getParH(i)->TempPress.kTemp << std::endl;
                std::cout << "getTemperatureInit = " << para->getTemperatureInit() << std::endl;
                std::cout << "getTemperatureBC = " << para->getTemperatureBC() << std::endl;
                //////////////////////////////////////////////////////////////////////////
                cudaMemoryManager->cudaAllocTempVeloBC(i);
                //cout << "nach alloc " << std::endl;
                //////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfVelocityNodes; m++)
                {
                    para->getParH(i)->TempVel.temp[m] = para->getTemperatureInit();
                    para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
                    para->getParH(i)->TempVel.velo[m] = para->getVelocity();
                    para->getParH(i)->TempVel.k[m] = para->getParH(i)->Qinflow.k[m];
                }
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor copy " << std::endl;
                cudaMemoryManager->cudaCopyTempVeloBCHD(i);
                //cout << "nach copy " << std::endl;
                //////////////////////////////////////////////////////////////////////////
            }
            cudaMemoryManager->cudaCopyVeloBC(i);
        }
    }


    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const int numberOfGeometryNodes = builder->getGeometrySize(i);
        std::cout << "size of GeomBoundaryQs, Level " << i << " : " << numberOfGeometryNodes << std::endl;

        para->getParH(i)->QGeom.kQ = numberOfGeometryNodes;
        para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
        if (numberOfGeometryNodes > 0)
        {
            //cout << "Groesse der Daten GeomBoundaryQs, Level:  " << i << " : " << numberOfGeometryNodes << "MyID: " << para->getMyID() << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //para->getParH(i)->QGeom.kQ = temp4;
            //para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cudaMemoryManager->cudaAllocGeomBC(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////////////////////
            //Indexarray
            builder->getGeometryIndices(para->getParH(i)->QGeom.k, i);
            //////////////////////////////////////////////////////////////////////////
            //preprocessing
            real* QQ = para->getParH(i)->QGeom.q27[0];
            unsigned int sizeQ = para->getParH(i)->QGeom.kQ;
            QforBoundaryConditions Q;
            Q.q27[dirE] = &QQ[dirE   *sizeQ];
            Q.q27[dirW] = &QQ[dirW   *sizeQ];
            Q.q27[dirN] = &QQ[dirN   *sizeQ];
            Q.q27[dirS] = &QQ[dirS   *sizeQ];
            Q.q27[dirT] = &QQ[dirT   *sizeQ];
            Q.q27[dirB] = &QQ[dirB   *sizeQ];
            Q.q27[dirNE] = &QQ[dirNE  *sizeQ];
            Q.q27[dirSW] = &QQ[dirSW  *sizeQ];
            Q.q27[dirSE] = &QQ[dirSE  *sizeQ];
            Q.q27[dirNW] = &QQ[dirNW  *sizeQ];
            Q.q27[dirTE] = &QQ[dirTE  *sizeQ];
            Q.q27[dirBW] = &QQ[dirBW  *sizeQ];
            Q.q27[dirBE] = &QQ[dirBE  *sizeQ];
            Q.q27[dirTW] = &QQ[dirTW  *sizeQ];
            Q.q27[dirTN] = &QQ[dirTN  *sizeQ];
            Q.q27[dirBS] = &QQ[dirBS  *sizeQ];
            Q.q27[dirBN] = &QQ[dirBN  *sizeQ];
            Q.q27[dirTS] = &QQ[dirTS  *sizeQ];
            Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
            Q.q27[dirTNE] = &QQ[dirTNE *sizeQ];
            Q.q27[dirTSW] = &QQ[dirTSW *sizeQ];
            Q.q27[dirTSE] = &QQ[dirTSE *sizeQ];
            Q.q27[dirTNW] = &QQ[dirTNW *sizeQ];
            Q.q27[dirBNE] = &QQ[dirBNE *sizeQ];
            Q.q27[dirBSW] = &QQ[dirBSW *sizeQ];
            Q.q27[dirBSE] = &QQ[dirBSE *sizeQ];
            Q.q27[dirBNW] = &QQ[dirBNW *sizeQ];
            //////////////////////////////////////////////////////////////////

            builder->getGeometryQs(Q.q27, i);
			//QDebugWriter::writeQValues(Q, para->getParH(i)->QGeom.k, para->getParH(i)->QGeom.kQ, "M:/TestGridGeneration/results/GeomGPU.dat");
            //////////////////////////////////////////////////////////////////
            for (int node_i = 0; node_i < numberOfGeometryNodes; node_i++)
            {
                Q.q27[dirZERO][node_i] = 0.0f;
            }
            //for(int test = 0; test < 3; test++)
            //{
            //	for (int tmp = 0; tmp < 27; tmp++)
            //	{
            //		cout <<"Kuhs: " << Q.q27[tmp][test]  << std::endl;
            //	}
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
                        para->getParH(i)->Temp.k[m] = para->getParH(i)->QGeom.k[m];
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


	std::cout << "-----finish BoundaryQs------" << std::endl;
}

void GridGenerator::allocArrays_OffsetScale()
{
    for (uint level = 0; level < builder->getNumberOfGridLevels() - 1; level++) 
    {
        const uint numberOfNodesPerLevelCF = builder->getNumberOfNodesCF(level);
        const uint numberOfNodesPerLevelFC = builder->getNumberOfNodesFC(level);

        std::cout << "number of nodes CF Level " << level << " : " << numberOfNodesPerLevelCF << std::endl;
        std::cout << "number of nodes FC level " << level << " : " << numberOfNodesPerLevelFC << std::endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size + memsize CF
        para->getParH(level)->K_CF = numberOfNodesPerLevelCF;
        para->getParD(level)->K_CF = para->getParH(level)->K_CF;
        para->getParH(level)->intCF.kCF = para->getParH(level)->K_CF;
        para->getParD(level)->intCF.kCF = para->getParH(level)->K_CF;
        para->getParH(level)->mem_size_kCF = sizeof(uint)* para->getParH(level)->K_CF;
        para->getParD(level)->mem_size_kCF = sizeof(uint)* para->getParD(level)->K_CF;
        para->getParH(level)->mem_size_kCF_off = sizeof(real)* para->getParH(level)->K_CF;
        para->getParD(level)->mem_size_kCF_off = sizeof(real)* para->getParD(level)->K_CF;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size + memsize FC
        para->getParH(level)->K_FC = numberOfNodesPerLevelFC;
        para->getParD(level)->K_FC = para->getParH(level)->K_FC;
        para->getParH(level)->intFC.kFC = para->getParH(level)->K_FC;
        para->getParD(level)->intFC.kFC = para->getParH(level)->K_FC;
        para->getParH(level)->mem_size_kFC = sizeof(uint)* para->getParH(level)->K_FC;
        para->getParD(level)->mem_size_kFC = sizeof(uint)* para->getParD(level)->K_FC;
        para->getParH(level)->mem_size_kFC_off = sizeof(real)* para->getParH(level)->K_FC;
        para->getParD(level)->mem_size_kFC_off = sizeof(real)* para->getParD(level)->K_FC;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //alloc
		cudaMemoryManager->cudaAllocInterfaceCF(level);
		cudaMemoryManager->cudaAllocInterfaceFC(level);
		cudaMemoryManager->cudaAllocInterfaceOffCF(level);
		cudaMemoryManager->cudaAllocInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //init
        builder->getOffsetCF(para->getParH(level)->offCF.xOffCF, para->getParH(level)->offCF.yOffCF, para->getParH(level)->offCF.zOffCF, level);
        builder->getOffsetFC(para->getParH(level)->offFC.xOffFC, para->getParH(level)->offFC.yOffFC, para->getParH(level)->offFC.zOffFC, level);
        builder->getGridInterfaceIndices(para->getParH(level)->intCF.ICellCFC, para->getParH(level)->intCF.ICellCFF, para->getParH(level)->intFC.ICellFCC, para->getParH(level)->intFC.ICellFCF, level);
        
        if (para->getUseStreams() || para->getNumprocs() > 1) {
            // split fine-to-coarse indices into border and bulk
            para->getParH(level)->intFCBorder.ICellFCC = para->getParH(level)->intFC.ICellFCC; 
            para->getParH(level)->intFCBorder.ICellFCF = para->getParH(level)->intFC.ICellFCF; 

            builder->getGridInterfaceIndicesFCBorderBulk(
                para->getParH(level)->intFCBorder.ICellFCC, para->getParH(level)->intFCBulk.ICellFCC,
                para->getParH(level)->intFCBorder.ICellFCF, para->getParH(level)->intFCBulk.ICellFCF,
                para->getParH(level)->intFCBorder.kFC, para->getParH(level)->intFCBulk.kFC, level);
            
            para->getParD(level)->intFCBorder.kFC = para->getParH(level)->intFCBorder.kFC;
            para->getParD(level)->intFCBulk.kFC = para->getParH(level)->intFCBulk.kFC;
            para->getParD(level)->intFCBorder.ICellFCC = para->getParD(level)->intFC.ICellFCC;
            para->getParD(level)->intFCBulk.ICellFCC = para->getParD(level)->intFCBorder.ICellFCC + para->getParD(level)->intFCBorder.kFC;
            para->getParD(level)->intFCBorder.ICellFCF = para->getParD(level)->intFC.ICellFCF;
            para->getParD(level)->intFCBulk.ICellFCF = para->getParD(level)->intFCBorder.ICellFCF + para->getParD(level)->intFCBorder.kFC;
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

    for (uint index = 0; index < para->getParH(level)->size_Mat_SP; index++)
        oss << verifyNeighborIndex(level, index, invalidNodes, stopperNodes, wrongNeighbors);


    oss << "invalid nodes found: " << invalidNodes << "\n";
    oss << "wrong neighbors found: " << wrongNeighbors << "\n";
    oss << "stopper nodes found : " << stopperNodes << "\n";
    oss << "---------report end---------\n";
    return oss.str();
}

std::string GridGenerator::verifyNeighborIndex(int level, int index , int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const
{
    std::ostringstream oss;

    const int geo = para->getParH(level)->geoSP[index];
    if (geo == 16)
    {
        stopperNodes++;
        return "";
    }

    real x = para->getParH(level)->coordX_SP[index];
    real y = para->getParH(level)->coordY_SP[index];
    real z = para->getParH(level)->coordZ_SP[index];

    real delta = para->getParH(level)->coordX_SP[2] - para->getParH(level)->coordX_SP[1];

    //std::cout << para->getParH(level)->coordX_SP[1] << ", " << para->getParH(level)->coordY_SP[1] << ", " << para->getParH(level)->coordZ_SP[1] << std::endl;
    //std::cout << para->getParH(level)->coordX_SP[para->getParH(level)->size_Mat_SP - 1] << ", " << para->getParH(level)->coordY_SP[para->getParH(level)->size_Mat_SP - 1] << ", " << para->getParH(level)->coordZ_SP[para->getParH(level)->size_Mat_SP - 1] << std::endl;
    
    real maxX = para->getParH(level)->coordX_SP[para->getParH(level)->size_Mat_SP - 1] - delta;
    real maxY = para->getParH(level)->coordY_SP[para->getParH(level)->size_Mat_SP - 1] - delta;
    real maxZ = para->getParH(level)->coordZ_SP[para->getParH(level)->size_Mat_SP - 1] - delta;
    real realNeighborX = vf::Math::lessEqual(x + delta, maxX) ? x + delta : para->getParH(level)->coordX_SP[1];
    real realNeighborY = vf::Math::lessEqual(y + delta, maxY) ? y + delta : para->getParH(level)->coordY_SP[1];
    real realNeighborZ = vf::Math::lessEqual(z + delta, maxZ) ? z + delta : para->getParH(level)->coordZ_SP[1];

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborX_SP[index], realNeighborX, y, z, "X");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborY_SP[index], x, realNeighborY, z, "Y");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ_SP[index], x, y, realNeighborZ, "Z");

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborY_SP[this->para->getParH(level)->neighborX_SP[index]], realNeighborX, realNeighborY, z, "XY");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ_SP[this->para->getParH(level)->neighborX_SP[index]], realNeighborX, y, realNeighborZ, "XZ");
    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ_SP[this->para->getParH(level)->neighborY_SP[index]], x, realNeighborY, realNeighborZ, "YZ");

    oss << checkNeighbor(level, x, y, z, index, wrongNeighbors, this->para->getParH(level)->neighborZ_SP[this->para->getParH(level)->neighborY_SP[this->para->getParH(level)->neighborX_SP[index]]], realNeighborX, realNeighborY, realNeighborZ, "XYZ");

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

    real neighborCoordX = para->getParH(level)->coordX_SP[neighborIndex];
    real neighborCoordY = para->getParH(level)->coordY_SP[neighborIndex];
    real neighborCoordZ = para->getParH(level)->coordZ_SP[neighborIndex];

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
