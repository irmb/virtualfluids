#include "GridGenerator.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GPU/CudaMemoryManager.h>

#include <sstream>
#include <iostream>
#include "utilities/math/Math.h"
#include "LBM/LB.h"
#include "Output/QDebugWriter.hpp"

#include "utilities/communication.h"



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

    if ((para->getNumprocs() > 1) /*&& (procNeighborsSendX.size() == procNeighborsRecvX.size())*/)
	{
		for (int direction = 0; direction < 6; direction++)
		{
            if( builder->getCommunicationProcess(direction) == INVALID_INDEX ) continue;

			for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
            {
                if( direction == CommunicationDirections::MX || direction == CommunicationDirections::PX )
                {
                    int j = (int)para->getParH(level)->sendProcessNeighborX.size();

		            para->getParH(level)->sendProcessNeighborX.emplace_back();
		            para->getParD(level)->sendProcessNeighborX.emplace_back();
		            para->getParH(level)->recvProcessNeighborX.emplace_back();
		            para->getParD(level)->recvProcessNeighborX.emplace_back();
		            if (para->getDiffOn()==true){
			            para->getParH(level)->sendProcessNeighborADX.emplace_back();
			            para->getParD(level)->sendProcessNeighborADX.emplace_back();
			            para->getParH(level)->recvProcessNeighborADX.emplace_back();
			            para->getParD(level)->recvProcessNeighborADX.emplace_back();
		            }

				    int tempSend = builder->getNumberOfSendIndices( direction, level );
				    int tempRecv = builder->getNumberOfReceiveIndices( direction, level );
				    if (tempSend > 0)
				    {
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //send
					    std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborX.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborX.back().numberOfNodes = tempSend;
					    para->getParD(level)->sendProcessNeighborX.back().numberOfNodes = tempSend;
					    para->getParH(level)->sendProcessNeighborX.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParD(level)->sendProcessNeighborX.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParH(level)->sendProcessNeighborX.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParD(level)->sendProcessNeighborX.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParH(level)->sendProcessNeighborX.back().memsizeFs = sizeof(real)     *tempSend;
					    para->getParD(level)->sendProcessNeighborX.back().memsizeFs = sizeof(real)     *tempSend;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //recv
					    std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborX.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborX.back().numberOfNodes = tempRecv;
					    para->getParD(level)->recvProcessNeighborX.back().numberOfNodes = tempRecv;
					    para->getParH(level)->recvProcessNeighborX.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParD(level)->recvProcessNeighborX.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParH(level)->recvProcessNeighborX.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParD(level)->recvProcessNeighborX.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParH(level)->recvProcessNeighborX.back().memsizeFs = sizeof(real)     *tempRecv;
					    para->getParD(level)->recvProcessNeighborX.back().memsizeFs = sizeof(real)     *tempRecv;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborX(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //init index arrays
                        builder->getSendIndices   (para->getParH(level)->sendProcessNeighborX[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborX[j].index, direction, level);
					    ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborXIndex(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
				    }
                }
                
                if( direction == CommunicationDirections::MY || direction == CommunicationDirections::PY )
                {
                    int j = (int)para->getParH(level)->sendProcessNeighborY.size();

		            para->getParH(level)->sendProcessNeighborY.emplace_back();
		            para->getParD(level)->sendProcessNeighborY.emplace_back();
		            para->getParH(level)->recvProcessNeighborY.emplace_back();
		            para->getParD(level)->recvProcessNeighborY.emplace_back();
		            if (para->getDiffOn()==true){
			            para->getParH(level)->sendProcessNeighborADY.emplace_back();
			            para->getParD(level)->sendProcessNeighborADY.emplace_back();
			            para->getParH(level)->recvProcessNeighborADY.emplace_back();
			            para->getParD(level)->recvProcessNeighborADY.emplace_back();
		            }

				    int tempSend = builder->getNumberOfSendIndices( direction, level );
				    int tempRecv = builder->getNumberOfReceiveIndices( direction, level );
				    if (tempSend > 0)
				    {
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //send
					    std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborY.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborY.back().numberOfNodes = tempSend;
					    para->getParD(level)->sendProcessNeighborY.back().numberOfNodes = tempSend;
					    para->getParH(level)->sendProcessNeighborY.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParD(level)->sendProcessNeighborY.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParH(level)->sendProcessNeighborY.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParD(level)->sendProcessNeighborY.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParH(level)->sendProcessNeighborY.back().memsizeFs = sizeof(real)     *tempSend;
					    para->getParD(level)->sendProcessNeighborY.back().memsizeFs = sizeof(real)     *tempSend;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //recv
					    std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborY.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborY.back().numberOfNodes = tempRecv;
					    para->getParD(level)->recvProcessNeighborY.back().numberOfNodes = tempRecv;
					    para->getParH(level)->recvProcessNeighborY.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParD(level)->recvProcessNeighborY.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParH(level)->recvProcessNeighborY.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParD(level)->recvProcessNeighborY.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParH(level)->recvProcessNeighborY.back().memsizeFs = sizeof(real)     *tempRecv;
					    para->getParD(level)->recvProcessNeighborY.back().memsizeFs = sizeof(real)     *tempRecv;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborY(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //init index arrays
                        builder->getSendIndices   (para->getParH(level)->sendProcessNeighborY[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborY[j].index, direction, level);
					    ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborYIndex(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
				    }
                }
                
                if( direction == CommunicationDirections::MZ || direction == CommunicationDirections::PZ )
                {
                    int j = (int)para->getParH(level)->sendProcessNeighborZ.size();

		            para->getParH(level)->sendProcessNeighborZ.emplace_back();
		            para->getParD(level)->sendProcessNeighborZ.emplace_back();
		            para->getParH(level)->recvProcessNeighborZ.emplace_back();
		            para->getParD(level)->recvProcessNeighborZ.emplace_back();
		            if (para->getDiffOn()==true){
			            para->getParH(level)->sendProcessNeighborADZ.emplace_back();
			            para->getParD(level)->sendProcessNeighborADZ.emplace_back();
			            para->getParH(level)->recvProcessNeighborADZ.emplace_back();
			            para->getParD(level)->recvProcessNeighborADZ.emplace_back();
		            }

				    int tempSend = builder->getNumberOfSendIndices( direction, level );
				    int tempRecv = builder->getNumberOfReceiveIndices( direction, level );
				    if (tempSend > 0)
				    {
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //send
					    std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborZ.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->sendProcessNeighborZ.back().numberOfNodes = tempSend;
					    para->getParD(level)->sendProcessNeighborZ.back().numberOfNodes = tempSend;
					    para->getParH(level)->sendProcessNeighborZ.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParD(level)->sendProcessNeighborZ.back().numberOfFs = para->getD3Qxx() * tempSend;
					    para->getParH(level)->sendProcessNeighborZ.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParD(level)->sendProcessNeighborZ.back().memsizeIndex = sizeof(unsigned int)*tempSend;
					    para->getParH(level)->sendProcessNeighborZ.back().memsizeFs = sizeof(real)     *tempSend;
					    para->getParD(level)->sendProcessNeighborZ.back().memsizeFs = sizeof(real)     *tempSend;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //recv
					    std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborZ.back().rankNeighbor = builder->getCommunicationProcess(direction);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    para->getParH(level)->recvProcessNeighborZ.back().numberOfNodes = tempRecv;
					    para->getParD(level)->recvProcessNeighborZ.back().numberOfNodes = tempRecv;
					    para->getParH(level)->recvProcessNeighborZ.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParD(level)->recvProcessNeighborZ.back().numberOfFs = para->getD3Qxx() * tempRecv;
					    para->getParH(level)->recvProcessNeighborZ.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParD(level)->recvProcessNeighborZ.back().memsizeIndex = sizeof(unsigned int)*tempRecv;
					    para->getParH(level)->recvProcessNeighborZ.back().memsizeFs = sizeof(real)     *tempRecv;
					    para->getParD(level)->recvProcessNeighborZ.back().memsizeFs = sizeof(real)     *tempRecv;
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //malloc on host and device
                        cudaMemoryManager->cudaAllocProcessNeighborZ(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
					    //init index arrays
                        builder->getSendIndices   (para->getParH(level)->sendProcessNeighborZ[j].index, direction, level);
                        builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborZ[j].index, direction, level);
					    ////////////////////////////////////////////////////////////////////////////////////////
                        cudaMemoryManager->cudaCopyProcessNeighborZIndex(level, j);
					    ////////////////////////////////////////////////////////////////////////////////////////
				    }
                }

			}
		}
	}


	// data exchange for F3 / G6
	if ((para->getNumprocs() > 1) && (para->getIsF3()) )
	{
		for (int direction = 0; direction < 6; direction++)
		{
			if (builder->getCommunicationProcess(direction) == INVALID_INDEX) continue;

			for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
			{
				if (direction == CommunicationDirections::MX || direction == CommunicationDirections::PX)
				{
                    int j = (int)para->getParH(level)->sendProcessNeighborF3X.size();

					para->getParH(level)->sendProcessNeighborF3X.emplace_back();
					para->getParD(level)->sendProcessNeighborF3X.emplace_back();
					para->getParH(level)->recvProcessNeighborF3X.emplace_back();
					para->getParD(level)->recvProcessNeighborF3X.emplace_back();

					int tempSend = builder->getNumberOfSendIndices(direction, level);
					int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3X.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3X.back().numberOfNodes = tempSend;
						para->getParD(level)->sendProcessNeighborF3X.back().numberOfNodes = tempSend;
						para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs = 6 * tempSend;
						para->getParD(level)->sendProcessNeighborF3X.back().numberOfGs = 6 * tempSend;
						para->getParH(level)->sendProcessNeighborF3X.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParD(level)->sendProcessNeighborF3X.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParH(level)->sendProcessNeighborF3X.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs;
						para->getParD(level)->sendProcessNeighborF3X.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3X.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3X.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3X.back().numberOfNodes = tempRecv;
						para->getParD(level)->recvProcessNeighborF3X.back().numberOfNodes = tempRecv;
						para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs = 6 * tempRecv;
						para->getParD(level)->recvProcessNeighborF3X.back().numberOfGs = 6 * tempRecv;
						para->getParH(level)->recvProcessNeighborF3X.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParD(level)->recvProcessNeighborF3X.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParH(level)->recvProcessNeighborF3X.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs;
						para->getParD(level)->recvProcessNeighborF3X.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3X.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						cudaMemoryManager->cudaAllocProcessNeighborF3X(level, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3X[j].index, direction, level);
						builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3X[j].index, direction, level);
						////////////////////////////////////////////////////////////////////////////////////////
						cudaMemoryManager->cudaCopyProcessNeighborF3XIndex(level, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}

				if (direction == CommunicationDirections::MY || direction == CommunicationDirections::PY)
				{
                    int j = (int)para->getParH(level)->sendProcessNeighborF3Y.size();

					para->getParH(level)->sendProcessNeighborF3Y.emplace_back();
					para->getParD(level)->sendProcessNeighborF3Y.emplace_back();
					para->getParH(level)->recvProcessNeighborF3Y.emplace_back();
					para->getParD(level)->recvProcessNeighborF3Y.emplace_back();

					int tempSend = builder->getNumberOfSendIndices(direction, level);
					int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3Y.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3Y.back().numberOfNodes = tempSend;
						para->getParD(level)->sendProcessNeighborF3Y.back().numberOfNodes = tempSend;
						para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs = 6 * tempSend;
						para->getParD(level)->sendProcessNeighborF3Y.back().numberOfGs = 6 * tempSend;
						para->getParH(level)->sendProcessNeighborF3Y.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParD(level)->sendProcessNeighborF3Y.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParH(level)->sendProcessNeighborF3Y.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs;
						para->getParD(level)->sendProcessNeighborF3Y.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3Y.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3Y.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3Y.back().numberOfNodes = tempRecv;
						para->getParD(level)->recvProcessNeighborF3Y.back().numberOfNodes = tempRecv;
						para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs = 6 * tempRecv;
						para->getParD(level)->recvProcessNeighborF3Y.back().numberOfGs = 6 * tempRecv;
						para->getParH(level)->recvProcessNeighborF3Y.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParD(level)->recvProcessNeighborF3Y.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParH(level)->recvProcessNeighborF3Y.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs;
						para->getParD(level)->recvProcessNeighborF3Y.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3Y.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						cudaMemoryManager->cudaAllocProcessNeighborF3Y(level, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3Y[j].index, direction, level);
						builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3Y[j].index, direction, level);
						////////////////////////////////////////////////////////////////////////////////////////
						cudaMemoryManager->cudaCopyProcessNeighborF3YIndex(level, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}

				if (direction == CommunicationDirections::MZ || direction == CommunicationDirections::PZ)
				{
                    int j = (int)para->getParH(level)->sendProcessNeighborF3Z.size();

					para->getParH(level)->sendProcessNeighborF3Z.emplace_back();
					para->getParD(level)->sendProcessNeighborF3Z.emplace_back();
					para->getParH(level)->recvProcessNeighborF3Z.emplace_back();
					para->getParD(level)->recvProcessNeighborF3Z.emplace_back();

					int tempSend = builder->getNumberOfSendIndices(direction, level);
					int tempRecv = builder->getNumberOfReceiveIndices(direction, level);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						std::cout << "size of Data for X send buffer, Level " << level << " : " << tempSend << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3Z.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->sendProcessNeighborF3Z.back().numberOfNodes = tempSend;
						para->getParD(level)->sendProcessNeighborF3Z.back().numberOfNodes = tempSend;
						para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs = 6 * tempSend;
						para->getParD(level)->sendProcessNeighborF3Z.back().numberOfGs = 6 * tempSend;
						para->getParH(level)->sendProcessNeighborF3Z.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParD(level)->sendProcessNeighborF3Z.back().memsizeIndex = sizeof(unsigned int) * tempSend;
						para->getParH(level)->sendProcessNeighborF3Z.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs;
						para->getParD(level)->sendProcessNeighborF3Z.back().memsizeGs = sizeof(real) * para->getParH(level)->sendProcessNeighborF3Z.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						std::cout << "size of Data for X receive buffer, Level " << level << " : " << tempRecv << std::endl;
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3Z.back().rankNeighbor = builder->getCommunicationProcess(direction);
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(level)->recvProcessNeighborF3Z.back().numberOfNodes = tempRecv;
						para->getParD(level)->recvProcessNeighborF3Z.back().numberOfNodes = tempRecv;
						para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs = 6 * tempRecv;
						para->getParD(level)->recvProcessNeighborF3Z.back().numberOfGs = 6 * tempRecv;
						para->getParH(level)->recvProcessNeighborF3Z.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParD(level)->recvProcessNeighborF3Z.back().memsizeIndex = sizeof(unsigned int) * tempRecv;
						para->getParH(level)->recvProcessNeighborF3Z.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs;
						para->getParD(level)->recvProcessNeighborF3Z.back().memsizeGs = sizeof(real) * para->getParH(level)->recvProcessNeighborF3Z.back().numberOfGs;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						cudaMemoryManager->cudaAllocProcessNeighborF3Z(level, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						builder->getSendIndices(para->getParH(level)->sendProcessNeighborF3Z[j].index, direction, level);
						builder->getReceiveIndices(para->getParH(level)->recvProcessNeighborF3Z[j].index, direction, level);
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
        
        if (para->getUseStreams()) {
            // split fine-to-coarse-coarse indices into border and bulk
            para->getParH(level)->intFCBorder.ICellFCC = para->getParH(level)->intFC.ICellFCC; 
            builder->getGridInterfaceIndicesFCCBorderBulk(para->getParH(level)->intFCBorder.ICellFCC, para->getParH(level)->intFCBorder.kFC, para->getParH(level)->intFCBulk.ICellFCC, para->getParH(level)->intFCBulk.kFC, level);
            
            para->getParD(level)->intFCBorder.kFC = para->getParH(level)->intFCBorder.kFC;
            para->getParD(level)->intFCBulk.kFC = para->getParH(level)->intFCBulk.kFC;
            para->getParD(level)->intFCBorder.ICellFCC = para->getParD(level)->intFC.ICellFCC;
            para->getParD(level)->intFCBulk.ICellFCC = para->getParD(level)->intFCBorder.ICellFCC + para->getParD(level)->intFCBorder.kFC;
        }
        std::cout << "sizeOld  " << para->getParH(level)->K_FC << std::endl;
        std::cout << "sizeNew  " << para->getParH(level)->intFCBorder.kFC + para->getParH(level)->intFCBulk.kFC
                  << " = border " << para->getParH(level)->intFCBorder.kFC << " + bulk "
                  << para->getParH(level)->intFCBulk.kFC << std::endl;
        std::cout << "first old  " << para->getParH(level)->intFC.ICellFCC[0] << std::endl;
        std::cout << "first new  " << para->getParH(level)->intFCBorder.ICellFCC[0]
                  << std::endl;
        //std::cout << "last border old  " << para->getParH(level)->intFC.ICellFCC[para->getParH(level)->intFCBorder.kFC - 1]
        //          << std::endl; //if (para->getParH(level)->intFCBorder.kFC > 0)
        //std::cout << "last border new  " << para->getParH(level)->intFCBorder.ICellFCC[para->getParH(level)->intFCBorder.kFC - 1]
        //          << std::endl;
        std::cout << "old pointer " << para->getParH(level)->intFC.ICellFCC << std::endl;
        std::cout << "border pointer (= old pointer) " << para->getParH(level)->intFCBorder.ICellFCC << std::endl;
        std::cout << "bulk pointer new " << para->getParH(level)->intFCBulk.ICellFCC << std::endl;
        std::cout << "first bulk old  "
                  << para->getParH(level)->intFC.ICellFCC[para->getParH(level)->intFCBorder.kFC] << std::endl;
        std::cout << "first bulk new  "
                  << para->getParH(level)->intFCBulk.ICellFCC[0] << std::endl;
        std::cout << "last bulk old  "
                  << para->getParH(level)->intFC.ICellFCC[para->getParH(level)->K_FC - 1] << std::endl;
        std::cout << "last bulk new  "
                  << para->getParH(level)->intFCBulk.ICellFCC[para->getParH(level)->intFCBulk.kFC - 1] << std::endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //copy
		cudaMemoryManager->cudaCopyInterfaceCF(level);
		cudaMemoryManager->cudaCopyInterfaceFC(level);
		cudaMemoryManager->cudaCopyInterfaceOffCF(level);
		cudaMemoryManager->cudaCopyInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        if (para->getUseStreams())
            cudaMemoryManager->cudaCheckInterfaceFCBulk(level);

        std::cout << "...Device " << std::endl;
        std::cout << "old pointer " << para->getParD(level)->intFC.ICellFCC << std::endl;
        std::cout << "border pointer (= old pointer) " << para->getParD(level)->intFCBorder.ICellFCC << std::endl;
        std::cout << "bulk pointer new " << para->getParD(level)->intFCBulk.ICellFCC << std::endl;
        std::cout << "sizeOld  " << para->getParD(level)->K_FC << std::endl;
        std::cout << "sizeNew  " << para->getParD(level)->intFCBorder.kFC + para->getParD(level)->intFCBulk.kFC
                  << " = border " << para->getParD(level)->intFCBorder.kFC << " + bulk "
                  << para->getParD(level)->intFCBulk.kFC << std::endl;
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
