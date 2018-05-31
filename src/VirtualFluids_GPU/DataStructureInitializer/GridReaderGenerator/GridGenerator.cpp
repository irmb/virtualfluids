#include "GridGenerator.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GPU/CudaMemoryManager.h>

#include <sstream>
#include <iostream>
#include "utilities/math/Math.h"
#include "LBM/LB.h"


std::shared_ptr<GridProvider> GridGenerator::make(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<GridProvider>(new GridGenerator(builder, para));
}

GridGenerator::GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para)
{
	this->builder = builder;
    this->para = para;
    this->cudaMemoryManager = CudaMemoryManager::make(para);
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

		builder->getNodeValues(
			para->getParH(level)->coordX_SP,
			para->getParH(level)->coordY_SP,
			para->getParH(level)->coordZ_SP,
			para->getParH(level)->neighborX_SP,
			para->getParH(level)->neighborY_SP,
			para->getParH(level)->neighborZ_SP,
			para->getParH(level)->geoSP,
			level);

		setInitalNodeValues(numberOfNodesPerLevel, level);

        cudaMemoryManager->cudaCopySP(level);
        cudaMemoryManager->cudaCopyCoord(level);

        //std::cout << verifyNeighborIndices(level);
	}
	std::cout << "Number of Nodes: " << numberOfNodesGlobal << std::endl;
	std::cout << "-----finish Coord, Neighbor, Geo------" << std::endl;
}

void GridGenerator::allocArrays_BoundaryValues()
{
	std::cout << "------read BoundaryValues------" << std::endl;

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfPressureValues = int(builder->getPressureSize(level));

        cout << "size pressure level " << level << " : " << numberOfPressureValues << endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(level)->QPress.kQ = numberOfPressureValues;
        para->getParD(level)->QPress.kQ = numberOfPressureValues;
        para->getParH(level)->kPressQread = numberOfPressureValues * para->getD3Qxx();
        para->getParD(level)->kPressQread = numberOfPressureValues * para->getD3Qxx();
        if (numberOfPressureValues > 1)
        {
            para->cudaAllocPress(level);
            builder->getPressureValues(para->getParH(level)->QPress.RhoBC, para->getParH(level)->QPress.k, para->getParH(level)->QPress.kN, level);
            para->cudaCopyPress(level);
        }
    }
    

    for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
        const auto numberOfVelocityValues = int(builder->getVelocitySize(level));
        cout << "size velocity level " << level << " : " << numberOfVelocityValues << endl;
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
            para->cudaAllocVeloBC(level);
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

            para->cudaCopyVeloBC(level);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn()==true){
            	//////////////////////////////////////////////////////////////////////////
            	para->getParH(level)->TempVel.kTemp = numberOfVelocityValues;
            	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
            	cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
            	cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
            	//////////////////////////////////////////////////////////////////////////
            	para->cudaAllocTempVeloBC(level);
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
            	para->cudaCopyTempVeloBCHD(level);
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
            cout << "size geometry values, Level " << i << " : " << numberOfGeometryValues << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->getParH(i)->QGeom.kQ = numberOfGeometryValues;
            para->getParD(i)->QGeom.kQ = numberOfGeometryValues;
            if (numberOfGeometryValues > 0)
            {

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                para->cudaAllocGeomValuesBC(i);
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
                ////Täst
                //for (int m = 0; m < temp4; m++)
                //{
                //	para->getParH(i)->QGeom.Vx[m] = para->getVelocity();//0.035f;
                //	para->getParH(i)->QGeom.Vy[m] = 0.0f;//para->getVelocity();//0.0f;
                //	para->getParH(i)->QGeom.Vz[m] = 0.0f;
                //}
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                para->cudaCopyGeomValuesBC(i);
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// advection - diffusion stuff
                //if (para->getDiffOn()==true){
                //	//////////////////////////////////////////////////////////////////////////
                //	para->getParH(i)->Temp.kTemp = temp4;
                //	cout << "Groesse kTemp = " << para->getParH(i)->Temp.kTemp << endl;
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


}


void GridGenerator::allocArrays_BoundaryQs()
{
	std::cout << "------read BoundaryQs-------" << std::endl;


    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        int numberOfPressureValues = (int)builder->getPressureSize(i);
        if (numberOfPressureValues > 0)
        {
            cout << "size Pressure:  " << i << " : " << numberOfPressureValues << endl;
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
                cout << "Groesse TempPress.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
                //////////////////////////////////////////////////////////////////////////
                para->cudaAllocTempPressBC(i);
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
                para->cudaCopyTempPressBCHD(i);
                //cout << "nach copy" << endl;
                //////////////////////////////////////////////////////////////////////////
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->cudaCopyPress(i);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }//ende if
    }//ende oberste for schleife



    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const auto numberOfVelocityNodes = int(builder->getVelocitySize(i));
        if (numberOfVelocityNodes > 0)
        {
            cout << "size velocity level " << i << " : " << numberOfVelocityNodes << endl;
            //cout << "Groesse velocity level:  " << i << " : " << temp3 << "MyID: " << para->getMyID() << endl;
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
                cout << "Groesse TempVel.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
                cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
                cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
                //////////////////////////////////////////////////////////////////////////
                para->cudaAllocTempVeloBC(i);
                //cout << "nach alloc " << endl;
                //////////////////////////////////////////////////////////////////////////
                for (int m = 0; m < numberOfVelocityNodes; m++)
                {
                    para->getParH(i)->TempVel.temp[m] = para->getTemperatureInit();
                    para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
                    para->getParH(i)->TempVel.velo[m] = para->getVelocity();
                    para->getParH(i)->TempVel.k[m] = para->getParH(i)->Qinflow.k[m];
                }
                //////////////////////////////////////////////////////////////////////////
                //cout << "vor copy " << endl;
                para->cudaCopyTempVeloBCHD(i);
                //cout << "nach copy " << endl;
                //////////////////////////////////////////////////////////////////////////
            }
            para->cudaCopyVeloBC(i);
        }
    }


    for (uint i = 0; i < builder->getNumberOfGridLevels(); i++) {
        const int numberOfGeometryNodes = builder->getGeometrySize(i);
        cout << "size of GeomBoundaryQs, Level " << i << " : " << numberOfGeometryNodes << endl;

        para->getParH(i)->QGeom.kQ = numberOfGeometryNodes;
        para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
        if (numberOfGeometryNodes > 0)
        {
            //cout << "Groesse der Daten GeomBoundaryQs, Level:  " << i << " : " << numberOfGeometryNodes << "MyID: " << para->getMyID() << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //para->getParH(i)->QGeom.kQ = temp4;
            //para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->cudaAllocGeomBC(i);
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

            //////////////////////////////////////////////////////////////////
            for (int i = 0; i < numberOfGeometryNodes; i++)
            {
                Q.q27[dirZERO][i] = 0.0f;
            }
            //for(int test = 0; test < 3; test++)
            //{
            //	for (int tmp = 0; tmp < 27; tmp++)
            //	{
            //		cout <<"Kuhs: " << Q.q27[tmp][test]  << endl;				
            //	}
            //}

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // advection - diffusion stuff
            if (para->getDiffOn() == true) {
                    //////////////////////////////////////////////////////////////////////////
                    para->getParH(i)->Temp.kTemp = numberOfGeometryNodes;
                    para->getParD(i)->Temp.kTemp = numberOfGeometryNodes;
                    cout << "Groesse Temp.kTemp = " << para->getParH(i)->Temp.kTemp << endl;
                    //////////////////////////////////////////////////////////////////////////
                    para->cudaAllocTempNoSlipBC(i);
                    //////////////////////////////////////////////////////////////////////////
                    for (int m = 0; m < numberOfGeometryNodes; m++)
                    {
                        para->getParH(i)->Temp.temp[m] = para->getTemperatureInit();
                        para->getParH(i)->Temp.k[m] = para->getParH(i)->QGeom.k[m];
                    }
                    //////////////////////////////////////////////////////////////////////////
                    para->cudaCopyTempNoSlipBCHD(i);
                    //////////////////////////////////////////////////////////////////////////
                }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->cudaCopyGeomBC(i);
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

        cout << "number of nodes CF Level " << level << " : " << numberOfNodesPerLevelCF << endl;
        cout << "number of nodes FC level " << level << " : " << numberOfNodesPerLevelFC << endl;

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
        para->cudaAllocInterfaceCF(level);
        para->cudaAllocInterfaceFC(level);
        para->cudaAllocInterfaceOffCF(level);
        para->cudaAllocInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //init
        builder->setOffsetCF(para->getParH(level)->offCF.xOffCF, para->getParH(level)->offCF.yOffCF, para->getParH(level)->offCF.zOffCF, level);
        builder->setOffsetFC(para->getParH(level)->offFC.xOffFC, para->getParH(level)->offFC.yOffFC, para->getParH(level)->offFC.zOffFC, level);
        builder->getGridInterfaceIndices(para->getParH(level)->intCF.ICellCFC, para->getParH(level)->intCF.ICellCFF, para->getParH(level)->intFC.ICellFCC, para->getParH(level)->intFC.ICellFCF, level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //copy
        para->cudaCopyInterfaceCF(level);
        para->cudaCopyInterfaceFC(level);
        para->cudaCopyInterfaceOffCF(level);
        para->cudaCopyInterfaceOffFC(level);
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
