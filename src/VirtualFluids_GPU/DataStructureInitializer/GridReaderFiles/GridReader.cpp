#include "GridReader.h"

#include <iostream>

#include "Parameter/Parameter.h"

#include "CoordNeighborGeoV.h"
#include "BoundaryQs.h"
#include "BoundaryValues.h"

#include <GPU/CudaMemoryManager.h>
#include "OffsetScale.h"

std::shared_ptr<GridProvider> GridReader::make(FileFormat format, std::shared_ptr<Parameter> para)
{
    return std::make_shared<GridReader>(format, para);
}

GridReader::GridReader(bool binaer, std::shared_ptr<Parameter> para)
{
    this->para = para;
    this->cudaMemoryManager = CudaMemoryManager::make(para);

	this->binaer=binaer;
	channelDirections.resize(6);
	channelBoundaryConditions.resize(6);
	BC_Values.resize(6);

	channelDirections[0] = "inlet";
	channelDirections[1] = "outlet";
	channelDirections[2] = "front";
	channelDirections[3] = "back";
	channelDirections[4] = "top";
	channelDirections[5] = "bottom";
}

GridReader::GridReader(FileFormat format, std::shared_ptr<Parameter> para)
{
    this->para = para;
    this->cudaMemoryManager = CudaMemoryManager::make(para);

    this->format = format;
    switch(format)
    {
    case FileFormat::ASCII:
        binaer = false; break;
    case FileFormat::BINARY:
        binaer = true; break;
    }

    channelDirections.resize(6);
    channelBoundaryConditions.resize(6);
    BC_Values.resize(6);

    channelDirections[0] = "inlet";
    channelDirections[1] = "outlet";
    channelDirections[2] = "front";
    channelDirections[3] = "back";
    channelDirections[4] = "top";
    channelDirections[5] = "bottom";
}

GridReader::~GridReader()
{

}


void rearrangeGeometry(Parameter* para, int lev)
{
    for (uint index = 0; index < para->getParH(lev)->size_Mat_SP; index++)
    {
        if (para->getParH(lev)->geoSP[index] == GEO_FLUID_OLD)
        {
            para->getParH(lev)->geoSP[index] = GEO_FLUID;
        }
    }
}

void GridReader::allocArrays_CoordNeighborGeo()
{
	std::cout << "-----Config Arrays Coord, Neighbor, Geo------" << std::endl;

    // Lenz: The first parameter in the CoordNeighborGeoV constructur is the file
	CoordNeighborGeoV coordX(para->getcoordX(), binaer, true);
	CoordNeighborGeoV coordY(para->getcoordY(), binaer, true);
	CoordNeighborGeoV coordZ(para->getcoordZ(), binaer, true);
	neighX = std::shared_ptr<CoordNeighborGeoV>(new CoordNeighborGeoV(para->getneighborX(), binaer, false));
	neighY = std::shared_ptr<CoordNeighborGeoV>(new CoordNeighborGeoV(para->getneighborY(), binaer, false));
	neighZ = std::shared_ptr<CoordNeighborGeoV>(new CoordNeighborGeoV(para->getneighborZ(), binaer, false));
	neighWSB = std::shared_ptr<CoordNeighborGeoV>(new CoordNeighborGeoV(para->getneighborWSB(), binaer, false));
	CoordNeighborGeoV geoV(para->getgeoVec(), binaer, false);

	int maxLevel = coordX.getLevel();
	std::cout << "Number of Level: " << maxLevel + 1 << std::endl;
	int numberOfNodesGlobal = 0;
	std::cout << "Number of Nodes: " << std::endl;

	for (int level = 0; level <= maxLevel; level++) 
	{		
		int numberOfNodesPerLevel = coordX.getSize(level) + 1;
		numberOfNodesGlobal += numberOfNodesPerLevel;
		std::cout << "Level " << level << " = " << numberOfNodesPerLevel << " Nodes" << std::endl;

		setNumberOfNodes(numberOfNodesPerLevel, level);

        cudaMemoryManager->cudaAllocCoord(level);
		cudaMemoryManager->cudaAllocSP(level);
		para->cudaAllocF3SP(level);
		cudaMemoryManager->cudaAllocNeighborWSB(level);

        if (para->getCalcMedian())
            para->cudaAllocMedianSP(level);
        if (para->getCalcParticle() || para->getUseWale())
            para->cudaAllocNeighborWSB(level);
        if (para->getUseWale())
            para->cudaAllocTurbulentViscosity(level);

		coordX.initalCoords(para->getParH(level)->coordX_SP, level);
		coordY.initalCoords(para->getParH(level)->coordY_SP, level);
		coordZ.initalCoords(para->getParH(level)->coordZ_SP, level);
		neighX->initalNeighbors(para->getParH(level)->neighborX_SP, level);
		neighY->initalNeighbors(para->getParH(level)->neighborY_SP, level);
		neighZ->initalNeighbors(para->getParH(level)->neighborZ_SP, level);
		neighWSB->initalNeighbors(para->getParH(level)->neighborWSB_SP, level);
		geoV.initalNeighbors(para->getParH(level)->geoSP, level);

        rearrangeGeometry(para.get(), level);
		setInitalNodeValues(numberOfNodesPerLevel, level);
        
        cudaMemoryManager->cudaCopyNeighborWSB(level);
        cudaMemoryManager->cudaCopySP(level);
        cudaMemoryManager->cudaCopyCoord(level);
	}
	std::cout << "Number of Nodes: " << numberOfNodesGlobal << std::endl;
	std::cout << "-----finish Coord, Neighbor, Geo------" <<std::endl;
}

void GridReader::allocArrays_BoundaryValues()
{
	std::cout << "------read BoundaryValues------" <<std::endl;

	this->makeReader(para);
	this->setChannelBoundaryCondition();
	int maxLevel = BC_Values[0]->getLevel();

	//initalValuesDomainDecompostion(maxLevel);

	vector<vector<vector<real> > > pressureV;
	pressureV.resize(maxLevel + 1);
	vector<vector<vector<real> > > velocityV;
	velocityV.resize(maxLevel + 1);
	vector<vector<vector<real> > > outflowV;
	outflowV.resize(maxLevel + 1);

	for (int i = 0; i <= maxLevel; i++) {
		pressureV[i].resize(2);
		velocityV[i].resize(3);
		outflowV[i].resize(2);
	}

	BoundaryValues *obj_geomV = new BoundaryValues(para->getgeomBoundaryBcValues(), para, "geo");

	////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	vector< BoundaryValues* > procNeighborsSendX, procNeighborsSendY, procNeighborsSendZ;
	vector< BoundaryValues* > procNeighborsRecvX, procNeighborsRecvY, procNeighborsRecvZ;
	vector< int >             neighborRankX, neighborRankY, neighborRankZ;

	procNeighborsSendX.resize(0);
	procNeighborsSendY.resize(0);
	procNeighborsSendZ.resize(0);
	procNeighborsRecvX.resize(0);
	procNeighborsRecvY.resize(0);
	procNeighborsRecvZ.resize(0);
	neighborRankX.resize(0);
	neighborRankY.resize(0);
	neighborRankZ.resize(0);

	if (para->getNumprocs() > 1)
	{
		for (int i = 0; i < para->getNumprocs(); i++)
		{
			BoundaryValues *pnXsend = new BoundaryValues(i, para, "send", "X");
			BoundaryValues *pnYsend = new BoundaryValues(i, para, "send", "Y");
			BoundaryValues *pnZsend = new BoundaryValues(i, para, "send", "Z");
			BoundaryValues *pnXrecv = new BoundaryValues(i, para, "recv", "X");
			BoundaryValues *pnYrecv = new BoundaryValues(i, para, "recv", "Y");
			BoundaryValues *pnZrecv = new BoundaryValues(i, para, "recv", "Z");
			if (para->getIsNeighborX())
			{
				procNeighborsSendX.push_back(pnXsend);
				procNeighborsRecvX.push_back(pnXrecv);
				neighborRankX.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankX: " << i << endl;
			}
			if (para->getIsNeighborY())
			{
				procNeighborsSendY.push_back(pnYsend);
				procNeighborsRecvY.push_back(pnYrecv);
				neighborRankY.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankY: " << i << endl;
			}
			if (para->getIsNeighborZ())
			{
				procNeighborsSendZ.push_back(pnZsend);
				procNeighborsRecvZ.push_back(pnZrecv);
				neighborRankZ.push_back(i);
				cout << "MyID: " << para->getMyID() << ", neighborRankZ: " << i << endl;
			}
		}
		cout << "MyID: " << para->getMyID() << ", size of neighborRankX: " << neighborRankX.size() << ", size of neighborRankY: " << neighborRankY.size() << ", size of neighborRankZ: " << neighborRankZ.size() << endl;
	}
	////////////////////////////////////////////////////////////////////////
	//3D domain decomposition convection diffusion
	vector< BoundaryValues* > procNeighborsSendADX, procNeighborsSendADY, procNeighborsSendADZ;
	vector< BoundaryValues* > procNeighborsRecvADX, procNeighborsRecvADY, procNeighborsRecvADZ;
	vector< int >             neighborRankADX, neighborRankADY, neighborRankADZ;

	procNeighborsSendADX.resize(0);
	procNeighborsSendADY.resize(0);
	procNeighborsSendADZ.resize(0);
	procNeighborsRecvADX.resize(0);
	procNeighborsRecvADY.resize(0);
	procNeighborsRecvADZ.resize(0);
	neighborRankADX.resize(0);
	neighborRankADY.resize(0);
	neighborRankADZ.resize(0);

	if (para->getDiffOn() == true && para->getNumprocs() > 1)
	{
		for (int i = 0; i < para->getNumprocs(); i++)
		{
			BoundaryValues *pnADXsend = new BoundaryValues(i, para, "send", "X");
			BoundaryValues *pnADYsend = new BoundaryValues(i, para, "send", "Y");
			BoundaryValues *pnADZsend = new BoundaryValues(i, para, "send", "Z");
			BoundaryValues *pnADXrecv = new BoundaryValues(i, para, "recv", "X");
			BoundaryValues *pnADYrecv = new BoundaryValues(i, para, "recv", "Y");
			BoundaryValues *pnADZrecv = new BoundaryValues(i, para, "recv", "Z");
			if (para->getIsNeighborX())
			{
				procNeighborsSendADX.push_back(pnADXsend);
				procNeighborsRecvADX.push_back(pnADXrecv);
				neighborRankADX.push_back(i);
			}
			if (para->getIsNeighborY())
			{
				procNeighborsSendADY.push_back(pnADYsend);
				procNeighborsRecvADY.push_back(pnADYrecv);
				neighborRankADY.push_back(i);
			}
			if (para->getIsNeighborZ())
			{
				procNeighborsSendADZ.push_back(pnADZsend);
				procNeighborsRecvADZ.push_back(pnADZrecv);
				neighborRankADZ.push_back(i);
			}
		}
	}


	for (int i = 0; i < channelBoundaryConditions.size(); i++) {
		if (this->channelBoundaryConditions[i] == "velocity") { BC_Values[i]->setBoundarys(velocityV); }
		else if (this->channelBoundaryConditions[i] == "pressure") { BC_Values[i]->setBoundarys(pressureV); }
		else if (this->channelBoundaryConditions[i] == "outflow") { BC_Values[i]->setBoundarys(outflowV); }
	}

	for (int i = 0; i <= maxLevel; i++) {
		int temp1 = (int)pressureV[i][0].size();
        cout << "size pressure level " << i << " : " << temp1 << endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(i)->QPress.kQ = temp1;
        para->getParD(i)->QPress.kQ = temp1;
        para->getParH(i)->kPressQread = temp1 * para->getD3Qxx();
        para->getParD(i)->kPressQread = temp1 * para->getD3Qxx();
		if (temp1 > 1)
		{
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocPress(i);
			///////////////////////////////
			////only for round of error test
			//para->cudaAllocTestRE(i, temp1);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int d = 0;
			int j = 0;
			int n = 0;

			for (vector<vector<vector<real> > >::iterator it = pressureV.begin(); it != pressureV.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							if (j == 0) para->getParH(i)->QPress.RhoBC[n] = *it3;
							if (j == 1) para->getParH(i)->QPress.kN[n] = (int)*it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp1; m++)
			{
				para->getParH(i)->QPress.RhoBC[m] = (para->getParH(i)->QPress.RhoBC[m] / para->getFactorPressBC() * (real)0.0);
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyPress(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////only for round of error test
			//para->cudaCopyTestREtoDevice(i,temp1);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//// advection - diffusion stuff
			////cout << "vor advec diff" << endl;
			//if (para->getDiffOn()==true){
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor setzen von kTemp" << endl;
			//	para->getParH(i)->TempPress.kTemp = temp1;
			//	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	para->cudaAllocTempPressBC(i);
			//	//cout << "nach alloc" << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	for (int m = 0; m < temp1; m++)
			//	{
			//		para->getParH(i)->TempPress.temp[m] = para->getTemperatureInit();
			//		para->getParH(i)->TempPress.velo[m] = (doubflo)0.0;
			//		para->getParH(i)->TempPress.k[m]    = para->getParH(i)->QPress.k[m];
			//	}
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor copy" << endl;
			//	para->cudaCopyTempPressBCHD(i);
			//	//cout << "nach copy" << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife




	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= maxLevel; i++) {
		int temp2 = (int)velocityV[i][0].size();
        cout << "size velocity level " << i << " : " << temp2 << endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int blocks = (temp2 / para->getParH(i)->numberofthreads) + 1;
        para->getParH(i)->Qinflow.kArray = blocks * para->getParH(i)->numberofthreads;
        para->getParD(i)->Qinflow.kArray = para->getParH(i)->Qinflow.kArray;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(i)->Qinflow.kQ = temp2;
        para->getParD(i)->Qinflow.kQ = temp2;
        para->getParH(i)->kInflowQ = temp2;
        para->getParD(i)->kInflowQ = temp2;
        para->getParH(i)->kInflowQread = temp2 * para->getD3Qxx();
        para->getParD(i)->kInflowQread = temp2 * para->getD3Qxx();
		if (temp2 > 1)
		{

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocVeloBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = velocityV.begin(); it != velocityV.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							if (j == 0) para->getParH(i)->Qinflow.Vx[n] = *it3;
							if (j == 1) para->getParH(i)->Qinflow.Vy[n] = *it3;
							if (j == 2) para->getParH(i)->Qinflow.Vz[n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
							 //cout << "n = " << n << endl;
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			//cout << "temp2 = " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp2; m++)
			{
				//para->getParH(i)->Qinflow.Vx[m] = para->getParH(i)->Qinflow.Vx[m] / para->getVelocityRatio();
				//para->getParH(i)->Qinflow.Vy[m] = para->getParH(i)->Qinflow.Vy[m] / para->getVelocityRatio();
				//para->getParH(i)->Qinflow.Vz[m] = para->getParH(i)->Qinflow.Vz[m] / para->getVelocityRatio();

                // Lenz: I do not know what this is, but it harms me ... comment out 
                // ==========================================================================
				//para->getParH(i)->Qinflow.Vx[m] = 0.0;//para->getVelocity();//0.035;
				//para->getParH(i)->Qinflow.Vy[m] = 0.0;//para->getVelocity();//0.0;
				//para->getParH(i)->Qinflow.Vz[m] = 0.0;
                // ==========================================================================

				//if (para->getParH(i)->Qinflow.Vz[m] > 0)
				//{
				//	cout << "velo Z = " << para->getParH(i)->Qinflow.Vz[m] << endl;
				//}
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyVeloBC(i);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//// advection - diffusion stuff
			//if (para->getDiffOn()==true){
			//	//////////////////////////////////////////////////////////////////////////
			//	para->getParH(i)->TempVel.kTemp = temp2;
			//	//cout << "Groesse kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
			//	cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
			//	cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	para->cudaAllocTempVeloBC(i);
			//	//cout << "nach alloc " << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//	for (int m = 0; m < temp2; m++)
			//	{
			//		para->getParH(i)->TempVel.temp[m]      = para->getTemperatureInit();
			//		para->getParH(i)->TempVel.tempPulse[m] = para->getTemperatureBC();
			//		para->getParH(i)->TempVel.velo[m]      = para->getVelocity();
			//		para->getParH(i)->TempVel.k[m]         = para->getParH(i)->Qinflow.k[m];
			//	}
			//	//////////////////////////////////////////////////////////////////////////
			//	//cout << "vor copy " << endl;
			//	para->cudaCopyTempVeloBCHD(i);
			//	//cout << "nach copy " << endl;
			//	//////////////////////////////////////////////////////////////////////////
			//}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife


	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= maxLevel; i++) {
		int temp = (int)outflowV[i][0].size();
        cout << "size outflow level " << i << " : " << temp << endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(i)->Qoutflow.kQ = temp;
        para->getParD(i)->Qoutflow.kQ = temp;
        para->getParH(i)->kOutflowQread = temp * para->getD3Qxx();
        para->getParD(i)->kOutflowQread = temp * para->getD3Qxx();
		if (temp > 1)
		{

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocOutflowBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int d = 0;
			int j = 0;
			int n = 0;

			for (vector<vector<vector<real> > >::iterator it = outflowV.begin(); it != outflowV.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							if (j == 0) para->getParH(i)->Qoutflow.RhoBC[n] = *it3;
							if (j == 1) para->getParH(i)->Qoutflow.kN[n] = (int)*it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int m = 0; m < temp; m++)
			{
				para->getParH(i)->Qoutflow.RhoBC[m] = (para->getParH(i)->Qoutflow.RhoBC[m] / para->getFactorPressBC()) * (real)0.0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyOutflowBC(i);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}//ende if
	}//ende oberste for schleife




	 //--------------------------------------------------------------------------//
	if (para->getIsGeometryValues()) {
		for (int i = 0; i <= maxLevel; i++) {
			int temp4 = obj_geomV->getSize(i);
            cout << "size obj_geomV, Level " << i << " : " << temp4 << endl;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            para->getParH(i)->QGeom.kQ = temp4;
            para->getParD(i)->QGeom.kQ = temp4;
			if (temp4 > 0)
			{

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaAllocGeomValuesBC(i);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_geomV->setValues(para->getParH(i)->QGeom.Vx, i, 0);
				obj_geomV->setValues(para->getParH(i)->QGeom.Vy, i, 1);
				obj_geomV->setValues(para->getParH(i)->QGeom.Vz, i, 2);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp4; m++)
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

	 //cout << "Test 2 " << endl;

	 ////--------------------------------------------------------------------------//
	 //if (para->getIsProp()){
	 //	BoundaryValues *obj_propV=new BoundaryValues(para->getpropellerValues(), para, "prop");
	 //	for (int i = 0; i <= maxLevel; i++) {
	 //		int temp4 = obj_propV->getSize(i);
	 //		if (temp4 > 0)
	 //		{
	 //			cout << "Groesse der Daten PropellerValues, Level " << i << " : " << temp4 << endl;
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //			para->getParH(i)->QPropeller.kQ = temp4;
	 //			para->getParD(i)->QPropeller.kQ = temp4;
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //			para->cudaAllocVeloPropeller(i);
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //			//Indexarray
	 //			obj_propV->initIndex(para->getParH(i)->QPropeller.k, i);
	 //			obj_propV->initArray(para->getParH(i)->QPropeller.Vx, i, 0);
	 //			obj_propV->initArray(para->getParH(i)->QPropeller.Vy, i, 1);
	 //			obj_propV->initArray(para->getParH(i)->QPropeller.Vz, i, 2);
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //			for (int m = 0; m < temp4; m++)
	 //			{
	 //				para->getParH(i)->QPropeller.Vx[m] = para->getParH(i)->QPropeller.Vx[m] / para->getVelocityRatio();
	 //				para->getParH(i)->QPropeller.Vy[m] = para->getParH(i)->QPropeller.Vy[m] / para->getVelocityRatio();
	 //				para->getParH(i)->QPropeller.Vz[m] = para->getParH(i)->QPropeller.Vz[m] / para->getVelocityRatio();
	 //			}
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //			para->cudaCopyVeloPropeller(i);
	 //			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //		}
	 //	}
	 //}//ende prop

	 //cout << "Test 3 " << endl;


	 //--------------------------------------------------------------------------//
	 //BoundaryValues *obj_cpTop=new BoundaryValues(para->getcpTop(), para, "cp");
	 //BoundaryValues *obj_cpBottom=new BoundaryValues(para->getcpBottom(), para, "cp");
	 //BoundaryValues *obj_cpBottom2=new BoundaryValues(para->getcpBottom2(), para, "cp");
	 //if (para->getIsCp()){
	 //	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //	//Top
	 //	for (int i = 0; i <= maxLevel; i++) {
	 //		int temp = obj_cpTop->getSize(i);
	 //		if (temp > 0)
	 //		{
	 //			cout << "Groesse der Daten CpTop, Level " << i << " : " << temp << endl;
	 //			////////////////////////////////////////////////////////////////////////////
	 //			para->getParH(i)->numberOfPointsCpTop = temp;
	 //			para->getParD(i)->numberOfPointsCpTop = temp;
	 //			////////////////////////////////////////////////////////////////////////////
	 //			para->cudaAllocCpTop(i);
	 //			////////////////////////////////////////////////////////////////////////////
	 //			//Indexarray
	 //			obj_cpTop->initIndex(para->getParH(i)->cpTopIndex, i);
	 //			////////////////////////////////////////////////////////////////////////////
	 //			for (int m = 0; m < temp; m++)
	 //			{
	 //				para->getParH(i)->cpPressTop[m] = 0.0;
	 //			}
	 //			////////////////////////////////////////////////////////////////////////////
	 //			para->cudaCopyCpTopInit(i);
	 //			////////////////////////////////////////////////////////////////////////////
	 //		}
	 //	}
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 ////Bottom
	 //for (int i = 0; i <= maxLevel; i++) {
	 //	int temp = obj_cpBottom->getSize(i);
	 //	if (temp > 0)
	 //	{
	 //		cout << "Groesse der Daten CpBottom, Level " << i << " : " << temp << endl;
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->getParH(i)->numberOfPointsCpBottom = temp;
	 //		para->getParD(i)->numberOfPointsCpBottom = temp;
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->cudaAllocCpBottom(i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //		//Indexarray
	 //		obj_cpBottom->initIndex(para->getParH(i)->cpBottomIndex, i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //		for (int m = 0; m < temp; m++)
	 //		{
	 //			para->getParH(i)->cpPressBottom[m] = 0.0;
	 //		}
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->cudaCopyCpBottomInit(i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //	}
	 //}
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 ////Bottom 2
	 //for (int i = 0; i <= maxLevel; i++) {
	 //	int temp = obj_cpBottom2->getSize(i);
	 //	if (temp > 0)
	 //	{
	 //		cout << "Groesse der Daten CpBottom2, Level " << i << " : " << temp << endl;
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->getParH(i)->numberOfPointsCpBottom2 = temp;
	 //		para->getParD(i)->numberOfPointsCpBottom2 = temp;
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->cudaAllocCpBottom2(i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //		//Indexarray
	 //		obj_cpBottom2->initIndex(para->getParH(i)->cpBottom2Index, i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //		for (int m = 0; m < temp; m++)
	 //		{
	 //			para->getParH(i)->cpPressBottom2[m] = 0.0;
	 //		}
	 //		////////////////////////////////////////////////////////////////////////////
	 //		para->cudaCopyCpBottom2Init(i);
	 //		////////////////////////////////////////////////////////////////////////////
	 //	}
	 //}
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //	delete obj_cpTop;
	 //	//delete obj_cpBottom;
	 //	//delete obj_cpBottom2;
	 //}//ende cp

	 //cout << "Test 4 " << endl;


	 //--------------------------------------------------------------------------//
	if (para->getConcFile()) {
		BoundaryValues *obj_Conc = new BoundaryValues(para->getConcentration());
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//concentration
		for (int i = 0; i <= maxLevel; i++) {
			int temp = obj_Conc->getSize(i);
            cout << "size Concentration, Level " << i << " : " << temp << endl;
            ////////////////////////////////////////////////////////////////////////////
            para->getParH(i)->numberOfPointsConc = temp;
            para->getParD(i)->numberOfPointsConc = temp;
			if (temp > 0)
			{

				////////////////////////////////////////////////////////////////////////////
				para->cudaAllocConcFile(i);
				////////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_Conc->initIndex(para->getParH(i)->concIndex, i);
				////////////////////////////////////////////////////////////////////////////
				para->cudaCopyConcFile(i);
				////////////////////////////////////////////////////////////////////////////
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		delete obj_Conc;
	}//end concentration

	 //cout << "Test 5 " << endl;



	 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //processor boundary (Martin Sch.) 
	 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //3D domain decomposition
	 // X
	if ((para->getNumprocs() > 1) && (procNeighborsSendX.size() == procNeighborsRecvX.size()))
	{
		for (int j = 0; j < procNeighborsSendX.size(); j++)
		{
			for (int i = 0; i <= maxLevel; i++) {
				int tempSend = procNeighborsSendX[j]->getSize(i);
				int tempRecv = procNeighborsRecvX[j]->getSize(i);
				cout << "maxlevel domains : " << maxLevel << ", tempsend : " << tempSend << endl;
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den X Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsSendX.size(), i, "send");
					para->getParH(i)->sendProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den X Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsRecvX.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborX(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendX[j]->initIndex(para->getParH(i)->sendProcessNeighborX[j].index, i);
					procNeighborsRecvX[j]->initIndex(para->getParH(i)->recvProcessNeighborX[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborXIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende X processor boundarys
	 //////////////////////////////////////////////////////////////////////////
	 // Y
	if ((para->getNumprocs() > 1) && (procNeighborsSendY.size() == procNeighborsRecvY.size()))
	{
		for (int j = 0; j < procNeighborsSendY.size(); j++)
		{
			for (int i = 0; i <= maxLevel; i++) {
				int tempSend = procNeighborsSendY[j]->getSize(i);
				int tempRecv = procNeighborsRecvY[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den Y Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsSendY.size(), i, "send");
					para->getParH(i)->sendProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den Y Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsRecvY.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborY(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendY[j]->initIndex(para->getParH(i)->sendProcessNeighborY[j].index, i);
					procNeighborsRecvY[j]->initIndex(para->getParH(i)->recvProcessNeighborY[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborYIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende Y processor boundarys
	 //////////////////////////////////////////////////////////////////////////
	 // Z
	if ((para->getNumprocs() > 1) && (procNeighborsSendZ.size() == procNeighborsRecvZ.size()))
	{
		for (int j = 0; j < procNeighborsSendZ.size(); j++)
		{
			for (int i = 0; i <= maxLevel; i++) {
				int tempSend = procNeighborsSendZ[j]->getSize(i);
				int tempRecv = procNeighborsRecvZ[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					cout << "Groesse der Daten für den Z Sendepuffer, Level " << i << " : " << tempSend << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsSendZ.size(), i, "send");
					para->getParH(i)->sendProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					cout << "Groesse der Daten für den Z Empfangspuffer, Level " << i << " : " << tempRecv << endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsRecvZ.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
					para->cudaAllocProcessNeighborZ(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendZ[j]->initIndex(para->getParH(i)->sendProcessNeighborZ[j].index, i);
					procNeighborsRecvZ[j]->initIndex(para->getParH(i)->recvProcessNeighborZ[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyProcessNeighborZIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende Z processor boundarys
	 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 //3D domain decomposition convection diffusion
	if (para->getDiffOn() == true) {
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// X
		if ((para->getNumprocs() > 1) && (procNeighborsSendADX.size() == procNeighborsRecvADX.size()))
		{
			for (int j = 0; j < procNeighborsSendADX.size(); j++)
			{
				for (int i = 0; i <= maxLevel; i++) {
					int tempSend = procNeighborsSendADX[j]->getSize(i);
					int tempRecv = procNeighborsRecvADX[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADX[j].rankNeighbor = neighborRankADX[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADX[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADX[j].memsizeFs = sizeof(real)     *tempSend;
						para->getParD(i)->sendProcessNeighborADX[j].memsizeFs = sizeof(real)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADX[j].rankNeighbor = neighborRankADX[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADX[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADX[j].memsizeFs = sizeof(real)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADX[j].memsizeFs = sizeof(real)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADX(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADX[j]->initIndex(para->getParH(i)->sendProcessNeighborADX[j].index, i);
						procNeighborsRecvADX[j]->initIndex(para->getParH(i)->recvProcessNeighborADX[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADXIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende X processor boundarys
		 //////////////////////////////////////////////////////////////////////////
		 // Y
		if ((para->getNumprocs() > 1) && (procNeighborsSendADY.size() == procNeighborsRecvADY.size()))
		{
			for (int j = 0; j < procNeighborsSendADY.size(); j++)
			{
				for (int i = 0; i <= maxLevel; i++) {
					int tempSend = procNeighborsSendADY[j]->getSize(i);
					int tempRecv = procNeighborsRecvADY[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADY[j].rankNeighbor = neighborRankADY[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADY[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADY[j].memsizeFs = sizeof(real)     *tempSend;
						para->getParD(i)->sendProcessNeighborADY[j].memsizeFs = sizeof(real)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADY[j].rankNeighbor = neighborRankADY[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADY[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADY[j].memsizeFs = sizeof(real)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADY[j].memsizeFs = sizeof(real)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADY(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADY[j]->initIndex(para->getParH(i)->sendProcessNeighborADY[j].index, i);
						procNeighborsRecvADY[j]->initIndex(para->getParH(i)->recvProcessNeighborADY[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADYIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende Y processor boundarys
		 //////////////////////////////////////////////////////////////////////////
		 // Z
		if ((para->getNumprocs() > 1) && (procNeighborsSendADZ.size() == procNeighborsRecvADZ.size()))
		{
			for (int j = 0; j < procNeighborsSendADZ.size(); j++)
			{
				for (int i = 0; i <= maxLevel; i++) {
					int tempSend = procNeighborsSendADZ[j]->getSize(i);
					int tempRecv = procNeighborsRecvADZ[j]->getSize(i);
					if (tempSend > 0)
					{
						////////////////////////////////////////////////////////////////////////////////////////
						//send
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADZ[j].rankNeighbor = neighborRankADZ[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->sendProcessNeighborADZ[j].numberOfNodes = tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].numberOfNodes = tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].numberOfFs = para->getD3Qxx() * tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
						para->getParH(i)->sendProcessNeighborADZ[j].memsizeFs = sizeof(real)     *tempSend;
						para->getParD(i)->sendProcessNeighborADZ[j].memsizeFs = sizeof(real)     *tempSend;
						////////////////////////////////////////////////////////////////////////////////////////
						//recv
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADZ[j].rankNeighbor = neighborRankADZ[j];
						////////////////////////////////////////////////////////////////////////////////////////
						para->getParH(i)->recvProcessNeighborADZ[j].numberOfNodes = tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].numberOfNodes = tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
						para->getParH(i)->recvProcessNeighborADZ[j].memsizeFs = sizeof(real)     *tempRecv;
						para->getParD(i)->recvProcessNeighborADZ[j].memsizeFs = sizeof(real)     *tempRecv;
						////////////////////////////////////////////////////////////////////////////////////////
						//malloc on host and device
						para->cudaAllocProcessNeighborADZ(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
						//init index arrays
						procNeighborsSendADZ[j]->initIndex(para->getParH(i)->sendProcessNeighborADZ[j].index, i);
						procNeighborsRecvADZ[j]->initIndex(para->getParH(i)->recvProcessNeighborADZ[j].index, i);
						////////////////////////////////////////////////////////////////////////////////////////
						para->cudaCopyProcessNeighborADZIndex(i, j);
						////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}//ende Z processor boundarys
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}













	//for (int level = 0; level < maxLevel; level++)
	//{
	//	setVelocitySizePerLevel(level, 0);
	//	setPressSizePerLevel(level, 0);
	//	setVelocitySizePerLevel(level, 0);
	//}

 //   for (int i = 0; i < channelBoundaryConditions.size(); i++)
 //   {
	//	if(this->channelBoundaryConditions[i] == "velocity")
	//		setVelocityValues(i);

	//	if (this->channelBoundaryConditions[i] == "pressure")
	//		setPressureValues(i);

	//	if (this->channelBoundaryConditions[i] == "outflow")
	//		setOutflowValues(i);
 //   }

}

void GridReader::allocArrays_OffsetScale()
{
    cout << "-----Config Arrays OffsetScale------" << endl;
    OffsetScale *obj_offCF = new OffsetScale(para->getscaleOffsetCF(), true);
    OffsetScale *obj_offFC = new OffsetScale(para->getscaleOffsetFC(), true);
    OffsetScale *obj_scaleCFC = new OffsetScale(para->getscaleCFC(), false);
    OffsetScale *obj_scaleCFF = new OffsetScale(para->getscaleCFF(), false);
    OffsetScale *obj_scaleFCC = new OffsetScale(para->getscaleFCC(), false);
    OffsetScale *obj_scaleFCF = new OffsetScale(para->getscaleFCF(), false);

    int maxLevel = obj_offCF->getLevel();

    int AnzahlKnotenGesCF = 0;
    int AnzahlKnotenGesFC = 0;

    for (int i = 0; i<maxLevel; i++) {
        unsigned int tempCF = obj_offCF->getSize(i);
        cout << "Groesse der Daten CF vom Level " << i << " : " << tempCF << endl;
        unsigned int tempFC = obj_offFC->getSize(i);
        cout << "Groesse der Daten FC vom Level " << i << " : " << tempFC << endl;

        AnzahlKnotenGesCF += tempCF;
        AnzahlKnotenGesFC += tempFC;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size + memsize CF
        para->getParH(i)->K_CF = tempCF;
        para->getParD(i)->K_CF = para->getParH(i)->K_CF;
        para->getParH(i)->intCF.kCF = para->getParH(i)->K_CF;
        para->getParD(i)->intCF.kCF = para->getParH(i)->K_CF;
        para->getParH(i)->mem_size_kCF = sizeof(unsigned int)* para->getParH(i)->K_CF;
        para->getParD(i)->mem_size_kCF = sizeof(unsigned int)* para->getParD(i)->K_CF;
        para->getParH(i)->mem_size_kCF_off = sizeof(real)* para->getParH(i)->K_CF;
        para->getParD(i)->mem_size_kCF_off = sizeof(real)* para->getParD(i)->K_CF;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //size + memsize FC
        para->getParH(i)->K_FC = tempFC;
        para->getParD(i)->K_FC = para->getParH(i)->K_FC;
        para->getParH(i)->intFC.kFC = para->getParH(i)->K_FC;
        para->getParD(i)->intFC.kFC = para->getParH(i)->K_FC;
        para->getParH(i)->mem_size_kFC = sizeof(unsigned int)* para->getParH(i)->K_FC;
        para->getParD(i)->mem_size_kFC = sizeof(unsigned int)* para->getParD(i)->K_FC;
        para->getParH(i)->mem_size_kFC_off = sizeof(real)* para->getParH(i)->K_FC;
        para->getParD(i)->mem_size_kFC_off = sizeof(real)* para->getParD(i)->K_FC;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //alloc
        para->cudaAllocInterfaceCF(i);
        para->cudaAllocInterfaceFC(i);
        para->cudaAllocInterfaceOffCF(i);
        para->cudaAllocInterfaceOffFC(i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //init
        obj_offCF->initArrayOffset(para->getParH(i)->offCF.xOffCF, para->getParH(i)->offCF.yOffCF, para->getParH(i)->offCF.zOffCF, i);
        obj_offFC->initArrayOffset(para->getParH(i)->offFC.xOffFC, para->getParH(i)->offFC.yOffFC, para->getParH(i)->offFC.zOffFC, i);
        obj_scaleCFC->initScale(para->getParH(i)->intCF.ICellCFC, i);
        obj_scaleCFF->initScale(para->getParH(i)->intCF.ICellCFF, i);
        obj_scaleFCC->initScale(para->getParH(i)->intFC.ICellFCC, i);
        obj_scaleFCF->initScale(para->getParH(i)->intFC.ICellFCF, i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //copy
        para->cudaCopyInterfaceCF(i);
        para->cudaCopyInterfaceFC(i);
        para->cudaCopyInterfaceOffCF(i);
        para->cudaCopyInterfaceOffFC(i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
    cout << "Gesamtanzahl Knoten CF = " << AnzahlKnotenGesCF << endl;
    cout << "Gesamtanzahl Knoten FC = " << AnzahlKnotenGesFC << endl;

    delete obj_offCF;
    delete obj_offFC;
    delete obj_scaleCFC;
    delete obj_scaleCFF;
    delete obj_scaleFCC;
    delete obj_scaleFCF;
    cout << "-----Ende OffsetScale------" << endl;
}


void GridReader::setPressureValues(int channelSide) const
{
	for (unsigned int level = 0; level <= BC_Values[channelSide]->getLevel(); level++)
	{
		int sizePerLevel = BC_Values[channelSide]->getSize(level);
        setPressSizePerLevel(level, sizePerLevel);

		if (sizePerLevel > 0)
		{
			std::cout << "size pressure level " << level << " : " << sizePerLevel << std::endl;

            cudaMemoryManager->cudaAllocPress(level);

			setPressRhoBC(sizePerLevel, level, channelSide);
            cudaMemoryManager->cudaCopyPress(level);
		}
	}
}

void GridReader::setPressRhoBC(int sizePerLevel, int level, int channelSide) const
{
	BC_Values[channelSide]->setPressValues(para->getParH(level)->QPress.RhoBC, para->getParH(level)->QPress.kN, level);
	for (int m = 0; m < sizePerLevel; m++)
		para->getParH(level)->QPress.RhoBC[m] = (para->getParH(level)->QPress.RhoBC[m] / para->getFactorPressBC());
}


void GridReader::setVelocityValues(int channelSide) const
{
	for (unsigned int level = 0; level <= BC_Values[channelSide]->getLevel(); level++)
	{
		int sizePerLevel = BC_Values[channelSide]->getSize(level);
        setVelocitySizePerLevel(level, sizePerLevel);

		if (sizePerLevel > 1)
		{
			std::cout << "size velocity level " << level << " : " << sizePerLevel << std::endl;

            cudaMemoryManager->cudaAllocVeloBC(level);

			setVelocity(level, sizePerLevel, channelSide);
            cudaMemoryManager->cudaCopyVeloBC(level);
		}
	}
}

void GridReader::setVelocity(int level, int sizePerLevel, int channelSide) const
{
	BC_Values[channelSide]->setVelocityValues(para->getParH(level)->Qinflow.Vx, para->getParH(level)->Qinflow.Vy, para->getParH(level)->Qinflow.Vz, level);

	for (int index = 0; index < sizePerLevel; index++)
	{
		//para->getParH(i)->Qinflow.Vx[m] = para->getParH(i)->Qinflow.Vx[m] / para->getVelocityRatio();
		//para->getParH(i)->Qinflow.Vy[m] = para->getParH(i)->Qinflow.Vy[m] / para->getVelocityRatio();
		//para->getParH(i)->Qinflow.Vz[m] = para->getParH(i)->Qinflow.Vz[m] / para->getVelocityRatio();
		para->getParH(level)->Qinflow.Vx[index] = 0.0;//para->getVelocity();//0.035;
		para->getParH(level)->Qinflow.Vy[index] = 0.0;//para->getVelocity();//0.0;
		para->getParH(level)->Qinflow.Vz[index] = 0.0;
	}
}


void GridReader::setOutflowValues(int channelSide) const
{
	for (unsigned int level = 0; level <= BC_Values[channelSide]->getLevel(); level++)
	{
		int sizePerLevel = BC_Values[channelSide]->getSize(level);
        setOutflowSizePerLevel(level, sizePerLevel);

		if (sizePerLevel > 1)
		{
			std::cout << "size outflow level " << level << " : " << sizePerLevel << std::endl;

            cudaMemoryManager->cudaAllocOutflowBC(level);

			setOutflow(level, sizePerLevel, channelSide);
            cudaMemoryManager->cudaCopyOutflowBC(level);

		}
	}
}

void GridReader::setOutflow(int level, int sizePerLevel, int channelSide) const
{
	BC_Values[channelSide]->setOutflowValues(para->getParH(level)->Qoutflow.RhoBC, para->getParH(level)->Qoutflow.kN, level);
	for (int index = 0; index < sizePerLevel; index++)
		para->getParH(level)->Qoutflow.RhoBC[index] = (para->getParH(level)->Qoutflow.RhoBC[index] / para->getFactorPressBC()) * (real)0.0;
}


void GridReader::initalValuesDomainDecompostion(int level)
{
	////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	std::vector< std::shared_ptr<BoundaryValues> > procNeighborsSendX, procNeighborsSendY, procNeighborsSendZ;
	std::vector< std::shared_ptr<BoundaryValues> > procNeighborsRecvX, procNeighborsRecvY, procNeighborsRecvZ;
	std::vector< int >             neighborRankX, neighborRankY, neighborRankZ;

	if (para->getNumprocs() > 1)
	{
		for (int process = 0; process < para->getNumprocs(); process++)
		{
			std::shared_ptr<BoundaryValues> pnXsend = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "send", "X"));
			std::shared_ptr<BoundaryValues> pnYsend = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "send", "Y"));
			std::shared_ptr<BoundaryValues> pnZsend = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "send", "Z"));
			std::shared_ptr<BoundaryValues> pnXrecv = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "recv", "X"));
			std::shared_ptr<BoundaryValues> pnYrecv = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "recv", "Y"));
			std::shared_ptr<BoundaryValues> pnZrecv = std::shared_ptr<BoundaryValues> (new BoundaryValues(process, para, "recv", "Z"));
			if (para->getIsNeighborX())
			{
				procNeighborsSendX.push_back(pnXsend);
				procNeighborsRecvX.push_back(pnXrecv);
				neighborRankX.push_back(process);
				std::cout << "MyID: " << para->getMyID() << ", neighborRankX: " << process << std::endl;
			}
			if (para->getIsNeighborY())
			{
				procNeighborsSendY.push_back(pnYsend);
				procNeighborsRecvY.push_back(pnYrecv);
				neighborRankY.push_back(process);
				std::cout << "MyID: " << para->getMyID() << ", neighborRankY: " << process << std::endl;
			}
			if (para->getIsNeighborZ())
			{
				procNeighborsSendZ.push_back(pnZsend);
				procNeighborsRecvZ.push_back(pnZrecv);
				neighborRankZ.push_back(process);
				std::cout << "MyID: " << para->getMyID() << ", neighborRankZ: " << process << std::endl;
			}
		}
		std::cout << "MyID: " << para->getMyID() << ", size of neighborRankX: " << neighborRankX.size() << ", size of neighborRankY: " << neighborRankY.size() << ", size of neighborRankZ: " << neighborRankZ.size() << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	// X
	if ((para->getNumprocs() > 1) && (procNeighborsSendX.size() == procNeighborsRecvX.size()))
	{
		for (int j = 0; j < procNeighborsSendX.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendX[j]->getSize(i);
				int tempRecv = procNeighborsRecvX[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					std::cout << "size of Data for X send buffer, Level " << i << " : " << tempSend << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsSendX.size(), i, "send");
					para->getParH(i)->sendProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborX[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborX[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					std::cout << "size of Data for X receive buffer, Level " << i << " : " << tempRecv << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsX((unsigned int)procNeighborsRecvX.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborX[j].rankNeighbor = neighborRankX[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborX[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborX[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
                    cudaMemoryManager->cudaAllocProcessNeighborX(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendX[j]->initIndex(para->getParH(i)->sendProcessNeighborX[j].index, i);
					procNeighborsRecvX[j]->initIndex(para->getParH(i)->recvProcessNeighborX[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
                    cudaMemoryManager->cudaCopyProcessNeighborXIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Y
	if ((para->getNumprocs() > 1) && (procNeighborsSendY.size() == procNeighborsRecvY.size()))
	{
		for (int j = 0; j < procNeighborsSendY.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendY[j]->getSize(i);
				int tempRecv = procNeighborsRecvY[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					std::cout << "size of Data for Y send buffer Level " << i << " : " << tempSend << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsSendY.size(), i, "send");
					para->getParH(i)->sendProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborY[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborY[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					std::cout << "size of Data for Y receive buffer, Level " << i << " : " << tempRecv << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsY((unsigned int)procNeighborsRecvY.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborY[j].rankNeighbor = neighborRankY[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborY[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborY[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
                    cudaMemoryManager->cudaAllocProcessNeighborY(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendY[j]->initIndex(para->getParH(i)->sendProcessNeighborY[j].index, i);
					procNeighborsRecvY[j]->initIndex(para->getParH(i)->recvProcessNeighborY[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
                    cudaMemoryManager->cudaCopyProcessNeighborYIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Z
	if ((para->getNumprocs() > 1) && (procNeighborsSendZ.size() == procNeighborsRecvZ.size()))
	{
		for (int j = 0; j < procNeighborsSendZ.size(); j++)
		{
			for (int i = 0; i <= level; i++) {
				int tempSend = procNeighborsSendZ[j]->getSize(i);
				int tempRecv = procNeighborsRecvZ[j]->getSize(i);
				if (tempSend > 0)
				{
					////////////////////////////////////////////////////////////////////////////////////////
					//send
					std::cout << "size of Data for Z send buffer, Level " << i << " : " << tempSend << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsSendZ.size(), i, "send");
					para->getParH(i)->sendProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfNodes = tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempSend;
					para->getParH(i)->sendProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempSend;
					para->getParD(i)->sendProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempSend;
					////////////////////////////////////////////////////////////////////////////////////////
					//recv
					std::cout << "size of Data for Z receive buffer, Level " << i << " : " << tempRecv << std::endl;
					////////////////////////////////////////////////////////////////////////////////////////
					para->setNumberOfProcessNeighborsZ((unsigned int)procNeighborsRecvZ.size(), i, "recv");
					para->getParH(i)->recvProcessNeighborZ[j].rankNeighbor = neighborRankZ[j];
					////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfNodes = tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].numberOfFs = para->getD3Qxx() * tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeIndex = sizeof(unsigned int)*tempRecv;
					para->getParH(i)->recvProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempRecv;
					para->getParD(i)->recvProcessNeighborZ[j].memsizeFs = sizeof(real)     *tempRecv;
					////////////////////////////////////////////////////////////////////////////////////////
					//malloc on host and device
                    cudaMemoryManager->cudaAllocProcessNeighborZ(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
					//init index arrays
					procNeighborsSendZ[j]->initIndex(para->getParH(i)->sendProcessNeighborZ[j].index, i);
					procNeighborsRecvZ[j]->initIndex(para->getParH(i)->recvProcessNeighborZ[j].index, i);
					////////////////////////////////////////////////////////////////////////////////////////
                    cudaMemoryManager->cudaCopyProcessNeighborZIndex(i, j);
					////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}
}

void GridReader::allocArrays_BoundaryQs()
{
	std::cout << "------read BoundaryQs-------" <<std::endl;

	std::vector<std::shared_ptr<BoundaryQs> > BC_Qs(channelDirections.size());
	this->makeReader(BC_Qs, para);

	int level = BC_Qs[0]->getLevel();
	BoundaryQs *obj_geomQ = new BoundaryQs(para->getgeomBoundaryBcQs(), para, "geo", false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Normals Geo
	BoundaryQs *obj_geomNormalX = new BoundaryQs(para->getgeomBoundaryNormalX(), para, "geoNormal", false);
	BoundaryQs *obj_geomNormalY = new BoundaryQs(para->getgeomBoundaryNormalY(), para, "geoNormal", false);
	BoundaryQs *obj_geomNormalZ = new BoundaryQs(para->getgeomBoundaryNormalZ(), para, "geoNormal", false);
	//////////////////////////////////////////////////////////////////////////
	//Normals Inflow
	BoundaryQs *obj_inflowNormalX = new BoundaryQs(para->getInflowBoundaryNormalX(), para, "inflowNormal", false);
	BoundaryQs *obj_inflowNormalY = new BoundaryQs(para->getInflowBoundaryNormalY(), para, "inflowNormal", false);
	BoundaryQs *obj_inflowNormalZ = new BoundaryQs(para->getInflowBoundaryNormalZ(), para, "inflowNormal", false);
	//////////////////////////////////////////////////////////////////////////
	//Normals Outflow
	BoundaryQs *obj_outflowNormalX = new BoundaryQs(para->getOutflowBoundaryNormalX(), para, "outflowNormal", false);
	BoundaryQs *obj_outflowNormalY = new BoundaryQs(para->getOutflowBoundaryNormalY(), para, "outflowNormal", false);
	BoundaryQs *obj_outflowNormalZ = new BoundaryQs(para->getOutflowBoundaryNormalZ(), para, "outflowNormal", false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//-----------------------------------------Vektoren Deklarationen-----------
	vector<vector<vector<real> > > noslipQs;
	vector<vector<unsigned int> > noslipIndex;
	vector<vector<vector<real> > > slipQs;
	vector<vector<unsigned int> > slipIndex;
	vector<vector<vector<real> > > pressureQs;
	vector<vector<unsigned int> > pressureIndex;
	vector<vector<vector<real> > > velocityQs;
	vector<vector<unsigned int> > velocityIndex;
	vector<vector<vector<real> > > outflowQs;
	vector<vector<unsigned int> > outflowIndex;
	//Geom-Daten werden mit dem Methode initArray() auf die Arrays kopiert. kein GetSet!

	//cout << "5: MyID: " << para->getMyID() << endl;
	noslipQs.resize(level + 1);
	slipQs.resize(level + 1);
	pressureQs.resize(level + 1);
	velocityQs.resize(level + 1);
	outflowQs.resize(level + 1);

	noslipIndex.resize(level + 1);
	slipIndex.resize(level + 1);
	pressureIndex.resize(level + 1);
	velocityIndex.resize(level + 1);
	outflowIndex.resize(level + 1);

	//cout << "6: MyID: " << para->getMyID() << endl;

	for (int i = 0; i <= level; i++) {
		noslipQs[i].resize(27);
		slipQs[i].resize(27);
		pressureQs[i].resize(27);
		velocityQs[i].resize(27);
		outflowQs[i].resize(27);
	}

	//cout << "7: MyID: " << para->getMyID() << endl;

	for (int i = 0; i < this->channelBoundaryConditions.size(); i++) {
		if (this->channelBoundaryConditions[i] == "noSlip") {  BC_Qs[i]->getQs(noslipQs);  BC_Qs[i]->getIndices(noslipIndex); }
		else if (this->channelBoundaryConditions[i] == "velocity") {  BC_Qs[i]->getQs(velocityQs);  BC_Qs[i]->getIndices(velocityIndex); }
		else if (this->channelBoundaryConditions[i] == "pressure") {  BC_Qs[i]->getQs(pressureQs);  BC_Qs[i]->getIndices(pressureIndex); }
		else if (this->channelBoundaryConditions[i] == "slip") {  BC_Qs[i]->getQs(slipQs);  BC_Qs[i]->getIndices(slipIndex); }
		else if (this->channelBoundaryConditions[i] == "outflow") {  BC_Qs[i]->getQs(outflowQs);  BC_Qs[i]->getIndices(outflowIndex); }
	}



	for (int i = 0; i <= level; i++) {
		int temp1 = (int)pressureIndex[i].size();
		if (temp1 > 0)
		{
			cout << "Groesse Pressure:  " << i << " : " << temp1 << endl;
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
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = pressureQs.begin(); it != pressureQs.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = pressureIndex.begin(); it != pressureIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QPress.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// advection - diffusion stuff
			//cout << "vor advec diff" << endl;
			if (para->getDiffOn() == true) {
				//////////////////////////////////////////////////////////////////////////
				//cout << "vor setzen von kTemp" << endl;
				para->getParH(i)->TempPress.kTemp = temp1;
				para->getParD(i)->TempPress.kTemp = temp1;
				cout << "Groesse TempPress.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
				//////////////////////////////////////////////////////////////////////////
				para->cudaAllocTempPressBC(i);
				//cout << "nach alloc" << endl;
				//////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp1; m++)
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



	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp3 = (int)velocityIndex[i].size();
		if (temp3 > 0)
		{
			cout << "Groesse velocity level " << i << " : " << temp3 << endl;
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
			//////////////////////////////////////////////////////////////////
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = velocityQs.begin(); it != velocityQs.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							Q.q27[j][n] = *it3;

							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = velocityIndex.begin(); it != velocityIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->Qinflow.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// advection - diffusion stuff
			if (para->getDiffOn() == true) {
				//////////////////////////////////////////////////////////////////////////
				para->getParH(i)->TempVel.kTemp = temp3;
				para->getParD(i)->TempVel.kTemp = temp3;
				cout << "Groesse TempVel.kTemp = " << para->getParH(i)->TempPress.kTemp << endl;
				cout << "getTemperatureInit = " << para->getTemperatureInit() << endl;
				cout << "getTemperatureBC = " << para->getTemperatureBC() << endl;
				//////////////////////////////////////////////////////////////////////////
				para->cudaAllocTempVeloBC(i);
				//cout << "nach alloc " << endl;
				//////////////////////////////////////////////////////////////////////////
				for (int m = 0; m < temp3; m++)
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
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyVeloBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getIsInflowNormal()) {
				int temp = obj_inflowNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten InflowBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QInflowNormalX.kQ = temp;
					para->getParH(i)->QInflowNormalY.kQ = temp;
					para->getParH(i)->QInflowNormalZ.kQ = temp;
					para->getParD(i)->QInflowNormalX.kQ = para->getParH(i)->QInflowNormalX.kQ;
					para->getParD(i)->QInflowNormalY.kQ = para->getParH(i)->QInflowNormalY.kQ;
					para->getParD(i)->QInflowNormalZ.kQ = para->getParH(i)->QInflowNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocInflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_inflowNormalX->setIndex(para->getParH(i)->QInflowNormalX.k, i);
					obj_inflowNormalY->setIndex(para->getParH(i)->QInflowNormalY.k, i);
					obj_inflowNormalZ->setIndex(para->getParH(i)->QInflowNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					real* QQX = para->getParH(i)->QInflowNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QInflowNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_inflowNormalX->setValues(QX.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Y
					real* QQY = para->getParH(i)->QInflowNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QInflowNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_inflowNormalY->setValues(QY.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Z
					real* QQZ = para->getParH(i)->QInflowNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QInflowNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_inflowNormalZ->setValues(QZ.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyInflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}//ende oberste for schleife


	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp = (int)outflowIndex[i].size();
		if (temp > 0)
		{
			cout << "Groesse Outflow:  " << i << " : " << temp << endl;
			//cout << "Groesse Pressure:  " << i << " : " << temp1 << "MyID: " << para->getMyID() << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			real* QQ = para->getParH(i)->Qoutflow.q27[0];
			unsigned int sizeQ = para->getParH(i)->Qoutflow.kQ;
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
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = outflowQs.begin(); it != outflowQs.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = outflowIndex.begin(); it != outflowIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->Qoutflow.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyOutflowBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (para->getIsOutflowNormal()) {
				int temp = obj_outflowNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten OutflowBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QOutflowNormalX.kQ = temp;
					para->getParH(i)->QOutflowNormalY.kQ = temp;
					para->getParH(i)->QOutflowNormalZ.kQ = temp;
					para->getParD(i)->QOutflowNormalX.kQ = para->getParH(i)->QOutflowNormalX.kQ;
					para->getParD(i)->QOutflowNormalY.kQ = para->getParH(i)->QOutflowNormalY.kQ;
					para->getParD(i)->QOutflowNormalZ.kQ = para->getParH(i)->QOutflowNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocOutflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_outflowNormalX->setIndex(para->getParH(i)->QOutflowNormalX.k, i);
					obj_outflowNormalY->setIndex(para->getParH(i)->QOutflowNormalY.k, i);
					obj_outflowNormalZ->setIndex(para->getParH(i)->QOutflowNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					real* QQX = para->getParH(i)->QOutflowNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QOutflowNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_outflowNormalX->setValues(QX.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Y
					real* QQY = para->getParH(i)->QOutflowNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QOutflowNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_outflowNormalY->setValues(QY.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Z
					real* QQZ = para->getParH(i)->QOutflowNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QOutflowNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_outflowNormalZ->setValues(QZ.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyOutflowNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}//ende if
	}//ende oberste for schleife





	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp2 = (int)noslipQs[i][0].size();
		para->getParH(i)->QWall.kQ = temp2;
		para->getParD(i)->QWall.kQ = para->getParH(i)->QWall.kQ;
		para->getParH(i)->kQ = para->getParH(i)->QWall.kQ;
		para->getParD(i)->kQ = para->getParH(i)->QWall.kQ;
		if (temp2 > 0)
		{
			cout << "Groesse NoSlip: " << i << " : " << temp2 << endl;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//para->getParH(i)->QWall.kQ = temp2;
			//para->getParD(i)->QWall.kQ = para->getParH(i)->QWall.kQ;
			//para->getParH(i)->kQ = para->getParH(i)->QWall.kQ;
			//para->getParD(i)->kQ = para->getParH(i)->QWall.kQ;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocWallBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			real* QQ = para->getParH(i)->QWall.q27[0];
			unsigned int sizeQ = para->getParH(i)->QWall.kQ;
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
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = noslipQs.begin(); it != noslipQs.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = noslipIndex.begin(); it != noslipIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QWall.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopyWallBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}//ende oberste for schleife



	 //--------------------------------------------------------------------------//
	for (int i = 0; i <= level; i++) {
		int temp2 = (int)slipQs[i][0].size();
        cout << "size Slip: " << i << " : " << temp2 << endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        para->getParH(i)->QSlip.kQ = temp2;
        para->getParD(i)->QSlip.kQ = para->getParH(i)->QSlip.kQ;
        para->getParH(i)->kSlipQ = para->getParH(i)->QSlip.kQ;
        para->getParD(i)->kSlipQ = para->getParH(i)->QSlip.kQ;
		if (temp2 > 0)
		{

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaAllocSlipBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//preprocessing
			real* QQ = para->getParH(i)->QSlip.q27[0];
			unsigned int sizeQ = para->getParH(i)->QSlip.kQ;
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
			int d = 0;
			int j = 0;
			int n = 0;
			for (vector<vector<vector<real> > >::iterator it = slipQs.begin(); it != slipQs.end(); it++) {
				if (i == d) {
					for (vector<vector<real> >::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						for (vector<real>::iterator it3 = it2->begin(); it3 != it2->end(); it3++) {
							Q.q27[j][n] = *it3;
							n++;
						}
						j++; // zaehlt die Spalte mit		
						n = 0;
					}
				}
				d++; // zaehlt das Level mit
				j = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			d = 0;
			n = 0;
			for (vector<vector<unsigned int> >::iterator it = slipIndex.begin(); it != slipIndex.end(); it++) {
				if (i == d) {
					for (vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
						para->getParH(i)->QSlip.k[n] = *it2;
						n++;
					}
				}
				d++;
				n = 0;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			para->cudaCopySlipBC(i);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}//ende oberste for schleife



	 //--------------------------------------------------------------------------//
	if (para->getIsGeo()) {
		for (int i = 0; i <= level; i++) {
			int temp4 = obj_geomQ->getSize(i);
			para->getParH(i)->QGeom.kQ = temp4;
			para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
			if (temp4 > 0)
			{
				cout << "Groesse der Daten GeomBoundaryQs, Level " << i << " : " << temp4 << endl;
				cout << "Groesse der Daten GeomBoundaryQs, Level:  " << i << " : " << temp4 << "MyID: " << para->getMyID() << endl;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//para->getParH(i)->QGeom.kQ = temp4;
				//para->getParD(i)->QGeom.kQ = para->getParH(i)->QGeom.kQ;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				para->cudaAllocGeomBC(i);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				//////////////////////////////////////////////////////////////////////////
				//Indexarray
				obj_geomQ->setIndex(para->getParH(i)->QGeom.k, i);
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
				for (int j = 0; j < 27; j++) {
					obj_geomQ->setValues(Q.q27, i);
				}//ende der for schleife
				 //////////////////////////////////////////////////////////////////
				for (int test = 0; test < temp4; test++)
				{
					Q.q27[dirZERO][test] = 0.0f;
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
					para->getParH(i)->Temp.kTemp = temp4;
					para->getParD(i)->Temp.kTemp = temp4;
					cout << "Groesse Temp.kTemp = " << para->getParH(i)->Temp.kTemp << endl;
					//////////////////////////////////////////////////////////////////////////
					para->cudaAllocTempNoSlipBC(i);
					//////////////////////////////////////////////////////////////////////////
					for (int m = 0; m < temp4; m++)
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
			if (para->getIsGeoNormal()) {
				int temp = obj_geomNormalX->getSize(i);
				if (temp > 0)
				{
					cout << "Groesse der Daten GeomBoundaryNormalsX, Level " << i << " : " << temp << endl;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->getParH(i)->QGeomNormalX.kQ = temp;
					para->getParH(i)->QGeomNormalY.kQ = temp;
					para->getParH(i)->QGeomNormalZ.kQ = temp;
					para->getParD(i)->QGeomNormalX.kQ = para->getParH(i)->QGeomNormalX.kQ;
					para->getParD(i)->QGeomNormalY.kQ = para->getParH(i)->QGeomNormalY.kQ;
					para->getParD(i)->QGeomNormalZ.kQ = para->getParH(i)->QGeomNormalZ.kQ;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaAllocGeomNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//Indexarray
					obj_geomNormalX->setIndex(para->getParH(i)->QGeomNormalX.k, i);
					obj_geomNormalY->setIndex(para->getParH(i)->QGeomNormalY.k, i);
					obj_geomNormalZ->setIndex(para->getParH(i)->QGeomNormalZ.k, i);
					//////////////////////////////////////////////////////////////////

					//////////////////////////////////////////////////////////////////////////
					//preprocessing X
					real* QQX = para->getParH(i)->QGeomNormalX.q27[0];
					unsigned int sizeQX = para->getParH(i)->QGeomNormalX.kQ;
					QforBoundaryConditions QX;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QX.q27[j] = &QQX[j * sizeQX];
						obj_geomNormalX->setValues(QX.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Y
					real* QQY = para->getParH(i)->QGeomNormalY.q27[0];
					unsigned int sizeQY = para->getParH(i)->QGeomNormalY.kQ;
					QforBoundaryConditions QY;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QY.q27[j] = &QQY[j * sizeQY];
						obj_geomNormalY->setValues(QY.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 //////////////////////////////////////////////////////////////////////////
					 //preprocessing Z
					real* QQZ = para->getParH(i)->QGeomNormalZ.q27[0];
					unsigned int sizeQZ = para->getParH(i)->QGeomNormalZ.kQ;
					QforBoundaryConditions QZ;
					//////////////////////////////////////////////////////////////////
					for (int j = 0; j < 27; j++) {
						QZ.q27[j] = &QQZ[j * sizeQZ];
						obj_geomNormalZ->setValues(QZ.q27, i);
					}//ende der for schleife
					 //////////////////////////////////////////////////////////////////

					 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					para->cudaCopyGeomNormals(i);
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}

			}
		}
	}//ende oberstes for
	 //--------------------------------------------------------------------------//





	//for (int i = 0; i < channelBoundaryConditions.size(); i++)
	//{
	//	if (this->channelBoundaryConditions[i] == "noSlip") { setNoSlipQs(BC_Qs[i]); }
	//	else if (this->channelBoundaryConditions[i] == "velocity") { setVelocityQs(BC_Qs[i]); }
	//	else if (this->channelBoundaryConditions[i] == "pressure") { setPressQs(BC_Qs[i]); }
	//	else if (this->channelBoundaryConditions[i] == "outflow") { setOutflowQs(BC_Qs[i]); }
	//}

	//std::shared_ptr<BoundaryQs> obj_geomQ = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getgeomBoundaryBcQs(), para, "geo", false));
	//if (para->getIsGeo())
	//	setGeoQs(obj_geomQ);

	std::cout << "-----finish BoundaryQs------" <<std::endl;
}


/*------------------------------------------------------------------------------------------------*/
/*----------------------------------------q setter methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridReader::setPressQs(std::shared_ptr<BoundaryQs> boundaryQ) const
{
	for (unsigned int level = 0; level <= boundaryQ->getLevel(); level++)
	{
		if (hasQs(boundaryQ, level))
		{
			this->printQSize("pressure", boundaryQ, level);
			this->initalQStruct(para->getParH(level)->QPress, boundaryQ, level);
            cudaMemoryManager->cudaCopyPress(level);
		}
	}
}

void GridReader::setVelocityQs(std::shared_ptr<BoundaryQs> boundaryQ) const
{
	for (unsigned int level = 0; level <= boundaryQ->getLevel(); level++)
	{
		if (hasQs(boundaryQ, level))
		{
			this->printQSize("velocity", boundaryQ, level);
			this->initalQStruct(para->getParH(level)->Qinflow, boundaryQ, level);
            cudaMemoryManager->cudaCopyVeloBC(level);
		}
	}
}

void GridReader::setOutflowQs(std::shared_ptr<BoundaryQs> boundaryQ) const
{
	for (unsigned int level = 0; level <= boundaryQ->getLevel(); level++)
	{
		if (hasQs(boundaryQ, level))
		{
			this->printQSize("outflow", boundaryQ, level);
			this->initalQStruct(para->getParH(level)->Qoutflow, boundaryQ, level);
            cudaMemoryManager->cudaCopyOutflowBC(level);
		}
	}
}

void GridReader::setNoSlipQs(std::shared_ptr<BoundaryQs> boundaryQ) const
{
	for (unsigned int level = 0; level <= boundaryQ->getLevel(); level++)
	{
		if (hasQs(boundaryQ, level))
		{
			this->printQSize("no slip", boundaryQ, level);
			this->setSizeNoSlip(boundaryQ, level);
			this->initalQStruct(para->getParH(level)->QWall, boundaryQ, level);
            cudaMemoryManager->cudaCopyWallBC(level);
		}
	}
}

void GridReader::setGeoQs(std::shared_ptr<BoundaryQs> boundaryQ) const
{
	for (unsigned int level = 0; level <= boundaryQ->getLevel(); level++)
	{
		if (hasQs(boundaryQ, level))
		{
			this->printQSize("geo Qs", boundaryQ, level);
			this->setSizeGeoQs(boundaryQ, level);
			this->initalQStruct(para->getParH(level)->QGeom, boundaryQ, level);

			modifyQElement(boundaryQ, level);

            cudaMemoryManager->cudaCopyGeomBC(level);
		}
	}
}

void GridReader::modifyQElement(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	QforBoundaryConditions Q;
	real* QQ = para->getParH(level)->QGeom.q27[0];
	Q.q27[dirZERO] = &QQ[dirZERO * para->getParH(level)->QGeom.kQ];
	for (unsigned int i = 0; i < boundaryQ->getSize(level); i++)
		Q.q27[dirZERO][i] = 0.0f;
}

/*------------------------------------------------------------------------------------------------*/
/*---------------------------------------private q methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridReader::initalQStruct(QforBoundaryConditions& Q, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	QforBoundaryConditions qTemp;
	this->setQ27Size(qTemp, Q.q27[0], Q.kQ);
	boundaryQ->setValues(qTemp.q27, level);
	boundaryQ->setIndex(Q.k, level);
}

bool GridReader::hasQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	return boundaryQ->getSize(level) > 0;
}

void GridReader::initalGridInformations()
{

}

void GridReader::setQ27Size(QforBoundaryConditions &Q, real* QQ, unsigned int sizeQ) const
{
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
}

void GridReader::setSizeNoSlip(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	para->getParH(level)->QWall.kQ = boundaryQ->getSize(level);
	para->getParD(level)->QWall.kQ = para->getParH(level)->QWall.kQ;
	para->getParH(level)->kQ = para->getParH(level)->QWall.kQ;
	para->getParD(level)->kQ = para->getParH(level)->QWall.kQ;
    cudaMemoryManager->cudaAllocWallBC(level);
}

void GridReader::setSizeGeoQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	para->getParH(level)->QGeom.kQ = boundaryQ->getSize(level);
	para->getParD(level)->QGeom.kQ = para->getParH(level)->QGeom.kQ;

    cudaMemoryManager->cudaAllocGeomBC(level);
}

void GridReader::printQSize(std::string bc, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const
{
	std::cout << "level " << level << ", " << bc << "-size: " << boundaryQ->getSize(level) << std::endl;
}


void GridReader::setDimensions()
{
	std::ifstream numberNodes;
	numberNodes.open(para->getnumberNodes().c_str(), std::ios::in);
	if (!numberNodes) {
		std::cerr << "can't open file NumberNodes: " << para->getnumberNodes() << std::endl;
		exit(1);
	}

	std::string buffer;
	int bufferInt;
	std::vector<int> localGridNX;
	std::vector<int> localGridNY;
	std::vector<int> localGridNZ;

	for (/*unsigned*/ int i = 0; i <= para->getMaxLevel(); i++) {
		numberNodes >> buffer;
		numberNodes >> bufferInt;
		localGridNX.push_back(bufferInt);
		numberNodes >> bufferInt;
		localGridNY.push_back(bufferInt);
		numberNodes >> bufferInt;
		localGridNZ.push_back(bufferInt);
	}
	para->setGridX(localGridNX);
	para->setGridY(localGridNY);
	para->setGridZ(localGridNZ);
}

void GridReader::setBoundingBox()
{
	std::ifstream numberNodes;
	numberNodes.open(para->getLBMvsSI().c_str(), std::ios::in);
	if (!numberNodes) {
		std::cerr << "can't open file LBMvsSI" << std::endl;
		exit(1);
	}
	real bufferreal;
	std::vector<real> minX, maxX, minY, maxY, minZ, maxZ;

	for (int i = 0; i <= para->getMaxLevel(); i++) {
		numberNodes >> bufferreal;
		minX.push_back(bufferreal);
		numberNodes >> bufferreal;
		minY.push_back(bufferreal);
		numberNodes >> bufferreal;
		minZ.push_back(bufferreal);
		numberNodes >> bufferreal;
		maxX.push_back(bufferreal);
		numberNodes >> bufferreal;
		maxY.push_back(bufferreal);
		numberNodes >> bufferreal;
		maxZ.push_back(bufferreal);
	}
	para->setMinCoordX(minX);
	para->setMinCoordY(minY);
	para->setMinCoordZ(minZ);
	para->setMaxCoordX(maxX);
	para->setMaxCoordY(maxY);
	para->setMaxCoordZ(maxZ);
}

void GridReader::initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex,  std::string boundaryCondition)
{
	std::vector<unsigned int>neighVec;
	std::vector<unsigned int>indexVec;
	
	int counter = 0;

	for(unsigned int i=0; i<neighX->getLevel();i++) {
		if(boundaryCondition =="periodic_y"){
			neighVec = neighY->getVec(i);
		} 
		else if(boundaryCondition =="periodic_x"){
			neighVec = neighX->getVec(i);
		}
		else if(boundaryCondition =="periodic_z"){
			neighVec = neighZ->getVec(i);
		}
		else {
			std::cout << "wrong String in periodicValue" << std::endl;
			exit(1);
		}

		for (std::vector<unsigned int>::iterator it = periodIndex[i].begin(); it != periodIndex[i].end(); it++) {
			if(periodV[i][0][counter] != 0) {
				neighVec[*it]=periodV[i][0][counter];
			}

			counter++;
		}


		if(boundaryCondition =="periodic_y"){
			neighY->setVec(i, neighVec);
		} 
		else if(boundaryCondition =="periodic_x"){
			neighX->setVec(i, neighVec);
		}
		else if(boundaryCondition =="periodic_z"){
			neighZ->setVec(i, neighVec);
		}

	}
}

void GridReader::makeReader(std::shared_ptr<Parameter> para)
{
	for (int i = 0; i < BC_Values.size(); i++)
	{
		if (channelDirections[i].compare("inlet") == 0){ BC_Values[i]  = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->getinletBcValues())); }
		if (channelDirections[i].compare("outlet") == 0){ BC_Values[i] = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->getoutletBcValues())); }
		if (channelDirections[i].compare("back") == 0){ BC_Values[i]   = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->getbackBcValues())); }
		if (channelDirections[i].compare("front") == 0){ BC_Values[i]  = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->getfrontBcValues())); }
		if (channelDirections[i].compare("top") == 0){ BC_Values[i]    = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->gettopBcValues())); }
		if (channelDirections[i].compare("bottom") == 0){ BC_Values[i] = std::shared_ptr<BoundaryValues>(new BoundaryValues(para->getbottomBcValues()));}
	}
}

void GridReader::makeReader(std::vector<std::shared_ptr<BoundaryQs> > &BC_Qs, std::shared_ptr<Parameter> para)
{
	for (int i = 0; i < BC_Qs.size(); i++)
	{
		if (channelDirections[i].compare("inlet") == 0){ BC_Qs[i]  = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getinletBcQs(), false)); }
		if (channelDirections[i].compare("outlet") == 0){ BC_Qs[i] = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getoutletBcQs(), false)); }
		if (channelDirections[i].compare("back") == 0){ BC_Qs[i]   = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getbackBcQs(), false)); }
		if (channelDirections[i].compare("front") == 0){ BC_Qs[i]  = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getfrontBcQs(), false)); }
		if (channelDirections[i].compare("top") == 0){ BC_Qs[i]    = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->gettopBcQs(), false)); }
		if (channelDirections[i].compare("bottom") == 0){ BC_Qs[i] = std::shared_ptr<BoundaryQs>(new BoundaryQs(para->getbottomBcQs(), false)); }
	}
}

void GridReader::setChannelBoundaryCondition()
{
	for (int i = 0; i < channelDirections.size(); i++)
	{
		this->channelBoundaryConditions[i] = BC_Values[i]->getBoundaryCondition();
		std::cout << this->channelDirections[i] << " Boundary: " << channelBoundaryConditions[i] << std::endl;
	}
}