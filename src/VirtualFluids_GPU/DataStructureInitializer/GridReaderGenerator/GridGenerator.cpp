#include "GridGenerator.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GPU/CudaMemoryManager.h>

#include <iostream>

GridGenerator::GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para)
{
	this->builder = builder;
    this->para = para;
    this->cudaMemoryManager = CudaMemoryManager::make(para);
}

GridGenerator::~GridGenerator()
{

}

void GridGenerator::setUnstructuredGridBuilder(std::shared_ptr<GridBuilder> builder)
{
	this->builder = builder;
}

void GridGenerator::allocArrays_CoordNeighborGeo()
{
	int numberOfLevels = 0;//coordX.getLevel();
	std::cout << "Number of Level: " << numberOfLevels + 1 << std::endl;
	int numberOfNodesGlobal = 0;
	std::cout << "Number of Nodes: " << std::endl;
	
	for (int level = 0; level <= numberOfLevels; level++) 
	{
		int numberOfNodesPerLevel = builder->getNumberOfNodes(level) + 1;
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
	}
	std::cout << "Number of Nodes: " << numberOfNodesGlobal << std::endl;
	std::cout << "-----finish Coord, Neighbor, Geo------" << std::endl;
}

void GridGenerator::allocArrays_BoundaryValues()
{
	std::cout << "------read BoundaryValues------" << std::endl;

	channelBoundaryConditions = builder->getTypeOfBoundaryConditions();

	for (int i = 0; i < channelBoundaryConditions.size(); i++)
	{
		if (this->channelBoundaryConditions[i] == "velocity") { setVelocityValues(i); }
		else if (this->channelBoundaryConditions[i] == "pressure") { setPressureValues(i); }
		else if (this->channelBoundaryConditions[i] == "outflow") { setOutflowValues(i); }
	}
}

void GridGenerator::setPressureValues(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		int sizePerLevel = builder->getBoundaryConditionSize(channelSide);
		if (sizePerLevel > 1)
		{
			std::cout << "size pressure level " << level << " : " << sizePerLevel << std::endl;

			setPressSizePerLevel(level, sizePerLevel);
            cudaMemoryManager->cudaAllocPress(level);

			setPressRhoBC(sizePerLevel, level, channelSide);
            cudaMemoryManager->cudaCopyPress(level);
		}
	}
}

void GridGenerator::setPressRhoBC(int sizePerLevel, int level, int channelSide) const
{
	builder->setPressValues(para->getParH(level)->QPress.RhoBC, para->getParH(level)->QPress.kN, channelSide, level);
	for (int m = 0; m < sizePerLevel; m++)
		para->getParH(level)->QPress.RhoBC[m] = (para->getParH(level)->QPress.RhoBC[m] / para->getFactorPressBC());
}


void GridGenerator::setVelocityValues(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
        int sizePerLevel = builder->getBoundaryConditionSize(channelSide);
        if (sizePerLevel > 1)
		{
			std::cout << "size velocity level " << level << " : " << sizePerLevel << std::endl;

			setVelocitySizePerLevel(level, sizePerLevel);
            cudaMemoryManager->cudaAllocVeloBC(level);

			setVelocity(level,  sizePerLevel, channelSide);
            cudaMemoryManager->cudaCopyVeloBC(level);
		}
	}
}

void GridGenerator::setVelocity(int level, int sizePerLevel, int channelSide) const
{
	builder->setVelocityValues(para->getParH(level)->Qinflow.Vx, para->getParH(level)->Qinflow.Vy, para->getParH(level)->Qinflow.Vz, channelSide,  level);

	for (int index = 0; index < sizePerLevel; index++)
	{
		//para->getParH(i)->Qinflow.Vx[m] = para->getParH(i)->Qinflow.Vx[m] / para->getVelocityRatio();
		//para->getParH(i)->Qinflow.Vy[m] = para->getParH(i)->Qinflow.Vy[m] / para->getVelocityRatio();
		//para->getParH(i)->Qinflow.Vz[m] = para->getParH(i)->Qinflow.Vz[m] / para->getVelocityRatio();
		para->getParH(level)->Qinflow.Vx[index] = para->getVelocity();//0.035;
		para->getParH(level)->Qinflow.Vy[index] = 0.0;//para->getVelocity();//0.0;
		para->getParH(level)->Qinflow.Vz[index] = 0.0;
	}
}

void GridGenerator::setOutflowValues(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		int sizePerLevel = builder->getBoundaryConditionSize(channelSide);
		if (sizePerLevel > 1)
		{
			std::cout << "size outflow level " << level << " : " << sizePerLevel << std::endl;

			setOutflowSizePerLevel(level, sizePerLevel);
            cudaMemoryManager->cudaAllocOutflowBC(level);

			setOutflow(level, sizePerLevel, channelSide);
            cudaMemoryManager->cudaCopyOutflowBC(level);

		}
	}
}

void GridGenerator::setOutflow(int level, int sizePerLevel, int channelSide) const
{
	builder->setOutflowValues(para->getParH(level)->Qoutflow.RhoBC, para->getParH(level)->Qoutflow.kN, channelSide, level);
	for (int index = 0; index < sizePerLevel; index++)
		para->getParH(level)->Qoutflow.RhoBC[index] = (para->getParH(level)->Qoutflow.RhoBC[index] / para->getFactorPressBC()) * (doubflo)0.0;
}



void GridGenerator::allocArrays_BoundaryQs()
{
	std::cout << "------read BoundaryQs-------" << std::endl;

	channelBoundaryConditions = builder->getTypeOfBoundaryConditions();

	for (int i = 0; i < channelBoundaryConditions.size(); i++)
	{
		if (this->channelBoundaryConditions[i] == "noSlip") { setNoSlipQs(i); }
		else if (this->channelBoundaryConditions[i] == "velocity") { setVelocityQs(i); }
		else if (this->channelBoundaryConditions[i] == "pressure") { setPressQs(i); }
		else if (this->channelBoundaryConditions[i] == "outflow") { setOutflowQs(i); }
	}

	if (para->getIsGeo())
		setGeoQs();

	std::cout << "-----finish BoundaryQs------" << std::endl;
}


/*------------------------------------------------------------------------------------------------*/
/*----------------------------------------q setter methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridGenerator::setPressQs(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		if (hasQs(channelSide, level))
		{
			this->printQSize("pressure", channelSide, level);
			this->initalQStruct(para->getParH(level)->QPress, channelSide, level);
            cudaMemoryManager->cudaCopyPress(level);
		}
	}
}

void GridGenerator::setVelocityQs(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		if (hasQs(channelSide, level))
		{
			this->printQSize("velocity", channelSide, level);
			this->initalQStruct(para->getParH(level)->Qinflow, channelSide, level);
            cudaMemoryManager->cudaCopyVeloBC(level);
		}
	}
}

void GridGenerator::setOutflowQs(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		if (hasQs(channelSide, level))
		{
			this->printQSize("outflow", channelSide, level);
			this->initalQStruct(para->getParH(level)->Qoutflow, channelSide, level);
            cudaMemoryManager->cudaCopyOutflowBC(level);
		}
	}
}

void GridGenerator::setNoSlipQs(int channelSide) const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		if (hasQs(channelSide, level))
		{
			this->printQSize("no slip", channelSide, level);
			this->setSizeNoSlip(channelSide, level);
			this->initalQStruct(para->getParH(level)->QWall, channelSide, level);
            cudaMemoryManager->cudaCopyWallBC(level);
		}
	}
}

void GridGenerator::setGeoQs() const
{
	for (unsigned int level = 0; level <= 0; level++)
	{
		if (hasQs(GEOMQS, level))
		{
			this->printQSize("geo Qs", GEOMQS, level);
			this->setSizeGeoQs(level);
			this->initalQStruct(para->getParH(level)->QGeom, GEOMQS, level);

			modifyQElement(GEOMQS, level);

            cudaMemoryManager->cudaCopyGeomBC(level);
		}
	}
}

void GridGenerator::modifyQElement(int channelSide,  unsigned int level) const
{
	QforBoundaryConditions Q;
	doubflo* QQ = para->getParH(level)->QGeom.q27[0];
	Q.q27[dirZERO] = &QQ[dirZERO * para->getParH(level)->QGeom.kQ];
	for (int i = 0; i < builder->getBoundaryConditionSize(channelSide); i++)
		Q.q27[dirZERO][i] = 0.0f;
}

/*------------------------------------------------------------------------------------------------*/
/*---------------------------------------private q methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridGenerator::initalQStruct(QforBoundaryConditions& Q, int channelSide, unsigned int level) const
{
	QforBoundaryConditions qTemp;
	this->setQ27Size(qTemp, Q.q27[0], Q.kQ);
	builder->setQs(qTemp.q27, Q.k, channelSide, level);
}

bool GridGenerator::hasQs(int channelSide, unsigned int level) const
{
	return builder->getBoundaryConditionSize(channelSide) > 0;
}

void GridGenerator::setQ27Size(QforBoundaryConditions &Q, doubflo* QQ, unsigned int sizeQ) const
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

void GridGenerator::setSizeNoSlip(int channelSide, unsigned int level) const
{
	para->getParH(level)->QWall.kQ = builder->getBoundaryConditionSize(channelSide);
	para->getParD(level)->QWall.kQ = para->getParH(level)->QWall.kQ;
	para->getParH(level)->kQ = para->getParH(level)->QWall.kQ;
	para->getParD(level)->kQ = para->getParH(level)->QWall.kQ;
    cudaMemoryManager->cudaAllocWallBC(level);
}

void GridGenerator::setSizeGeoQs(unsigned int level) const
{
	para->getParH(level)->QGeom.kQ = builder->getBoundaryConditionSize(GEOMQS);
	para->getParD(level)->QGeom.kQ = para->getParH(level)->QGeom.kQ;

    cudaMemoryManager->cudaAllocGeomBC(level);
}

void GridGenerator::printQSize(std::string bc,int channelSide, unsigned int level) const
{
	std::cout << "level " << level << ", " << bc << "-size: " << builder->getBoundaryConditionSize(channelSide) << std::endl;
}


void GridGenerator::setDimensions()
{
	std::vector<int> localGridNX(1);
	std::vector<int> localGridNY(1);
	std::vector<int> localGridNZ(1);

	builder->getDimensions(localGridNX[0], localGridNY[0], localGridNZ[0], 0);

	para->setGridX(localGridNX);
	para->setGridY(localGridNY);
	para->setGridZ(localGridNZ);
}

void GridGenerator::setBoundingBox()
{
	std::vector<int> localGridNX(1);
	std::vector<int> localGridNY(1);
	std::vector<int> localGridNZ(1);
	builder->getDimensions(localGridNX[0], localGridNY[0], localGridNZ[0], 0);

	std::vector<doubflo> minX, maxX, minY, maxY, minZ, maxZ;
	minX.push_back(0);
	minY.push_back(0);
	minZ.push_back(0);

	maxX.push_back((doubflo)localGridNX[0]);
	maxY.push_back((doubflo)localGridNY[0]);
	maxZ.push_back((doubflo)localGridNZ[0]);

	para->setMinCoordX(minX);
	para->setMinCoordY(minY);
	para->setMinCoordZ(minZ);
	para->setMaxCoordX(maxX);
	para->setMaxCoordY(maxY);
	para->setMaxCoordZ(maxZ);
}

void GridGenerator::initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way)
{

}

