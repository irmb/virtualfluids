#include "GridGenerator.h"

#include "Parameter/Parameter.h"
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GPU/CudaMemoryManager.h>

#include <sstream>
#include <iostream>
#include "utilities/math/Math.h"

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

		builder->getNodeValues(
			para->getParH(level)->coordX_SP,
			para->getParH(level)->coordY_SP,
			para->getParH(level)->coordZ_SP,
			para->getParH(level)->neighborX_SP,
			para->getParH(level)->neighborY_SP,
			para->getParH(level)->neighborZ_SP,
			para->getParH(level)->geoSP,
			level);

   /*     if (true) {
            const real delta = para->getParH(level)->coordX_SP[2] - para->getParH(level)->coordX_SP[1];
            for (uint i = 0; i < numberOfNodesPerLevel; i++)
            {
                const real coordX = para->getParH(level)->coordX_SP[i];
                const real coordY = para->getParH(level)->coordY_SP[i];
                const real coordZ = para->getParH(level)->coordZ_SP[i];

                const uint neighborX = para->getParH(level)->neighborX_SP[i];
                const uint neighborY = para->getParH(level)->neighborY_SP[i];
                const uint neighborZ = para->getParH(level)->neighborZ_SP[i];

                const uint type = para->getParH(level)->geoSP[i];

                if (coordX == 40 || coordY ==1 || coordY == 40 || coordZ == 1 || coordZ == 40)
                    continue;

                real expectedXNeighborX;
                if (level == 0 && coordX == 39)
                    expectedXNeighborX = 2;
                else
                    expectedXNeighborX = coordX + delta;

                const real expectedXNeighborY = coordY;
                const real expectedXNeighborZ = coordZ;

                const real actualXNeighborX = para->getParH(level)->coordX_SP[neighborX];
                const real actualXNeighborY = para->getParH(level)->coordY_SP[neighborX];
                const real actualXNeighborZ = para->getParH(level)->coordZ_SP[neighborX];

                const bool equal = (expectedXNeighborX == actualXNeighborX) && (expectedXNeighborY == actualXNeighborY) && (expectedXNeighborZ == actualXNeighborZ);
                if(!equal)
                {
                    printf("coordinate: %2.2f, %2.2f, %2.2f   expectedX neighbor: %2.2f, %2.2f, %2.2f    actual X neighbor: %2.2f, %2.2f, %2.2f;   type: %d\n", coordX, coordY, coordZ, expectedXNeighborX, expectedXNeighborY, expectedXNeighborZ, actualXNeighborX, actualXNeighborY, actualXNeighborZ, type);
                }
            }
        }*/



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

	channelBoundaryConditions = builder->getTypeOfBoundaryConditions();

	for (int i = 0; i < channelBoundaryConditions.size(); i++)
	{
        setVelocityValues(i);
        setPressureValues(i); 
        setOutflowValues(i); 
	}
}

void GridGenerator::setPressureValues(int channelSide) const
{
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
	{

        int sizePerLevel = 0;
        if ((this->channelBoundaryConditions[channelSide] == "pressure"))
		    sizePerLevel = builder->getBoundaryConditionSize(channelSide);

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

void GridGenerator::setPressRhoBC(int sizePerLevel, int level, int channelSide) const
{
	builder->setPressValues(para->getParH(level)->QPress.RhoBC, para->getParH(level)->QPress.kN, channelSide, level);
	for (int m = 0; m < sizePerLevel; m++)
		para->getParH(level)->QPress.RhoBC[m] = (para->getParH(level)->QPress.RhoBC[m] / para->getFactorPressBC());
}


void GridGenerator::setVelocityValues(int channelSide) const
{
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
	{
        int sizePerLevel = 0;
        if ((this->channelBoundaryConditions[channelSide] == "velocity"))
            sizePerLevel = builder->getBoundaryConditionSize(channelSide);

		setVelocitySizePerLevel(level, sizePerLevel);

        if (sizePerLevel > 0)
        {
            std::cout << "size velocity level " << level << " : " << sizePerLevel << std::endl;

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
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
	{
        int sizePerLevel = 0;
        if ((this->channelBoundaryConditions[channelSide] == "outflow"))
            sizePerLevel = builder->getBoundaryConditionSize(channelSide);

			setOutflowSizePerLevel(level, sizePerLevel);
            if (sizePerLevel > 0)
            {
                std::cout << "size outflow level " << level << " : " << sizePerLevel << std::endl;

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
		para->getParH(level)->Qoutflow.RhoBC[index] = (para->getParH(level)->Qoutflow.RhoBC[index] / para->getFactorPressBC()) * (real)0.0;
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
		cudaMemoryManager->cudaAllocInterfaceCF(level);
		cudaMemoryManager->cudaAllocInterfaceFC(level);
		cudaMemoryManager->cudaAllocInterfaceOffCF(level);
		cudaMemoryManager->cudaAllocInterfaceOffFC(level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //init
        builder->setOffsetCF(para->getParH(level)->offCF.xOffCF, para->getParH(level)->offCF.yOffCF, para->getParH(level)->offCF.zOffCF, level);
        builder->setOffsetFC(para->getParH(level)->offFC.xOffFC, para->getParH(level)->offFC.yOffFC, para->getParH(level)->offFC.zOffFC, level);
        builder->getGridInterfaceIndices(para->getParH(level)->intCF.ICellCFC, para->getParH(level)->intCF.ICellCFF, para->getParH(level)->intFC.ICellFCC, para->getParH(level)->intFC.ICellFCF, level);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //copy
		cudaMemoryManager->cudaCopyInterfaceCF(level);
		cudaMemoryManager->cudaCopyInterfaceFC(level);
		cudaMemoryManager->cudaCopyInterfaceOffCF(level);
		cudaMemoryManager->cudaCopyInterfaceOffFC(level);
    }
}


/*------------------------------------------------------------------------------------------------*/
/*----------------------------------------q setter methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridGenerator::setPressQs(int channelSide) const
{
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
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
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
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
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
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
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
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
	for (uint level = 0; level < builder->getNumberOfGridLevels(); level++)
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

void GridGenerator::modifyQElement(int channelSide,  uint level) const
{
	QforBoundaryConditions Q;
	real* QQ = para->getParH(level)->QGeom.q27[0];
	Q.q27[dirZERO] = &QQ[dirZERO * para->getParH(level)->QGeom.kQ];
	for (int i = 0; i < builder->getBoundaryConditionSize(channelSide); i++)
		Q.q27[dirZERO][i] = 0.0f;
}

/*------------------------------------------------------------------------------------------------*/
/*---------------------------------------private q methods----------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
void GridGenerator::initalQStruct(QforBoundaryConditions& Q, int channelSide, uint level) const
{
	QforBoundaryConditions qTemp;
	this->setQ27Size(qTemp, Q.q27[0], Q.kQ);
	builder->setQs(qTemp.q27, Q.k, channelSide, level);
}

bool GridGenerator::hasQs(int channelSide, uint level) const
{
	return builder->getBoundaryConditionSize(channelSide) > 0;
}


void GridGenerator::setQ27Size(QforBoundaryConditions &Q, real* QQ, uint sizeQ) const
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

void GridGenerator::setSizeNoSlip(int channelSide, uint level) const
{
	para->getParH(level)->QWall.kQ = builder->getBoundaryConditionSize(channelSide);
	para->getParD(level)->QWall.kQ = para->getParH(level)->QWall.kQ;
	para->getParH(level)->kQ = para->getParH(level)->QWall.kQ;
	para->getParD(level)->kQ = para->getParH(level)->QWall.kQ;
    cudaMemoryManager->cudaAllocWallBC(level);
}

void GridGenerator::setSizeGeoQs(uint level) const
{
	para->getParH(level)->QGeom.kQ = builder->getBoundaryConditionSize(GEOMQS);
	para->getParD(level)->QGeom.kQ = para->getParH(level)->QGeom.kQ;

    cudaMemoryManager->cudaAllocGeomBC(level);
}

void GridGenerator::printQSize(std::string bc,int channelSide, uint level) const
{
	std::cout << "level " << level << ", " << bc << "-size: " << builder->getBoundaryConditionSize(channelSide) << std::endl;
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
