#define _CRT_SECURE_NO_DEPRECATE
#include "SimulationFileWriter.h"

#include <iostream>
#include <iomanip>
#include <omp.h>
#include "stdint.h"

#include <GridGenerator/utilities/Transformator/Transformator.h>
#include "SimulationFileNames.h"

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/geometries/Vertex/Vertex.cuh>

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>


/*#################################################################################*/
/*---------------------------------public methods----------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::writeSimulationFiles(std::string folder, std::shared_ptr<GridBuilder> builder, bool binaer, std::shared_ptr<Transformator> trans)
{
    std::cout << "Write Simulation Files ... \n";
    SimulationFileWriter::folder = folder;
    clock_t  begin = clock();

    openFiles();
    writeLevelAndLevelSize((int)builder->getNumberOfNodes(0), builder->getQsValues());
    writeCoordFiles(binaer, builder, trans);
    writeBoundaryQsFile(builder->getQsValues());
    closeFiles();

    clock_t  end = clock();
    real time = real(end - begin) / CLOCKS_PER_SEC;
    std::cout << "time write files: " << time << " sec" << std::endl;
    std::cout << "... finish writing Simulation Files!\n";
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::writeCoordFiles(bool binaer, std::shared_ptr<GridBuilder> builder, std::shared_ptr<Transformator> trans)
{
    for (uint i = 0; i < builder->getNumberOfNodes(0); i++) {
        writeCoordsNeighborsGeo(i, binaer, builder, trans);
    }
}

void SimulationFileWriter::writeBoundaryQsFile(std::vector<std::vector<std::vector<real> > > qFiles)
{
    for (int rb = 0; rb < QFILES; rb++) {
        for (int index = 0; index < qFiles[rb].size(); index++) {
            writeBoundary(qFiles[rb][index], rb);
        }
    }
}

void SimulationFileWriter::openFiles()
{
    std::string path = folder;
    xCoordFile.open((path + simulationFileNames::coordX).c_str(), std::ios::out | std::ios::binary);
    yCoordFile.open((path + simulationFileNames::coordY).c_str(), std::ios::out | std::ios::binary);
    zCoordFile.open((path + simulationFileNames::coordZ).c_str(), std::ios::out | std::ios::binary);
    xNeighborFile.open((path + simulationFileNames::neighborX).c_str(), std::ios::out | std::ios::binary);
    yNeighborFile.open((path + simulationFileNames::neighborY).c_str(), std::ios::out | std::ios::binary);
    zNeighborFile.open((path + simulationFileNames::neighborZ).c_str(), std::ios::out | std::ios::binary);
    geoVecFile.open((path + simulationFileNames::geoVec).c_str(), std::ios::out | std::ios::binary);

    std::vector<std::string> qNames;
    qNames.push_back(path + simulationFileNames::inletBoundaryQ);
    qNames.push_back(path + simulationFileNames::outletBoundaryQ);
    qNames.push_back(path + simulationFileNames::topBoundaryQ);
    qNames.push_back(path + simulationFileNames::bottomBoundaryQ);
    qNames.push_back(path + simulationFileNames::frontBoundaryQ);
    qNames.push_back(path + simulationFileNames::backBoundaryQ);
	qNames.push_back(path + simulationFileNames::geomBoundaryQ);

    std::vector<std::string> valueNames;
    valueNames.push_back(path + simulationFileNames::inletBoundaryValues);
    valueNames.push_back(path + simulationFileNames::outletBoundaryValues);
    valueNames.push_back(path + simulationFileNames::topBoundaryValues);
    valueNames.push_back(path + simulationFileNames::bottomBoundaryValues);
    valueNames.push_back(path + simulationFileNames::frontBoundaryValues);
    valueNames.push_back(path + simulationFileNames::backBoundaryValues);
	valueNames.push_back(path + simulationFileNames::geomBoundaryValues);

    for (int i = 0; i < QFILES; i++){
        std::shared_ptr<std::ofstream> outQ(new std::ofstream);
        outQ->open(qNames[i].c_str(), std::ios::out | std::ios::binary);
        qStreams.push_back(outQ);

        std::shared_ptr<std::ofstream> outV(new std::ofstream);
        outV->open(valueNames[i].c_str(), std::ios::out | std::ios::binary);
        valueStreams.push_back(outV);
    }
}

void SimulationFileWriter::writeLevelAndLevelSize(int sizeCoords, std::vector<std::vector<std::vector<real> > > qFiles)
{
    std::string zeroIndex = "0 ";
    std::string zeroGeo = "16 ";
    std::string level = "0";
    
    xCoordFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    yCoordFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    zCoordFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    xNeighborFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    yNeighborFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    zNeighborFile << level << "\n" << sizeCoords << "\n" << zeroIndex;
    geoVecFile << level << "\n" << sizeCoords << "\n" << zeroGeo;

    std::string geoRB = "noSlip\n";

    for (int rb = 0; rb < QFILES; rb++) {
        *qStreams[rb] << level << "\n" << qFiles[rb].size() << "\n";
        *valueStreams[rb] << geoRB << level << "\n" << qFiles[rb].size() << "\n";
    }
}

void SimulationFileWriter::writeCoordsNeighborsGeo(const int& i, bool binaer, std::shared_ptr<GridBuilder> builder, std::shared_ptr<Transformator> trans)
{
    Grid grid = *builder->getGrid(0, 0).get();
    int index = grid.matrixIndex[i];

    int type = grid.field[index] == SOLID ? 16 : 1;
    real x, y, z;
    grid.transIndexToCoords(index, x, y, z);

    Vertex v = Vertex(x, y, z);
    trans->transformGridToWorld(v);
    double xWorld = v.x;
    double yWorld = v.y;
    double zWorld = v.z;

    if (binaer) 
	{
        xCoordFile.write((char*)&xWorld, sizeof(double));
        yCoordFile.write((char*)&yWorld, sizeof(double));
        zCoordFile.write((char*)&zWorld, sizeof(double));

        xNeighborFile.write((char*)(&grid.neighborIndexX[index] + 1), sizeof(unsigned int));
        yNeighborFile.write((char*)(&grid.neighborIndexY[index] + 1), sizeof(unsigned int));
        zNeighborFile.write((char*)(&grid.neighborIndexZ[index] + 1), sizeof(unsigned int));

        geoVecFile.write((char*)&type, sizeof(unsigned int));
    }
    else 
	{
        xCoordFile << xWorld << " ";
        yCoordFile << yWorld << " ";
        zCoordFile << zWorld << " ";

        xNeighborFile << (grid.neighborIndexX[index] + 1) << " ";
        yNeighborFile << (grid.neighborIndexY[index] + 1) << " ";
        zNeighborFile << (grid.neighborIndexZ[index] + 1) << " ";

        geoVecFile << type << " ";
    }
}

void SimulationFileWriter::writeBoundary(std::vector<real> boundary, int rb)
{
    uint32_t key = *((uint32_t*)&boundary[boundary.size() - 2]);
    int index = (int)boundary[boundary.size() - 1];

    *qStreams[rb] << (index + 1) << " " << key;

    for (int i = 0; i < boundary.size() - 2; i++) {
        *qStreams[rb] << " " << std::fixed << std::setprecision(16) << boundary[i];
    }
    *valueStreams[rb] << (index+1) << " 0 0 0";

    *qStreams[rb] << "\n";
    *valueStreams[rb] << "\n";
}

void SimulationFileWriter::closeFiles()
{
    xCoordFile.close();
    yCoordFile.close();
    zCoordFile.close();
    xNeighborFile.close();
    yNeighborFile.close();
    zNeighborFile.close();
    geoVecFile.close();

    for (int rb = 0; rb < QFILES; rb++) {
        qStreams[rb]->close();
        valueStreams[rb]->close();
    }
}


/*#################################################################################*/
/*-------------------------------static attributes---------------------------------*/
/*---------------------------------------------------------------------------------*/
std::vector<std::tr1::shared_ptr<std::ofstream> > SimulationFileWriter::qStreams;
std::vector<std::tr1::shared_ptr<std::ofstream> > SimulationFileWriter::valueStreams;

std::ofstream SimulationFileWriter::xCoordFile;
std::ofstream SimulationFileWriter::yCoordFile;
std::ofstream SimulationFileWriter::zCoordFile;
std::ofstream SimulationFileWriter::xNeighborFile;
std::ofstream SimulationFileWriter::yNeighborFile;
std::ofstream SimulationFileWriter::zNeighborFile;
std::ofstream SimulationFileWriter::geoVecFile;
std::string SimulationFileWriter::folder;
