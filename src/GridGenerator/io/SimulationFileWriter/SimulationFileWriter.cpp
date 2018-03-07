#define _CRT_SECURE_NO_DEPRECATE
#include "SimulationFileWriter.h"

#include <iostream>
#include <iomanip>
#include <omp.h>
#include "stdint.h"

#include "SimulationFileNames.h"

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/Grid.h>

#include <GridGenerator/grid/GridBuilder/GridBuilder.h>


/*#################################################################################*/
/*---------------------------------public methods----------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::write(std::string folder, SPtr<GridBuilder> builder, FILEFORMAT format)
{
    SimulationFileWriter::folder = folder;

    std::cout << "Write Simulation Files ... \n";
    clock_t  begin = clock();

    write(builder, format);

    clock_t  end = clock();
    real time = real(end - begin) / CLOCKS_PER_SEC;
    std::cout << "time write files: " << time << " sec" << std::endl;
    std::cout << "... finish writing Simulation Files!\n";
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::write(SPtr<GridBuilder> builder, FILEFORMAT format)
{
    const uint numberOfLevel = builder->getNumberOfGridLevels();
    openFiles();
    writeLevel(numberOfLevel);
    for (uint level = 0; level < numberOfLevel; level++)
    {
        writeLevelSize(builder->getNumberOfNodes(level), builder->getQsValues());
        writeCoordFiles(builder, level, format);

        if (level < numberOfLevel - 1)
        {
            writeLevelSizeGridInterface(builder->getNumberOfNodesCF(level), builder->getNumberOfNodesFC(level));
            writeGridInterfaceToFile(builder, level);
        }
        writeBoundaryQsFile(builder->getQsValues());
    }
    closeFiles();
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

    scaleCF_coarse_File.open((path + simulationFileNames::scaleCFC).c_str(), std::ios::out | std::ios::binary);
    scaleCF_fine_File.open((path + simulationFileNames::scaleCFF).c_str(), std::ios::out | std::ios::binary);
    scaleFC_coarse_File.open((path + simulationFileNames::scaleFCC).c_str(), std::ios::out | std::ios::binary);
    scaleFC_fine_File.open((path + simulationFileNames::scaleFCF).c_str(), std::ios::out | std::ios::binary);

    offsetVecCF_File.open((path + simulationFileNames::offsetVecCF).c_str(), std::ios::out | std::ios::binary);
    offsetVecFC_File.open((path + simulationFileNames::offsetVecFC).c_str(), std::ios::out | std::ios::binary);


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
        SPtr<std::ofstream> outQ(new std::ofstream);
        outQ->open(qNames[i].c_str(), std::ios::out | std::ios::binary);
        qStreams.push_back(outQ);

        SPtr<std::ofstream> outV(new std::ofstream);
        outV->open(valueNames[i].c_str(), std::ios::out | std::ios::binary);
        valueStreams.push_back(outV);
    }
}

void SimulationFileWriter::writeLevel(uint numberOfLevels)
{
    const std::string level = std::to_string(numberOfLevels - 1);
    
    xCoordFile << level << "\n";
    yCoordFile << level << "\n";
    zCoordFile << level << "\n";
    xNeighborFile << level << "\n";
    yNeighborFile << level << "\n";
    zNeighborFile << level << "\n";
    geoVecFile << level << "\n";

    scaleCF_coarse_File << level << "\n";
    scaleCF_fine_File << level << "\n";
    scaleFC_coarse_File << level << "\n";
    scaleFC_fine_File << level << "\n";

    offsetVecCF_File << level << "\n";
    offsetVecFC_File << level << "\n";

    const std::string geoRB = "noSlip\n";

    for (int rb = 0; rb < QFILES; rb++) {
        *qStreams[rb] << level << "\n";
        *valueStreams[rb] << geoRB << level << "\n";
    }
}

void SimulationFileWriter::writeLevelSize(uint numberOfNodes, std::vector<std::vector<std::vector<real> > > qFiles)
{
    const std::string zeroIndex = "0 ";
    const std::string zeroGeo = "16 ";

    xCoordFile << numberOfNodes << "\n" << zeroIndex;
    yCoordFile << numberOfNodes << "\n" << zeroIndex;
    zCoordFile << numberOfNodes << "\n" << zeroIndex;
    xNeighborFile << numberOfNodes << "\n" << zeroIndex;
    yNeighborFile << numberOfNodes << "\n" << zeroIndex;
    zNeighborFile << numberOfNodes << "\n" << zeroIndex;
    geoVecFile << numberOfNodes << "\n" << zeroGeo;

    const std::string geoRB = "noSlip\n";

    for (int rb = 0; rb < QFILES; rb++) {
        *qStreams[rb] << qFiles[rb].size() << "\n";
        *valueStreams[rb] << geoRB << qFiles[rb].size() << "\n";
    }
}

void SimulationFileWriter::writeLevelSizeGridInterface(uint sizeCF, uint sizeFC)
{
    scaleCF_coarse_File << sizeCF << "\n";
    scaleCF_fine_File << sizeCF << "\n";
    scaleFC_coarse_File << sizeFC << "\n";
    scaleFC_fine_File << sizeFC << "\n";

    offsetVecCF_File << sizeCF << "\n";
    offsetVecFC_File << sizeFC << "\n";
}

void SimulationFileWriter::writeCoordFiles(SPtr<GridBuilder> builder, uint level, FILEFORMAT format)
{
    for (uint index = 0; index < builder->getNumberOfNodes(level); index++)
        writeCoordsNeighborsGeo(builder, index, level, format);
}

void SimulationFileWriter::writeCoordsNeighborsGeo(SPtr<GridBuilder> builder, int index, uint level, FILEFORMAT format)
{
    SPtr<Grid> grid = builder->getGrid(level);
    if (grid->getSparseIndex(index) == -1)
        return;

    int type = grid->getFieldEntry(index) == SOLID ? 16 : 19;
    real x, y, z;
    grid->transIndexToCoords(index, x, y, z);

    if (format == FILEFORMAT::BINARY)
	{
        xCoordFile.write((char*)&x, sizeof(double));
        yCoordFile.write((char*)&y, sizeof(double));
        zCoordFile.write((char*)&z, sizeof(double));

        xNeighborFile.write((char*)(&grid->getNeighborsX()[index] + 1), sizeof(unsigned int));
        yNeighborFile.write((char*)(&grid->getNeighborsY()[index] + 1), sizeof(unsigned int));
        zNeighborFile.write((char*)(&grid->getNeighborsZ()[index] + 1), sizeof(unsigned int));

        geoVecFile.write((char*)&type, sizeof(unsigned int));
    }
    else 
	{
        xCoordFile << x << " ";
        yCoordFile << y << " ";
        zCoordFile << z << " ";

        xNeighborFile << (grid->getNeighborsX()[index] + 1) << " ";
        yNeighborFile << (grid->getNeighborsY()[index] + 1) << " ";
        zNeighborFile << (grid->getNeighborsZ()[index] + 1) << " ";

        geoVecFile << type << " ";
    }
}

void SimulationFileWriter::writeGridInterfaceToFile(SPtr<GridBuilder> builder, uint level)
{
    const uint numberOfNodesCF = builder->getNumberOfNodesCF(level);
    if(numberOfNodesCF > 0)
    {
        uint* cf_coarse = builder->getCF_coarse(level);
        uint* cf_fine = builder->getCF_fine(level);

        writeGridInterfaceToFile(numberOfNodesCF, scaleCF_coarse_File, cf_coarse, scaleCF_fine_File, cf_fine, offsetVecCF_File);
    }

    const uint numberOfNodesFC = builder->getNumberOfNodesFC(level);
    if (numberOfNodesFC > 0)
    {
        uint* fc_coarse = builder->getFC_coarse(level);
        uint* fc_fine = builder->getFC_fine(level);

        writeGridInterfaceToFile(numberOfNodesFC, scaleFC_coarse_File, fc_coarse, scaleFC_fine_File, fc_fine, offsetVecFC_File);
    }
}

void SimulationFileWriter::writeGridInterfaceToFile(const uint numberOfNodes, std::ofstream& coarseFile, uint* coarse, std::ofstream& fineFile, uint* fine, std::ofstream& offsetFile)
{
    for (uint index = 0; index < numberOfNodes; index++)
    {
        coarseFile << coarse[index] + 1 << " ";
        fineFile << fine[index] + 1 << " ";
        offsetFile << 0 << " " << 0 << " " << 0 << " ";
    }
    coarseFile << "\n";
    fineFile << "\n";
    offsetFile << "\n";
}

void SimulationFileWriter::writeBoundaryQsFile(std::vector<std::vector<std::vector<real> > > qFiles)
{
    for (int rb = 0; rb < QFILES; rb++) {
        for (int index = 0; index < qFiles[rb].size(); index++) {
            writeBoundary(qFiles[rb][index], rb);
        }
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

    scaleCF_coarse_File.close();
    scaleCF_fine_File.close();
    scaleFC_coarse_File.close();
    scaleFC_fine_File.close();

    offsetVecCF_File.close();
    offsetVecFC_File.close();

    for (int rb = 0; rb < QFILES; rb++) {
        qStreams[rb]->close();
        valueStreams[rb]->close();
    }
}


/*#################################################################################*/
/*-------------------------------static attributes---------------------------------*/
/*---------------------------------------------------------------------------------*/
std::vector<SPtr<std::ofstream> > SimulationFileWriter::qStreams;
std::vector<SPtr<std::ofstream> > SimulationFileWriter::valueStreams;

std::ofstream SimulationFileWriter::xCoordFile;
std::ofstream SimulationFileWriter::yCoordFile;
std::ofstream SimulationFileWriter::zCoordFile;
std::ofstream SimulationFileWriter::xNeighborFile;
std::ofstream SimulationFileWriter::yNeighborFile;
std::ofstream SimulationFileWriter::zNeighborFile;
std::ofstream SimulationFileWriter::geoVecFile;
std::string SimulationFileWriter::folder;

std::ofstream SimulationFileWriter::scaleCF_coarse_File;
std::ofstream SimulationFileWriter::scaleCF_fine_File;
std::ofstream SimulationFileWriter::scaleFC_coarse_File;
std::ofstream SimulationFileWriter::scaleFC_fine_File;

std::ofstream SimulationFileWriter::offsetVecCF_File;
std::ofstream SimulationFileWriter::offsetVecFC_File;
