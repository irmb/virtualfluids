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
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include "SimulationFileWriter.h"

#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cmath>

#include "Timer/Timer.h"

#include "grid/NodeValues.h"
#include "grid/Grid.h"
#include "grid/GridBuilder/GridBuilder.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"

#include "io/SimulationFileWriter/SimulationFileNames.h"

#include "utilities/communication.h"

using namespace vf::gpu;

/*#################################################################################*/
/*---------------------------------public methods----------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::write(const std::string& folder, SPtr<GridBuilder> builder, FILEFORMAT format)
{
    SimulationFileWriter::folder = folder;

    VF_LOG_INFO("Start writing simulation files to {}", folder);
    vf::basics::Timer timer;
    timer.start();

    write(builder, format);

    VF_LOG_INFO("    Time writing files: {} sec", timer.getCurrentRuntimeInSeconds());
    VF_LOG_INFO("Done writing simulation Files!");
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void SimulationFileWriter::write(SPtr<GridBuilder> builder, FILEFORMAT format)
{
    const uint numberOfLevel = builder->getNumberOfGridLevels();
    openFiles(builder);
    writeLevel(numberOfLevel);
    //auto qs = createBCVectors(builder->getGrid(0));

    VF_LOG_INFO("   Coordinate and neighbor files:");
    for (uint level = 0; level < numberOfLevel; level++)
    {
        writeNumberNodes(builder, level);
        writeLBMvsSI(builder, level);

        writeLevelSize(builder->getNumberOfNodes(level), format);
        writeCoordFiles(builder, level, format);

        if (level < numberOfLevel - 1)
        {
            writeLevelSizeGridInterface(builder->getNumberOfNodesCF(level), builder->getNumberOfNodesFC(level));
            writeGridInterfaceToFile(builder, level);
        }
    }
    
    VF_LOG_INFO("   Boundary Condition files:");
    writeBoundaryQsFile(builder);
    
    VF_LOG_INFO("    Communication files:");
    writeCommunicationFiles(builder);

    closeFiles();
}

void SimulationFileWriter::openFiles(SPtr<GridBuilder> builder)
{
    std::string path = folder;
    xCoordFile.open((      path + simulationFileNames::coordX).c_str(),      std::ios::out | std::ios::binary);
    yCoordFile.open((      path + simulationFileNames::coordY).c_str(),      std::ios::out | std::ios::binary);
    zCoordFile.open((      path + simulationFileNames::coordZ).c_str(),      std::ios::out | std::ios::binary);
    xNeighborFile.open((   path + simulationFileNames::neighborX).c_str(),   std::ios::out | std::ios::binary);
    yNeighborFile.open((   path + simulationFileNames::neighborY).c_str(),   std::ios::out | std::ios::binary);
    zNeighborFile.open((   path + simulationFileNames::neighborZ).c_str(),   std::ios::out | std::ios::binary);
    wsbNeighborFile.open(( path + simulationFileNames::neighborWSB).c_str(), std::ios::out | std::ios::binary);
    geoVecFile.open((      path + simulationFileNames::geoVec).c_str(),      std::ios::out | std::ios::binary);

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

    if(builder->getCommunicationProcess(CommunicationDirections::MX) != INVALID_INDEX) sendFiles   [CommunicationDirections::MX].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MX)) + "Xs.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PX) != INVALID_INDEX) sendFiles   [CommunicationDirections::PX].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PX)) + "Xs.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::MY) != INVALID_INDEX) sendFiles   [CommunicationDirections::MY].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MY)) + "Ys.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PY) != INVALID_INDEX) sendFiles   [CommunicationDirections::PY].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PY)) + "Ys.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::MZ) != INVALID_INDEX) sendFiles   [CommunicationDirections::MZ].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MZ)) + "Zs.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PZ) != INVALID_INDEX) sendFiles   [CommunicationDirections::PZ].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PZ)) + "Zs.dat").c_str() );

    if(builder->getCommunicationProcess(CommunicationDirections::MX) != INVALID_INDEX) receiveFiles[CommunicationDirections::MX].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MX)) + "Xr.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PX) != INVALID_INDEX) receiveFiles[CommunicationDirections::PX].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PX)) + "Xr.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::MY) != INVALID_INDEX) receiveFiles[CommunicationDirections::MY].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MY)) + "Yr.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PY) != INVALID_INDEX) receiveFiles[CommunicationDirections::PY].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PY)) + "Yr.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::MZ) != INVALID_INDEX) receiveFiles[CommunicationDirections::MZ].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::MZ)) + "Zr.dat").c_str() );
    if(builder->getCommunicationProcess(CommunicationDirections::PZ) != INVALID_INDEX) receiveFiles[CommunicationDirections::PZ].open( (path + std::to_string(builder->getCommunicationProcess(CommunicationDirections::PZ)) + "Zr.dat").c_str() );

    numberNodes_File.open((path + simulationFileNames::numberNodes).c_str(), std::ios::out | std::ios::binary);
    LBMvsSI_File.open((path + simulationFileNames::LBMvsSI).c_str(), std::ios::out | std::ios::binary);
}

void SimulationFileWriter::writeNumberNodes(SPtr<GridBuilder> builder, uint level)
{
    SPtr<Grid> grid = builder->getGrid(level);
    numberNodes_File << level << '\n';

    numberNodes_File << grid->getNumberOfNodesX() << ' ';
    numberNodes_File << grid->getNumberOfNodesY() << ' ';
    numberNodes_File << grid->getNumberOfNodesZ() << ' ';
    numberNodes_File << '\n';
}

void SimulationFileWriter::writeLBMvsSI(SPtr<GridBuilder> builder, uint level)
{
    SPtr<Grid> grid = builder->getGrid(level);

    LBMvsSI_File << grid->getStartX() << ' ';
    LBMvsSI_File << grid->getStartY() << ' ';
    LBMvsSI_File << grid->getStartZ() << ' ';

    LBMvsSI_File << grid->getEndX() << ' ';
    LBMvsSI_File << grid->getEndY() << ' ';
    LBMvsSI_File << grid->getEndZ() << ' ';
    LBMvsSI_File << '\n';
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
    wsbNeighborFile << level << "\n";
    geoVecFile << level << "\n";

    scaleCF_coarse_File << level << "\n";
    scaleCF_fine_File << level << "\n";
    scaleFC_coarse_File << level << "\n";
    scaleFC_fine_File << level << "\n";

    offsetVecCF_File << level << "\n";
    offsetVecFC_File << level << "\n";

    //const std::string geoRB = "noSlip\n";

    //for (int rb = 0; rb < QFILES; rb++) {
    //    *qStreams[rb] << level << "\n";
    //    *valueStreams[rb] << geoRB << level << "\n";
    //}
}

void SimulationFileWriter::writeLevelSize(uint numberOfNodes, FILEFORMAT format)
{
    const std::string zeroIndex = "0 ";
    const std::string zeroGeo = "16 ";

    if (format == FILEFORMAT::BINARY)
    {
        //const uint zeroIndex = 0;
        //const uint zeroGeo   = 16;

        //xCoordFile.write((char*)&zeroIndex, sizeof(double));
        //yCoordFile.write((char*)&zeroIndex, sizeof(double));
        //zCoordFile.write((char*)&zeroIndex, sizeof(double));

        //// + 1 for numbering shift between GridGenerator and VF_GPU
        //xNeighborFile.write((char*)(&zeroIndex), sizeof(unsigned int));
        //yNeighborFile.write((char*)(&zeroIndex), sizeof(unsigned int));
        //zNeighborFile.write((char*)(&zeroIndex), sizeof(unsigned int));
        //wsbNeighborFile.write((char*)(&zeroIndex), sizeof(unsigned int));

        //geoVecFile.write((char*)&zeroGeo, sizeof(unsigned int));
        
        xCoordFile      << numberOfNodes << "\n" << zeroIndex;
        yCoordFile      << numberOfNodes << "\n" << zeroIndex;
        zCoordFile      << numberOfNodes << "\n" << zeroIndex;
        xNeighborFile   << numberOfNodes << "\n" << zeroIndex;
        yNeighborFile   << numberOfNodes << "\n" << zeroIndex;
        zNeighborFile   << numberOfNodes << "\n" << zeroIndex;
        wsbNeighborFile << numberOfNodes << "\n" << zeroIndex;
        geoVecFile      << numberOfNodes << "\n" << zeroGeo  ;
    }
    else 
    {
        xCoordFile      << numberOfNodes << "\n" << zeroIndex << "\n";
        yCoordFile      << numberOfNodes << "\n" << zeroIndex << "\n";
        zCoordFile      << numberOfNodes << "\n" << zeroIndex << "\n";
        xNeighborFile   << numberOfNodes << "\n" << zeroIndex << "\n";
        yNeighborFile   << numberOfNodes << "\n" << zeroIndex << "\n";
        zNeighborFile   << numberOfNodes << "\n" << zeroIndex << "\n";
        wsbNeighborFile << numberOfNodes << "\n" << zeroIndex << "\n";
        geoVecFile      << numberOfNodes << "\n" << zeroGeo   << "\n";
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
    for (uint index = 0; index < builder->getGrid(level)->getSize(); index++){
        writeCoordsNeighborsGeo(builder, index, level, format);
    }

    xCoordFile << "\n";
    yCoordFile << "\n";
    zCoordFile << "\n";

    xNeighborFile << "\n";
    yNeighborFile << "\n";
    zNeighborFile << "\n";
    wsbNeighborFile << "\n";

    geoVecFile << "\n";
}

void SimulationFileWriter::writeCoordsNeighborsGeo(SPtr<GridBuilder> builder, int index, uint level, FILEFORMAT format)
{
    SPtr<Grid> grid = builder->getGrid(level);
    if (grid->getSparseIndex(index) == -1)
        return;

    // Lenz: in the GPU code all nodes that perform collisions have to be fluid = 19
    bool isStopper = grid->getFieldEntry(index) == STOPPER_SOLID                || 
                     grid->getFieldEntry(index) == STOPPER_OUT_OF_GRID          || 
                     grid->getFieldEntry(index) == STOPPER_OUT_OF_GRID_BOUNDARY || 
                     grid->getFieldEntry(index) == STOPPER_COARSE_UNDER_FINE;
    int type = !isStopper ? 19 : 16;

    // old code from Soeren P.
    //int type = grid->getFieldEntry(index) == FLUID ? 19 : 16;

    real x, y, z;
    grid->transIndexToCoords(index, x, y, z);

    if (format == FILEFORMAT::BINARY)
    {
        double tmpX = (double)x;
        double tmpY = (double)y;
        double tmpZ = (double)z;

        xCoordFile.write((char*)&tmpX, sizeof(double));
        yCoordFile.write((char*)&tmpY, sizeof(double));
        zCoordFile.write((char*)&tmpZ, sizeof(double));

        // + 1 for numbering shift between GridGenerator and VF_GPU
        int tmpNeighborX        = grid->getNeighborsX()[index]        + 1;
        int tmpNeighborY        = grid->getNeighborsY()[index]        + 1;
        int tmpNeighborZ        = grid->getNeighborsZ()[index]        + 1;
        int tmpNeighborNegative = grid->getNeighborsNegative()[index] + 1;

        xNeighborFile.write  ((char*)(&tmpNeighborX       ), sizeof(unsigned int));
        yNeighborFile.write  ((char*)(&tmpNeighborY       ), sizeof(unsigned int));
        zNeighborFile.write  ((char*)(&tmpNeighborZ       ), sizeof(unsigned int));
        wsbNeighborFile.write((char*)(&tmpNeighborNegative), sizeof(unsigned int));

        geoVecFile.write((char*)&type, sizeof(unsigned int));
    }
    else 
    {
        xCoordFile << x << "\n";
        yCoordFile << y << "\n";
        zCoordFile << z << "\n";

        // + 1 for numbering shift between GridGenerator and VF_GPU
        xNeighborFile << (grid->getNeighborsX()[index] + 1) << "\n";
        yNeighborFile << (grid->getNeighborsY()[index] + 1) << "\n";
        zNeighborFile << (grid->getNeighborsZ()[index] + 1) << "\n";
        wsbNeighborFile << (grid->getNeighborsNegative()[index] + 1) << "\n";

        geoVecFile << type << "\n";
    }
}

void SimulationFileWriter::writeGridInterfaceToFile(SPtr<GridBuilder> builder, uint level)
{
    const uint numberOfNodesCF = builder->getNumberOfNodesCF(level);
    const uint numberOfNodesFC = builder->getNumberOfNodesFC(level);

    {
        uint* cf_coarse = new uint[numberOfNodesCF];
        uint* cf_fine = new uint[numberOfNodesCF];
        uint* fc_coarse = new uint[numberOfNodesFC];
        uint* fc_fine = new uint[numberOfNodesFC];

        builder->getGridInterfaceIndices(cf_coarse, cf_fine, fc_coarse, fc_fine, level);

        if (numberOfNodesCF > 0)
        {
            writeGridInterfaceToFile(numberOfNodesCF, scaleCF_coarse_File, cf_coarse, scaleCF_fine_File, cf_fine);
        }

        if (numberOfNodesFC > 0)
        {
            writeGridInterfaceToFile(numberOfNodesFC, scaleFC_coarse_File, fc_coarse, scaleFC_fine_File, fc_fine);
        }

        delete [] cf_coarse;
        delete [] cf_fine;
        delete [] fc_coarse;
        delete [] fc_fine;
    }

    {
        real* cf_offset_X = new real[numberOfNodesCF];
        real* cf_offset_Y = new real[numberOfNodesCF];
        real* cf_offset_Z = new real[numberOfNodesCF];

        real* fc_offset_X = new real[numberOfNodesFC];
        real* fc_offset_Y = new real[numberOfNodesFC];
        real* fc_offset_Z = new real[numberOfNodesFC];

        builder->getOffsetCF(cf_offset_X, cf_offset_Y, cf_offset_Z, level);
        builder->getOffsetFC(fc_offset_X, fc_offset_Y, fc_offset_Z, level);

        if (numberOfNodesCF > 0)
        {
            writeGridInterfaceOffsetToFile(numberOfNodesCF, offsetVecCF_File, cf_offset_X, cf_offset_Y, cf_offset_Z);
        }

        if (numberOfNodesFC > 0)
        {
            writeGridInterfaceOffsetToFile(numberOfNodesFC, offsetVecFC_File, fc_offset_X, fc_offset_Y, fc_offset_Z);
        }

        delete[] cf_offset_X;
        delete[] cf_offset_Y;
        delete[] cf_offset_Z;

        delete[] fc_offset_X;
        delete[] fc_offset_Y;
        delete[] fc_offset_Z;
    }
}

void SimulationFileWriter::writeGridInterfaceToFile(uint numberOfNodes, std::ofstream &coarseFile, uint *coarse,
                                                    std::ofstream &fineFile, uint *fine)
{
    for (uint index = 0; index < numberOfNodes; index++) {
        coarseFile << coarse[index] << " \n";
        fineFile << fine[index] << " \n";
    }
    coarseFile << "\n";
    fineFile << "\n";
}

void SimulationFileWriter::writeGridInterfaceOffsetToFile(uint numberOfNodes, std::ofstream &offsetFile, real *offset_X,
                                                          real *offset_Y, real *offset_Z)
{
    for (uint index = 0; index < numberOfNodes; index++) {
        offsetFile << offset_X[index] << " " << offset_Y[index] << " " << offset_Z[index] << " \n";
    }
    offsetFile << "\n";
}

/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
std::vector<std::vector<std::vector<real> > > SimulationFileWriter::createBCVectors(SPtr<Grid> grid)
{
    std::vector<std::vector<std::vector<real> > > qs;
    qs.resize(QFILES);
    for (uint i = 0; i < grid->getSize(); i++)
    {
        real x, y, z;
        grid->transIndexToCoords(i, x, y, z);

        if (grid->getFieldEntry(i) == BC_SOLID) addShortQsToVector(i, qs, grid); //addQsToVector(i, qs, grid);
        //if (x == 0 && y < grid->getNumberOfNodesY() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 0, -1, INLETQS, qs, grid);
        //if (x == grid->getNumberOfNodesX() - 2 && y < grid->getNumberOfNodesY() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 0, 1, OUTLETQS, qs, grid);

        //if (z == grid->getNumberOfNodesZ() - 2 && x < grid->getNumberOfNodesX() - 1 && y < grid->getNumberOfNodesY() - 1) fillRBForNode(i, 2, 1, TOPQS, qs, grid);
        //if (z == 0 && x < grid->getNumberOfNodesX() - 1 && y < grid->getNumberOfNodesY() - 1) fillRBForNode(i, 2, -1, BOTTOMQS, qs, grid);

        //if (y == 0 && x < grid->getNumberOfNodesX() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 1, -1, FRONTQS, qs, grid);
        //if (y == grid->getNumberOfNodesY() - 2 && x < grid->getNumberOfNodesX() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 1, 1, BACKQS, qs, grid);
    }
    return qs;
}

void SimulationFileWriter::addShortQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        /*int qIndex = i * grid->getSize() + grid->getSparseIndex(index);
        real q = grid->getDistribution()[qIndex];*/
        real q = grid->getQValue(index, i);
        if (q > 0) {
            //printf("Q%d (old:%d, new:%d), : %2.8f \n", i, coordsVec[index].matrixIndex, index, grid.d.f[i * grid.size + coordsVec[index].matrixIndex]);
            qKey += (uint32_t)pow(2, 26 - i);
            qNode.push_back(q);
        }
    }
    if (qKey > 0) {
        float transportKey = *((float*)&qKey);
        qNode.push_back((real)transportKey);
        qNode.push_back((real)index);
        qs[GEOMQS].push_back(qNode);
    }
    qNode.clear();
}

void SimulationFileWriter::addQsToVector(int index, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid)
{
    std::vector<real> qNode;
    qNode.push_back((real)index);

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        //int qIndex = i * grid->getSize() + grid->getSparseIndex(index);
        //real q = grid->getDistribution()[qIndex];
        real q = grid->getQValue(index, i);
        qNode.push_back(q);
        if (q > 0)
            printf("Q= %f; Index = %d \n", q, index);
            //qNode.push_back(q);
  //      else
  //          qNode.push_back(-1);
    }
    qs[GEOMQS].push_back(qNode);
    qNode.clear();
}

void SimulationFileWriter::fillRBForNode(int index, int direction, int directionSign, int rb, std::vector<std::vector<std::vector<real> > > &qs, SPtr<Grid> grid)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        if (grid->getDirection()[i * DIMENSION + direction] != directionSign)
            continue;

        qKey += (uint32_t)pow(2, 26 - i);
        qNode.push_back(0.5f);
    }
    if (qKey > 0) {
        float transportKey = *((float*)&qKey);
        qNode.push_back((real)transportKey);
        qNode.push_back((real)index);
        qs[rb].push_back(qNode);
    }
    qNode.clear();
}


void SimulationFileWriter::writeBoundaryQsFile(SPtr<GridBuilder> builder)
{
    // old code of Soeren P.
  //  for (int rb = 0; rb < QFILES; rb++) {
  //      for (int index = 0; index < qFiles[rb].size(); index++) {
  //          //writeBoundary(qFiles[rb][index], rb);
        //    writeBoundaryShort(qFiles[rb][index], rb);
        //}
  //  }

    SideType sides[] = {SideType::MX, SideType::PX, SideType::PZ, SideType::MZ, SideType::MY, SideType::PY, SideType::GEOMETRY};

    for (uint side = 0; side < QFILES; side++) {

        for (uint level = 0; level < builder->getNumberOfGridLevels(); level++) {
            
            auto bc = builder->getBoundaryCondition( sides[side], level );

            if( level == 0 ){
            
                if     ( !bc )                          *valueStreams[side] << "noSlip\n";
                else if( bc->getType() == BC_PRESSURE ) *valueStreams[side] << "pressure\n";
                else if( bc->getType() == BC_VELOCITY ) *valueStreams[side] << "velocity\n";
                else if( bc->getType() == BC_SOLID    ) *valueStreams[side] << "noSlip\n";
                else if( bc->getType() == BC_SLIP     ) *valueStreams[side] << "slip\n";
                else if( bc->getType() == BC_OUTFLOW  ) *valueStreams[side] << "outflow\n";

                *valueStreams[side] << builder->getNumberOfGridLevels() - 1 << "\n";
                *qStreams[side]     << builder->getNumberOfGridLevels() - 1 << "\n";
            }

            if( !bc ){
                *valueStreams[side] << "0\n";
                *qStreams[side]     << "0\n";
            }
            else{
                writeBoundaryShort( builder->getGrid(level), bc, side );
            }
        }
    }
}

void SimulationFileWriter::writeBoundary(std::vector<real> boundary, int rb)
{
    int index = (int)boundary[0];

    *qStreams[rb] << (index + 1);

    for (std::size_t i = 1; i < boundary.size(); i++) {
        *qStreams[rb] << " " << std::fixed << std::setprecision(16) << boundary[i];
    }
    *valueStreams[rb] << (index+1) << " 0 0 0";

    *qStreams[rb] << "\n";
    *valueStreams[rb] << "\n";
}

void SimulationFileWriter::writeBoundaryShort(std::vector<real> boundary, int rb)
{
    uint32_t key = *((uint32_t*)&boundary[boundary.size() - 2]);
    int index = (int)boundary[boundary.size() - 1];

    *qStreams[rb] << (index + 1) << " " << key;

    for (std::size_t i = 0; i < boundary.size() - 2; i++) {
        *qStreams[rb] << " " << std::fixed << std::setprecision(16) << boundary[i];
    }
    *valueStreams[rb] << (index + 1) << " 0 0 0";

    *qStreams[rb] << "\n";
    *valueStreams[rb] << "\n";
}

void SimulationFileWriter::writeBoundaryShort(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, uint side)
{
    uint numberOfBoundaryNodes = (uint)boundaryCondition->indices.size();

    *valueStreams[side] << numberOfBoundaryNodes << "\n";
    *qStreams[side]     << numberOfBoundaryNodes << "\n";

    for( uint index = 0; index < numberOfBoundaryNodes; index++ ){
    
        // + 1 for numbering shift between GridGenerator and VF_GPU
        *valueStreams[side] << grid->getSparseIndex( boundaryCondition->indices[index] ) + 1 << " ";
        *qStreams[side]     << grid->getSparseIndex( boundaryCondition->indices[index] ) + 1 << " ";

        {
            uint32_t key = 0;

            for (int dir = 26; dir >= 0; dir--)
            {
                real q = boundaryCondition->getQ(index,dir);
                if (q > 0) {
                    key += (uint32_t)pow(2, 26 - dir);
                }
            }

            *qStreams[side] << key << " ";

            for (int dir = 26; dir >= 0; dir--)
            {
                real q = boundaryCondition->getQ(index,dir);
                if (q > 0) {
                    *qStreams[side] << std::fixed << std::setprecision(16) << q << " ";
                }
            }

            *qStreams[side] << "\n";
        }

        if( boundaryCondition->getType() == BC_PRESSURE )
        {
            auto bcPressure = dynamic_cast< PressureBoundaryCondition* >( boundaryCondition.get() );
            if(bcPressure != nullptr) {
                *valueStreams[side] << bcPressure->getRho() << " ";
                // + 1 for numbering shift between GridGenerator and VF_GPU
                *valueStreams[side] << grid->getSparseIndex( bcPressure->neighborIndices[index] ) + 1 << " ";
            }
        }

        if( boundaryCondition->getType() == BC_VELOCITY )
        {
            auto bcVelocity = dynamic_cast< VelocityBoundaryCondition* >( boundaryCondition.get() );
            if(bcVelocity != nullptr) {
                *valueStreams[side] << bcVelocity->getVx(index) << " ";
                *valueStreams[side] << bcVelocity->getVy(index) << " ";
                *valueStreams[side] << bcVelocity->getVz(index) << " ";
            }
        }

        if( boundaryCondition->getType() == BC_SOLID )
        {
            auto bcGeometry = dynamic_cast< GeometryBoundaryCondition* >( boundaryCondition.get() );
            if(bcGeometry != nullptr) {
                *valueStreams[side] << bcGeometry->getVx(index) << " ";
                *valueStreams[side] << bcGeometry->getVy(index) << " ";
                *valueStreams[side] << bcGeometry->getVz(index) << " ";
            }
        }

        *valueStreams[side] << "\n";
    }
    
}

void SimulationFileWriter::writeCommunicationFiles(SPtr<GridBuilder> builder)
{
    const uint numberOfLevel = builder->getNumberOfGridLevels();
    
    for( uint direction = 0; direction < 6; direction++ ){
        
        if(builder->getCommunicationProcess(direction) == INVALID_INDEX) continue;

        sendFiles[direction] << "processor\n";
        sendFiles[direction] << builder->getNumberOfGridLevels() - 1 << "\n";

        receiveFiles[direction] << "processor\n";
        receiveFiles[direction] << builder->getNumberOfGridLevels() - 1 << "\n";

        for (uint level = 0; level < numberOfLevel; level++){
        
            uint numberOfSendNodes    = builder->getGrid(level)->getNumberOfSendNodes(direction);
            uint numberOfReceiveNodes = builder->getGrid(level)->getNumberOfReceiveNodes(direction);
        
            sendFiles[direction] <<    numberOfSendNodes    << "\n";
            receiveFiles[direction] << numberOfReceiveNodes << "\n";

            for( uint index = 0; index < numberOfSendNodes; index++ )
                // + 1 for numbering shift between GridGenerator and VF_GPU
                sendFiles[direction]    << builder->getGrid(level)->getSparseIndex( builder->getGrid(level)->getSendIndex   (direction, index) ) + 1 << "\n";

            for( uint index = 0; index < numberOfReceiveNodes; index++ )
                // + 1 for numbering shift between GridGenerator and VF_GPU
                receiveFiles[direction] << builder->getGrid(level)->getSparseIndex( builder->getGrid(level)->getReceiveIndex(direction, index) ) + 1 << "\n";

        }
    }
}

void SimulationFileWriter::closeFiles()
{
    xCoordFile.close();
    yCoordFile.close();
    zCoordFile.close();
    xNeighborFile.close();
    yNeighborFile.close();
    zNeighborFile.close();
    wsbNeighborFile.close();
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

    sendFiles[CommunicationDirections::MX].close();
    sendFiles[CommunicationDirections::PX].close();
    sendFiles[CommunicationDirections::MY].close();
    sendFiles[CommunicationDirections::PY].close();
    sendFiles[CommunicationDirections::MZ].close();
    sendFiles[CommunicationDirections::PZ].close();
    
    receiveFiles[CommunicationDirections::MX].close();
    receiveFiles[CommunicationDirections::PX].close();
    receiveFiles[CommunicationDirections::MY].close();
    receiveFiles[CommunicationDirections::PY].close();
    receiveFiles[CommunicationDirections::MZ].close();
    receiveFiles[CommunicationDirections::PZ].close();

    numberNodes_File.close();
    LBMvsSI_File.close();
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
std::ofstream SimulationFileWriter::wsbNeighborFile;
std::ofstream SimulationFileWriter::geoVecFile;
std::string SimulationFileWriter::folder;

std::ofstream SimulationFileWriter::scaleCF_coarse_File;
std::ofstream SimulationFileWriter::scaleCF_fine_File;
std::ofstream SimulationFileWriter::scaleFC_coarse_File;
std::ofstream SimulationFileWriter::scaleFC_fine_File;

std::ofstream SimulationFileWriter::offsetVecCF_File;
std::ofstream SimulationFileWriter::offsetVecFC_File;

std::ofstream SimulationFileWriter::numberNodes_File;
std::ofstream SimulationFileWriter::LBMvsSI_File;

std::array<std::ofstream, 6> SimulationFileWriter::sendFiles;
std::array<std::ofstream, 6> SimulationFileWriter::receiveFiles;
//! \}
