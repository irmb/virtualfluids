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
//! \addtogroup cpu_SimulationObservers SimulationObservers
//! \ingroup cpu_core core
//! \{
//! \author Alena Karanchuk
//=======================================================================================
#include "MPIIOSimulationObserver.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "MPIIODataStructures.h"
#include "MemoryUtil.h"
#include "UbFileInputASCII.h"
#include "UbFileOutputASCII.h"
#include "UbLogger.h"
#include "UbScheduler.h"

using namespace MPIIODataStructures;

MPIIOSimulationObserver::MPIIOSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                   std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), comm(comm)
{
    UbSystem::makeDirectory(path + "/mpi_io_cp");

    //-------------------------   define MPI types  ---------------------------------

    MPI_Datatype typesGP[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
    int blocksGP[3]         = { 34, 6, 5 };
    MPI_Aint offsetsGP[3], lbGP, extentGP;

    offsetsGP[0] = 0;
    MPI_Type_get_extent(MPI_DOUBLE, &lbGP, &extentGP);
    offsetsGP[1] = blocksGP[0] * extentGP;

    MPI_Type_get_extent(MPI_INT, &lbGP, &extentGP);
    offsetsGP[2] = offsetsGP[1] + blocksGP[1] * extentGP;

    MPI_Type_create_struct(3, blocksGP, offsetsGP, typesGP, &gridParamType);
    MPI_Type_commit(&gridParamType);

    //-----------------------------------------------------------------------

    MPI_Datatype typesBlock[2] = { MPI_INT, MPI_CHAR };
    int blocksBlock[2]         = { 13, 1 };
    MPI_Aint offsetsBlock[2], lbBlock, extentBlock;

    offsetsBlock[0] = 0;
    MPI_Type_get_extent(MPI_INT, &lbBlock, &extentBlock);
    offsetsBlock[1] = blocksBlock[0] * extentBlock;

    MPI_Type_create_struct(2, blocksBlock, offsetsBlock, typesBlock, &block3dType);
    MPI_Type_commit(&block3dType);

    //---------------------------------------

    MPI_Type_contiguous(7, MPI_INT, &dataSetParamType);
    MPI_Type_commit(&dataSetParamType);

    //-----------------------------------------------------------------------

    MPI_Datatype typesBC[3] = { MPI_LONG_LONG_INT, MPI_FLOAT, MPI_CHAR };
    int blocksBC[3]         = { 5, 34, 1 };
    MPI_Aint offsetsBC[3], lbBC, extentBC;

    offsetsBC[0] = 0;
    MPI_Type_get_extent(MPI_LONG_LONG_INT, &lbBC, &extentBC);
    offsetsBC[1] = blocksBC[0] * extentBC;

    MPI_Type_get_extent(MPI_FLOAT, &lbBC, &extentBC);
    offsetsBC[2] = offsetsBC[1] + blocksBC[1] * extentBC;

    MPI_Type_create_struct(3, blocksBC, offsetsBC, typesBC, &boundCondType);
    MPI_Type_commit(&boundCondType);

    //---------------------------------------

    MPI_Type_contiguous(9, MPI_CHAR, &arrayPresenceType);
    MPI_Type_commit(&arrayPresenceType);
}

MPIIOSimulationObserver::~MPIIOSimulationObserver()
{
    MPI_Type_free(&gridParamType);
    MPI_Type_free(&block3dType);
    MPI_Type_free(&dataSetParamType);
    MPI_Type_free(&boundCondType);
    MPI_Type_free(&arrayPresenceType);
}

void MPIIOSimulationObserver::writeBlocks(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
    size = 1;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIOSimulationObserver::writeBlocksToFile start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    int blocksCount  = 0; // quantity of all the blocks in the grid, max 2147483648 blocks!
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();

    std::vector<SPtr<Block3D>> blocksVector[25]; // max 25 levels
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        // grid->getBlocks(level, rank, blockVector[level]);
        grid->getBlocks(level, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    GridParam *gridParameters      = new GridParam;
    gridParameters->trafoParams[0] = grid->getCoordinateTransformator()->Tx1;
    gridParameters->trafoParams[1] = grid->getCoordinateTransformator()->Tx2;
    gridParameters->trafoParams[2] = grid->getCoordinateTransformator()->Tx3;
    gridParameters->trafoParams[3] = grid->getCoordinateTransformator()->Sx1;
    gridParameters->trafoParams[4] = grid->getCoordinateTransformator()->Sx2;
    gridParameters->trafoParams[5] = grid->getCoordinateTransformator()->Sx3;
    gridParameters->trafoParams[6] = grid->getCoordinateTransformator()->alpha;
    gridParameters->trafoParams[7] = grid->getCoordinateTransformator()->beta;
    gridParameters->trafoParams[8] = grid->getCoordinateTransformator()->gamma;

    gridParameters->trafoParams[9]  = grid->getCoordinateTransformator()->toX1factorX1;
    gridParameters->trafoParams[10] = grid->getCoordinateTransformator()->toX1factorX2;
    gridParameters->trafoParams[11] = grid->getCoordinateTransformator()->toX1factorX3;
    gridParameters->trafoParams[12] = grid->getCoordinateTransformator()->toX1delta;
    gridParameters->trafoParams[13] = grid->getCoordinateTransformator()->toX2factorX1;
    gridParameters->trafoParams[14] = grid->getCoordinateTransformator()->toX2factorX2;
    gridParameters->trafoParams[15] = grid->getCoordinateTransformator()->toX2factorX3;
    gridParameters->trafoParams[16] = grid->getCoordinateTransformator()->toX2delta;
    gridParameters->trafoParams[17] = grid->getCoordinateTransformator()->toX3factorX1;
    gridParameters->trafoParams[18] = grid->getCoordinateTransformator()->toX3factorX2;
    gridParameters->trafoParams[19] = grid->getCoordinateTransformator()->toX3factorX3;
    gridParameters->trafoParams[20] = grid->getCoordinateTransformator()->toX3delta;

    gridParameters->trafoParams[21] = grid->getCoordinateTransformator()->fromX1factorX1;
    gridParameters->trafoParams[22] = grid->getCoordinateTransformator()->fromX1factorX2;
    gridParameters->trafoParams[23] = grid->getCoordinateTransformator()->fromX1factorX3;
    gridParameters->trafoParams[24] = grid->getCoordinateTransformator()->fromX1delta;
    gridParameters->trafoParams[25] = grid->getCoordinateTransformator()->fromX2factorX1;
    gridParameters->trafoParams[26] = grid->getCoordinateTransformator()->fromX2factorX2;
    gridParameters->trafoParams[27] = grid->getCoordinateTransformator()->fromX2factorX3;
    gridParameters->trafoParams[28] = grid->getCoordinateTransformator()->fromX2delta;
    gridParameters->trafoParams[29] = grid->getCoordinateTransformator()->fromX3factorX1;
    gridParameters->trafoParams[30] = grid->getCoordinateTransformator()->fromX3factorX2;
    gridParameters->trafoParams[31] = grid->getCoordinateTransformator()->fromX3factorX3;
    gridParameters->trafoParams[32] = grid->getCoordinateTransformator()->fromX3delta;

    gridParameters->active         = grid->getCoordinateTransformator()->active;
    gridParameters->transformation = grid->getCoordinateTransformator()->transformation;

    gridParameters->deltaX     = grid->getDeltaX(minInitLevel);
    UbTupleInt3 blocknx        = grid->getBlockNX();
    gridParameters->blockNx1   = val<1>(blocknx);
    gridParameters->blockNx2   = val<2>(blocknx);
    gridParameters->blockNx3   = val<3>(blocknx);
    gridParameters->nx1        = grid->getNX1();
    gridParameters->nx2        = grid->getNX2();
    gridParameters->nx3        = grid->getNX3();
    gridParameters->periodicX1 = grid->isPeriodicX1();
    gridParameters->periodicX2 = grid->isPeriodicX2();
    gridParameters->periodicX3 = grid->isPeriodicX3();

    //----------------------------------------------------------------------

    Block3d *block3dArray = new Block3d[blocksCount];
    int ic                = 0;
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //    all the blocks of the current level
        {
            // save data describing the block
            block3dArray[ic].x1                  = block->getX1();
            block3dArray[ic].x2                  = block->getX2();
            block3dArray[ic].x3                  = block->getX3();
            block3dArray[ic].bundle              = block->getBundle();
            block3dArray[ic].rank                = block->getRank();
            block3dArray[ic].lrank               = block->getLocalRank();
            block3dArray[ic].part                = block->getPart();
            block3dArray[ic].globalID            = block->getGlobalID();
            block3dArray[ic].localID             = block->getLocalID();
            block3dArray[ic].level               = block->getLevel();
            block3dArray[ic].interpolationFlagCF = block->getCollectionOfInterpolationFlagCF();
            block3dArray[ic].interpolationFlagFC = block->getCollectionOfInterpolationFlagFC();
            block3dArray[ic].counter             = block->getMaxGlobalID();
            block3dArray[ic].active              = block->isActive();

            ic++;
        }
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIOSimulationObserver::writeBlocksToFile start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    MPI_File file_handler;
    //   MPI_Info info = MPI_INFO_NULL;

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                           &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    real start {0.};
    real finish {0.};
    MPI_Offset write_offset = (MPI_Offset)(size * sizeof(int));

    if (comm->isRoot()) {
        start = MPI_Wtime();

        // each process writes the quantity of it's blocks
        MPI_File_write_at(file_handler, 0 /*rank*sizeof(int)*/, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
        // each process writes parameters of the grid
        MPI_File_write_at(file_handler, write_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
        // each process writes it's blocks
        MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(GridParam)), &block3dArray[0], blocksCount,
                          block3dType, MPI_STATUS_IGNORE);
        // MPI_File_sync(file_handler);
    }
    MPI_File_close(&file_handler);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOSimulationObserver::writeBlocksToFile time: " << finish - start << " s");
    }

    delete[] block3dArray;
    delete gridParameters;
}

void MPIIOSimulationObserver::readBlocks(int step)
{
    int rank;
    //   int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //   size = 1;
    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIOSimulationObserver::readBlocks start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    MPI_File_read_at(file_handler, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    Block3d *block3dArray = new Block3d[blocksCount];

    GridParam *gridParameters = new GridParam;

    // calculate the read offset
    MPI_Offset read_offset = (MPI_Offset)(sizeof(int));

    // read parameters of the grid
    MPI_File_read_at(file_handler, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
    // read all the blocks
    if (comm->isRoot())
        MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(GridParam)), &block3dArray[0], blocksCount,
                         block3dType, MPI_STATUS_IGNORE);

    MPI_Bcast(block3dArray, blocksCount, block3dType, comm->getRoot(), MPI_COMM_WORLD);

    MPI_File_close(&file_handler);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOSimulationObserver::readBlocks time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOSimulationObserver::readBlocks start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // clear the grid
    grid->deleteBlocks();

    // restore the grid
    SPtr<CoordinateTransformation3D> trafo(new CoordinateTransformation3D());
    trafo->Tx1   = gridParameters->trafoParams[0];
    trafo->Tx2   = gridParameters->trafoParams[1];
    trafo->Tx3   = gridParameters->trafoParams[2];
    trafo->Sx1   = gridParameters->trafoParams[3];
    trafo->Sx2   = gridParameters->trafoParams[4];
    trafo->Sx3   = gridParameters->trafoParams[5];
    trafo->alpha = gridParameters->trafoParams[6];
    trafo->beta  = gridParameters->trafoParams[7];
    trafo->gamma = gridParameters->trafoParams[8];

    trafo->toX1factorX1 = gridParameters->trafoParams[9];
    trafo->toX1factorX2 = gridParameters->trafoParams[10];
    trafo->toX1factorX3 = gridParameters->trafoParams[11];
    trafo->toX1delta    = gridParameters->trafoParams[12];
    trafo->toX2factorX1 = gridParameters->trafoParams[13];
    trafo->toX2factorX2 = gridParameters->trafoParams[14];
    trafo->toX2factorX3 = gridParameters->trafoParams[15];
    trafo->toX2delta    = gridParameters->trafoParams[16];
    trafo->toX3factorX1 = gridParameters->trafoParams[17];
    trafo->toX3factorX2 = gridParameters->trafoParams[18];
    trafo->toX3factorX3 = gridParameters->trafoParams[19];
    trafo->toX3delta    = gridParameters->trafoParams[20];

    trafo->fromX1factorX1 = gridParameters->trafoParams[21];
    trafo->fromX1factorX2 = gridParameters->trafoParams[22];
    trafo->fromX1factorX3 = gridParameters->trafoParams[23];
    trafo->fromX1delta    = gridParameters->trafoParams[24];
    trafo->fromX2factorX1 = gridParameters->trafoParams[25];
    trafo->fromX2factorX2 = gridParameters->trafoParams[26];
    trafo->fromX2factorX3 = gridParameters->trafoParams[27];
    trafo->fromX2delta    = gridParameters->trafoParams[28];
    trafo->fromX3factorX1 = gridParameters->trafoParams[29];
    trafo->fromX3factorX2 = gridParameters->trafoParams[30];
    trafo->fromX3factorX3 = gridParameters->trafoParams[31];
    trafo->fromX3delta    = gridParameters->trafoParams[32];

    trafo->active         = gridParameters->active;
    trafo->transformation = gridParameters->transformation;

    grid->setCoordinateTransformator(trafo);

    grid->setDeltaX(gridParameters->deltaX);
    grid->setBlockNX(gridParameters->blockNx1, gridParameters->blockNx2, gridParameters->blockNx3);
    grid->setNX1(gridParameters->nx1);
    grid->setNX2(gridParameters->nx2);
    grid->setNX3(gridParameters->nx3);
    grid->setPeriodicX1(gridParameters->periodicX1);
    grid->setPeriodicX2(gridParameters->periodicX2);
    grid->setPeriodicX3(gridParameters->periodicX3);

    // regenerate blocks
    for (int n = 0; n < blocksCount; n++) {
        SPtr<Block3D> block(
            new Block3D(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, block3dArray[n].level));
        block->setActive(block3dArray[n].active);
        block->setBundle(block3dArray[n].bundle);
        block->setRank(block3dArray[n].rank);
        block->setLocalRank(block3dArray[n].lrank);
        block->setGlobalID(block3dArray[n].globalID);
        block->setLocalID(block3dArray[n].localID);
        block->setPart(block3dArray[n].part);
        block->setLevel(block3dArray[n].level);
        block->setCollectionOfInterpolationFlagCF(block3dArray[n].interpolationFlagCF);
        block->setCollectionOfInterpolationFlagFC(block3dArray[n].interpolationFlagFC);

        grid->addBlock(block);
    }

    delete gridParameters;
    delete[] block3dArray;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIOSimulationObserver::readBlocks end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}

void MPIIOSimulationObserver::clearAllFiles(int step)
{
    MPI_File file_handler;
    MPI_Info info       = MPI_INFO_NULL;
    MPI_Offset new_size = 0;

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

    std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
    int rc1 = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                            &file_handler);
    if (rc1 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename1);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename21 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetF.bin";
    int rc21 = MPI_File_open(MPI_COMM_WORLD, filename21.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc21 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename21);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename22 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
    int rc22 = MPI_File_open(MPI_COMM_WORLD, filename22.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc22 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename22);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename23 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH2.bin";
    int rc23 = MPI_File_open(MPI_COMM_WORLD, filename23.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc23 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename23);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename3 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
    int rc3 = MPI_File_open(MPI_COMM_WORLD, filename3.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc3 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename3);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename4 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
    // MPI_File_delete(filename4.c_str(), info);
    int rc4 = MPI_File_open(MPI_COMM_WORLD, filename4.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc4 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename4);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename5 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
    // MPI_File_delete(filename5.c_str(), info);
    int rc5 = MPI_File_open(MPI_COMM_WORLD, filename5.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc5 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename5);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename6 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
    // MPI_File_delete(filename6.c_str(), info);
    int rc6 = MPI_File_open(MPI_COMM_WORLD, filename6.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc6 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename6);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename7 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
    // MPI_File_delete(filename7.c_str(), info);
    int rc7 = MPI_File_open(MPI_COMM_WORLD, filename7.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc7 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename7);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename8 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
    // MPI_File_delete(filename8.c_str(), info);
    int rc8 = MPI_File_open(MPI_COMM_WORLD, filename8.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc8 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename8);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename9 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
    // MPI_File_delete(filename9.c_str(), info);
    int rc9 = MPI_File_open(MPI_COMM_WORLD, filename9.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc9 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename9);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpPhaseField1.bin";
    int rc10 = MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc10 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename10);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename11 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpPhaseField2.bin";
    int rc11 = MPI_File_open(MPI_COMM_WORLD, filename11.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc11 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename11);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename12 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpPressureField.bin";
    int rc12 = MPI_File_open(MPI_COMM_WORLD, filename12.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc12 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename12);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

}

void MPIIOSimulationObserver::writeCpTimeStep(int step)
{
    if (comm->isRoot()) {
        UbFileOutputASCII f(path + "/mpi_io_cp/cp.txt");
        f.writeInteger(step);
    }
}
//////////////////////////////////////////////////////////////////////////
int MPIIOSimulationObserver::readCpTimeStep()
{
    UbFileInputASCII f(path + "/mpi_io_cp/cp.txt");
    int step = f.readInteger();
    return step;
}
//////////////////////////////////////////////////////////////////////////

//! \}
