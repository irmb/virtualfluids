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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MPIIOMigrationSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Alena Karanchuk
//=======================================================================================
#include "MPIIOMigrationSimulationObserver.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include <parallel/Communicator.h>
#include "CoordinateTransformation3D.h"
#include "EsoSplit.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "Grid3DVisitor.h"
#include "PointerDefinitions.h"
#include "RenumberBlockVisitor.h"
#include "UbFileInputASCII.h"
#include "UbFileOutputASCII.h"
#include "UbScheduler.h"
#include "WbWriter.h"
#include <MemoryUtil.h>
#include <UbSystem.h>

using namespace MPIIODataStructures;

MPIIOMigrationSimulationObserver::MPIIOMigrationSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Grid3DVisitor> mV, const std::string &path, std::shared_ptr<vf::parallel::Communicator> comm)
    : MPIIOSimulationObserver(grid, s, path, comm)
{
    memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));
    metisVisitor = mV;

    //-------------------------   define MPI types  ---------------------------------

    MPI_Datatype typesDataSet[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
    int blocksDataSet[3]         = { 5, 2, 2 };
    MPI_Aint offsetsDatatSet[3], lbDataSet, extentDataSet;

    offsetsDatatSet[0] = 0;
    MPI_Type_get_extent(MPI_DOUBLE, &lbDataSet, &extentDataSet);
    offsetsDatatSet[1] = blocksDataSet[0] * extentDataSet;

    MPI_Type_get_extent(MPI_INT, &lbDataSet, &extentDataSet);
    offsetsDatatSet[2] = offsetsDatatSet[1] + blocksDataSet[1] * extentDataSet;

    MPI_Type_create_struct(3, blocksDataSet, offsetsDatatSet, typesDataSet, &dataSetType);
    MPI_Type_commit(&dataSetType);

    //-----------------------------------------------------------------------

    MPI_Type_contiguous(1, MPI_INT, &dataSetSmallType);
    MPI_Type_commit(&dataSetSmallType);

    //-----------------------------------------------------------------------

    MPI_Type_contiguous(4, MPI_INT, &boundCondParamType);
    MPI_Type_commit(&boundCondParamType);

    //---------------------------------------

    MPI_Type_contiguous(3, MPI_INT, &boundCondTypeAdd);
    MPI_Type_commit(&boundCondTypeAdd);
}
//////////////////////////////////////////////////////////////////////////
MPIIOMigrationSimulationObserver::~MPIIOMigrationSimulationObserver()
{
    MPI_Type_free(&dataSetType);
    MPI_Type_free(&dataSetSmallType);
    MPI_Type_free(&boundCondParamType);
    MPI_Type_free(&boundCondTypeAdd);
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationSimulationObserver::update(real step)
{
    if (scheduler->isDue(step)) 
    {
        if (comm->isRoot())
            UBLOG(logINFO, "MPIIOMigrationSimulationObserver save step: " << step);
        if (comm->isRoot())
            UBLOG(logINFO, "Save check point - start");
        /*if (comm->isRoot())*/ clearAllFiles((int)step);

        writeBlocks((int)step);
        writeDataSet((int)step);
        writeBoundaryConds((int)step);

        writeCpTimeStep((int)step);

        if (comm->isRoot())
            UBLOG(logINFO, "Save check point - end");
    }
}

void MPIIOMigrationSimulationObserver::clearAllFiles(int step)
{
    MPI_File file_handler;
    MPI_Info info       = MPI_INFO_NULL;
    MPI_Offset new_size = 0;

    MPIIOSimulationObserver::clearAllFiles(step);

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

    std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc10 = MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc10 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename10);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);
}

void MPIIOMigrationSimulationObserver::writeBlocks(int step)
{
    grid->renumberBlockIDs();
    MPIIOSimulationObserver::writeBlocks(step);
}

void MPIIOMigrationSimulationObserver::writeDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks, that belong to this process

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;
    DataSetMigration *dataSetArray = new DataSetMigration[blocksCount];
    std::vector<real> doubleValuesArrayF; // real-values (arrays of f's) in all blocks  Fdistribution
    std::vector<real> doubleValuesArrayH1; // real-values (arrays of f's) in all blocks  H1distribution
    std::vector<real> doubleValuesArrayH2; // real-values (arrays of f's) in all blocks  H2distribution

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeDataSet start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool multiPhase1 = false;
    bool multiPhase2 = false;
    DSArraysPresence arrPresence;
    memset(&arrPresence, 0, sizeof(arrPresence));
    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<EsoSplit> D3Q27EsoTwist3DSplittedVectorPtrF = 0, D3Q27EsoTwist3DSplittedVectorPtrH1 = 0, D3Q27EsoTwist3DSplittedVectorPtrH2 = 0;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF = 0, localDistributionsH1 = 0, localDistributionsH2 = 0;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF = 0, nonLocalDistributionsH1 = 0, nonLocalDistributionsH2 = 0;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsF = 0, zeroDistributionsH1 = 0, zeroDistributionsH2 = 0;

    SPtr<LBMKernel> kernel;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            kernel = dynamicPointerCast<LBMKernel>(block->getKernel());

            dataSetArray[ic].globalID = block->getGlobalID(); // id of the block needed to find it while regenerating the grid
            dataSetArray[ic].ghostLayerWidth = kernel->getGhostLayerWidth();
            dataSetArray[ic].collFactor = kernel->getCollisionFactor();
            dataSetArray[ic].deltaT = kernel->getDeltaT();
            dataSetArray[ic].compressible = kernel->getCompressible();
            dataSetArray[ic].withForcing = kernel->getWithForcing();
            dataSetArray[ic].collFactorL = kernel->getCollisionFactorL();
            dataSetArray[ic].collFactorG = kernel->getCollisionFactorG();
            dataSetArray[ic].densityRatio = kernel->getDensityRatio();

            D3Q27EsoTwist3DSplittedVectorPtrF = dynamicPointerCast<EsoSplit>(block->getKernel()->getDataSet()->getFdistributions());
            localDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getLocalDistributions();
            nonLocalDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getNonLocalDistributions();
            zeroDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getZeroDistributions();

            D3Q27EsoTwist3DSplittedVectorPtrH1 = dynamicPointerCast<EsoSplit>(block->getKernel()->getDataSet()->getHdistributions());
            if (D3Q27EsoTwist3DSplittedVectorPtrH1 != 0)
            {
                multiPhase1 = true;
                localDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getLocalDistributions();
                nonLocalDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getNonLocalDistributions();
                zeroDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getZeroDistributions();
            }

            D3Q27EsoTwist3DSplittedVectorPtrH2 = dynamicPointerCast<EsoSplit>(block->getKernel()->getDataSet()->getH2distributions());
            if (D3Q27EsoTwist3DSplittedVectorPtrH2 != 0)
            {
                multiPhase2 = true;
                localDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getLocalDistributions();
                nonLocalDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getNonLocalDistributions();
                zeroDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getZeroDistributions();
            }

            if (firstBlock) // && block->getKernel()) // when first (any) valid block...
            {
                if (localDistributionsF)
                {
                    dataSetParamStr1.nx[0] = static_cast<int>(localDistributionsF->getNX1());
                    dataSetParamStr1.nx[1] = static_cast<int>(localDistributionsF->getNX2());
                    dataSetParamStr1.nx[2] = static_cast<int>(localDistributionsF->getNX3());
                    dataSetParamStr1.nx[3] = static_cast<int>(localDistributionsF->getNX4());
                }

                if (nonLocalDistributionsF)
                {
                    dataSetParamStr2.nx[0] = static_cast<int>(nonLocalDistributionsF->getNX1());
                    dataSetParamStr2.nx[1] = static_cast<int>(nonLocalDistributionsF->getNX2());
                    dataSetParamStr2.nx[2] = static_cast<int>(nonLocalDistributionsF->getNX3());
                    dataSetParamStr2.nx[3] = static_cast<int>(nonLocalDistributionsF->getNX4());
                }
                if (zeroDistributionsF)
                {
                    dataSetParamStr3.nx[0] = static_cast<int>(zeroDistributionsF->getNX1());
                    dataSetParamStr3.nx[1] = static_cast<int>(zeroDistributionsF->getNX2());
                    dataSetParamStr3.nx[2] = static_cast<int>(zeroDistributionsF->getNX3());
                    dataSetParamStr3.nx[3] = 1;
                }

                // ... than save some parameters that are equal in all blocks
                dataSetParamStr1.nx1 = dataSetParamStr2.nx1 = dataSetParamStr3.nx1 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX1());
                dataSetParamStr1.nx2 = dataSetParamStr2.nx2 = dataSetParamStr3.nx2 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX2());
                dataSetParamStr1.nx3 = dataSetParamStr2.nx3 = dataSetParamStr3.nx3 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX3());

                doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
                    dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
                    dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

                SPtr<CbArray4D<real, IndexerX4X3X2X1>> averageDensityArray = block->getKernel()->getDataSet()->getAverageDensity();
                if (averageDensityArray)
                    arrPresence.isAverageDensityArrayPresent = true;
                else
                    arrPresence.isAverageDensityArrayPresent = false;

                SPtr<CbArray4D<real, IndexerX4X3X2X1>> AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
                if (AverageVelocityArray3DPtr)
                    arrPresence.isAverageVelocityArrayPresent = true;
                else
                    arrPresence.isAverageVelocityArrayPresent = false;

                SPtr<CbArray4D<real, IndexerX4X3X2X1>> AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
                if (AverageFluctArray3DPtr)
                    arrPresence.isAverageFluktuationsArrayPresent = true;
                else
                    arrPresence.isAverageFluktuationsArrayPresent = false;

                SPtr<CbArray4D<real, IndexerX4X3X2X1>> AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
                if (AverageTripleArray3DPtr)
                    arrPresence.isAverageTripleArrayPresent = true;
                else
                    arrPresence.isAverageTripleArrayPresent = false;

                SPtr<CbArray4D<real, IndexerX4X3X2X1>> ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
                if (ShearStressValArray3DPtr)
                    arrPresence.isShearStressValArrayPresent = true;
                else
                    arrPresence.isShearStressValArrayPresent = false;

                SPtr<CbArray3D<real, IndexerX3X2X1>> relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
                if (relaxationFactor3DPtr)
                    arrPresence.isRelaxationFactorPresent = true;
                else
                    arrPresence.isRelaxationFactorPresent = false;

                SPtr<CbArray3D<real, IndexerX3X2X1>> phaseField3DPtr1 = block->getKernel()->getDataSet()->getPhaseField();
                if (phaseField3DPtr1)
                    arrPresence.isPhaseField1Present = true;
                else
                    arrPresence.isPhaseField1Present = false;

                SPtr<CbArray3D<real, IndexerX3X2X1>> phaseField3DPtr2 = block->getKernel()->getDataSet()->getPhaseField2();
                if (phaseField3DPtr2)
                    arrPresence.isPhaseField2Present = true;
                else
                    arrPresence.isPhaseField2Present = false;

                SPtr<CbArray3D<real, IndexerX3X2X1>> pressureFieldPtr = block->getKernel()->getDataSet()->getPressureField();
                if (pressureFieldPtr)
                    arrPresence.isPressureFieldPresent = true;
                else
                    arrPresence.isPressureFieldPresent = false;

                firstBlock = false;
            }

            if (localDistributionsF && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), localDistributionsF->getDataVector().begin(), localDistributionsF->getDataVector().end());
            if (nonLocalDistributionsF && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), nonLocalDistributionsF->getDataVector().begin(), nonLocalDistributionsF->getDataVector().end());
            if (zeroDistributionsF && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), zeroDistributionsF->getDataVector().begin(), zeroDistributionsF->getDataVector().end());

            if (multiPhase1)
            {
                if (localDistributionsH1 && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), localDistributionsH1->getDataVector().begin(), localDistributionsH1->getDataVector().end());
                if (nonLocalDistributionsH1 && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), nonLocalDistributionsH1->getDataVector().begin(), nonLocalDistributionsH1->getDataVector().end());
                if (zeroDistributionsH1 && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), zeroDistributionsH1->getDataVector().begin(), zeroDistributionsH1->getDataVector().end());
            }

            if (multiPhase2)
            {
                if (localDistributionsH2 && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                    doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), localDistributionsH2->getDataVector().begin(), localDistributionsH2->getDataVector().end());
                if (nonLocalDistributionsH2 && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                    doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), nonLocalDistributionsH2->getDataVector().begin(), nonLocalDistributionsH2->getDataVector().end());
                if (zeroDistributionsH2 && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                    doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), zeroDistributionsH2->getDataVector().begin(), zeroDistributionsH2->getDataVector().end());
            }
            
            ic++;
        }
    }

    // register new MPI-type depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

    // write to the file
    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetF.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_write_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
    
    MPI_Offset write_offset;
    size_t sizeofOneDataSet = sizeof(DataSetMigration) + doubleCountInBlock * sizeof(real);

    for (int nb = 0; nb < blocksCount; nb++) 
    {
        write_offset = (MPI_Offset)(3 * sizeof(dataSetParam) + dataSetArray[nb].globalID * sizeofOneDataSet);
        MPI_File_write_at(file_handler, write_offset, &dataSetArray[nb], 1, dataSetType, MPI_STATUS_IGNORE);
        MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetMigration)), &doubleValuesArrayF[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
    }

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);

    //-------------------------------- H1 ----------------------------------------------------
    if (multiPhase1)
    {
        filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
        rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
        if (rc != MPI_SUCCESS)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);

        sizeofOneDataSet = doubleCountInBlock * sizeof(real);

        for (int nb = 0; nb < blocksCount; nb++) 
        {
            write_offset = (MPI_Offset)(dataSetArray[nb].globalID * sizeofOneDataSet);
            MPI_File_write_at(file_handler, write_offset, &doubleValuesArrayH1[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
        }

        MPI_File_sync(file_handler);
        MPI_File_close(&file_handler);
    }

    //-------------------------------- H2 ----------------------------------------------------
    if (multiPhase2)
    {
        filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH2.bin";
        rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
        if (rc != MPI_SUCCESS)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);

        sizeofOneDataSet = doubleCountInBlock * sizeof(real);

        for (int nb = 0; nb < blocksCount; nb++) 
        {
            write_offset = (MPI_Offset)(dataSetArray[nb].globalID * sizeofOneDataSet);
            MPI_File_write_at(file_handler, write_offset, &doubleValuesArrayH2[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
        }

        MPI_File_sync(file_handler);
        MPI_File_close(&file_handler);
    }
    //--------------------------------

    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeDataSet time: " << finish - start << " s");
    }

    delete[] dataSetArray;

    MPI_File file_handler1;
    std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler1);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename1);
    MPI_File_write_at(file_handler1, (MPI_Offset)0, &arrPresence, 1, arrayPresenceType, MPI_STATUS_IGNORE);
    MPI_File_sync(file_handler1);
    MPI_File_close(&file_handler1);

    if (arrPresence.isAverageDensityArrayPresent)
        write4DArray(step, AverageDensity, std::string("/cpAverageDensityArray.bin"));

    if (arrPresence.isAverageVelocityArrayPresent)
        write4DArray(step, AverageVelocity, std::string("/cpAverageVelocityArray.bin"));

    if (arrPresence.isAverageFluktuationsArrayPresent)
        write4DArray(step, AverageFluktuations, std::string("/cpAverageFluktuationsArray.bin"));

    if (arrPresence.isAverageTripleArrayPresent)
        write4DArray(step, AverageTriple, std::string("/cpAverageTripleArray.bin"));

    if (arrPresence.isShearStressValArrayPresent)
        write4DArray(step, ShearStressVal, std::string("/cpShearStressValArray.bin"));

    if (arrPresence.isRelaxationFactorPresent)
        write3DArray(step, RelaxationFactor, std::string("/cpRelaxationFactor.bin"));

    if (arrPresence.isPhaseField1Present)
        write3DArray(step, PhaseField1, std::string("/cpPhaseField1.bin"));

    if (arrPresence.isPhaseField2Present)
        write3DArray(step, PhaseField2, std::string("/cpPhaseField2.bin"));

    if (arrPresence.isPressureFieldPresent)
        write3DArray(step, PressureField, std::string("/cpPressureField.bin"));

}

void MPIIOMigrationSimulationObserver::write4DArray(int step, Arrays arrayType, std::string fname)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallMigration *dataSetSmallArray = new DataSetSmallMigration[blocksCount];
    std::vector<real> doubleValuesArray; // real-values of the AverageDensityArray in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write4DArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<CbArray4D<real, IndexerX4X3X2X1>> ___Array;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            dataSetSmallArray[ic].globalID = block->getGlobalID(); // id of the block needed to find it while regenerating the grid

            switch (arrayType) 
            {
                case AverageDensity:
                    ___Array = block->getKernel()->getDataSet()->getAverageDensity();
                    break;
                case AverageVelocity:
                    ___Array = block->getKernel()->getDataSet()->getAverageVelocity();
                    break;
                case AverageFluktuations:
                    ___Array = block->getKernel()->getDataSet()->getAverageFluctuations();
                    break;
                case AverageTriple:
                    ___Array = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
                    break;
                case ShearStressVal:
                    ___Array = block->getKernel()->getDataSet()->getShearStressValues();
                    break;
                default:
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationSimulationObserver::write4DArray : 4D array type does not exist!"));
                    break;
            }

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(___Array->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(___Array->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(___Array->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(___Array->getNX4());
                doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if (___Array && (dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) && (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), ___Array->getDataVector().begin(), ___Array->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write4DArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

    MPI_Offset write_offset;
    size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(real);

    for (int nb = 0; nb < blocksCount; nb++) 
    {
        write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
        MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
        MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)),
                          &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
    }

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write4DArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIOMigrationSimulationObserver::write3DArray(int step, Arrays arrayType, std::string fname)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallMigration *dataSetSmallArray = new DataSetSmallMigration[blocksCount];
    std::vector<real> doubleValuesArray; // real-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write3DArray start collect data to file = " << fname);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<CbArray3D<real, IndexerX3X2X1>> ___Array;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            dataSetSmallArray[ic].globalID = block->getGlobalID(); // id of the block needed to find it while regenerating the grid

            switch (arrayType) 
            {
                case RelaxationFactor:
                    ___Array = block->getKernel()->getDataSet()->getRelaxationFactor();
                    break;
                case PhaseField1:
                    ___Array = block->getKernel()->getDataSet()->getPhaseField();
                    break;
                case PhaseField2:
                    ___Array = block->getKernel()->getDataSet()->getPhaseField2();
                    break;
                case PressureField:
                    ___Array = block->getKernel()->getDataSet()->getPressureField();
                    break;
                default:
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationSimulationObserver::write3DArray : 3D array type does not exist!"));
                    break;
            }

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(___Array->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(___Array->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(___Array->getNX3());
                dataSetParamStr.nx[3] = 1;
                doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if (___Array && (dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), ___Array->getDataVector().begin(), ___Array->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write3DArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

    size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(real);

    MPI_Offset write_offset;
    for (int nb = 0; nb < blocksCount; nb++) 
    {
        write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
        MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
        MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)),
                          &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
    }

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::write3DArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIOMigrationSimulationObserver::writeBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeBoundaryConds start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    int blocksCount          = 0; // quantity of blocks, that belong to this process
    size_t allBytesCount     = 0; // quantity of bytes, that one process writes to the file

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    BCAddMigration *bcAddArray = new BCAddMigration[blocksCount];
    size_t *bytesCount         = new size_t[blocksCount]; // quantity of bytes, that each block writes to the file
    std::vector<BoundaryCondition> *bcVector = new std::vector<BoundaryCondition>[blocksCount];
    std::vector<int> *bcindexmatrixVector    = new std::vector<int>[blocksCount];
    std::vector<int> *indexContainerVector   = new std::vector<int>[blocksCount];

    bool bcindexmatrixCountNotInit = true;
    int ic                         = 0;
    SPtr<BCArray3D> bcArr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) // all the blocks of the current level
        {
            bcArr = block->getKernel()->getBCSet()->getBCArray();

            bcAddArray[ic].globalID = block->getGlobalID();                // id of the block needed to find it while regenerating the grid
            bcAddArray[ic].boundCond_count      = 0; // how many BoundaryConditions in this block
            bcAddArray[ic].indexContainer_count = 0; // how many indexContainer-values in this block
            bytesCount[ic] = sizeof(BCAddMigration);
            bcVector[ic].resize(0);
            bcindexmatrixVector[ic].resize(0);
            indexContainerVector[ic].resize(0);

            for (std::size_t bc = 0; bc < bcArr->getBCVectorSize(); bc++) 
            {
                BoundaryCondition *bouCond = new BoundaryCondition();
                if (bcArr->bcvector[bc] == NULL) 
                    memset(bouCond, 0, sizeof(BoundaryCondition));
                else 
                {
                    bouCond->noslipBoundaryFlags    = bcArr->bcvector[bc]->getNoSlipBoundary();
                    bouCond->slipBoundaryFlags      = bcArr->bcvector[bc]->getSlipBoundary();
                    bouCond->velocityBoundaryFlags  = bcArr->bcvector[bc]->getVelocityBoundary();
                    bouCond->densityBoundaryFlags   = bcArr->bcvector[bc]->getDensityBoundary();
                    bouCond->wallModelBoundaryFlags = bcArr->bcvector[bc]->getWallModelBoundary();
                    bouCond->bcVelocityX1           = bcArr->bcvector[bc]->getBoundaryVelocityX1();
                    bouCond->bcVelocityX2           = bcArr->bcvector[bc]->getBoundaryVelocityX2();
                    bouCond->bcVelocityX3           = bcArr->bcvector[bc]->getBoundaryVelocityX3();
                    bouCond->bcDensity              = bcArr->bcvector[bc]->getBoundaryDensity();
                    bouCond->bcPhaseField           = bcArr->bcvector[bc]->getBoundaryPhaseField();
                    bouCond->nx1                    = bcArr->bcvector[bc]->nx1;
                    bouCond->nx2                    = bcArr->bcvector[bc]->nx2;
                    bouCond->nx3                    = bcArr->bcvector[bc]->nx3;
                    for (int iq = 0; iq < 26; iq++)
                        bouCond->q[iq] = bcArr->bcvector[bc]->getQ(iq);
                    bouCond->bcStrategyKey = bcArr->bcvector[bc]->getBCStrategyKey();
                }

                bcVector[ic].push_back(*bouCond);
                bcAddArray[ic].boundCond_count++;
                bytesCount[ic] += sizeof(BoundaryCondition);
             }

            if (bcindexmatrixCountNotInit) 
            {
                boundCondParamStr.nx1                = static_cast<int>(bcArr->bcindexmatrix.getNX1());
                boundCondParamStr.nx2                = static_cast<int>(bcArr->bcindexmatrix.getNX2());
                boundCondParamStr.nx3                = static_cast<int>(bcArr->bcindexmatrix.getNX3());
                boundCondParamStr.bcindexmatrixCount = static_cast<int>(bcArr->bcindexmatrix.getDataVector().size());
                bcindexmatrixCountNotInit            = false;
            }
            bcindexmatrixVector[ic].insert(bcindexmatrixVector[ic].begin(), bcArr->bcindexmatrix.getDataVector().begin(), bcArr->bcindexmatrix.getDataVector().end());
            bytesCount[ic] += boundCondParamStr.bcindexmatrixCount * sizeof(int);

            indexContainerVector[ic].insert(indexContainerVector[ic].begin(), bcArr->indexContainer.begin(), bcArr->indexContainer.end());
            bcAddArray[ic].indexContainer_count = static_cast<int>(bcArr->indexContainer.size());
            bytesCount[ic] += bcAddArray[ic].indexContainer_count * sizeof(int);

            allBytesCount += bytesCount[ic];

            ic++;
        }
    }

    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "<< Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;
    // MPI_Info_create (&info);
    // MPI_Info_set(info,"romio_cb_write","enable");
    // MPI_Info_set(info,"cb_buffer_size","4194304");
    // MPI_Info_set(info,"striping_unit","4194304");

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset = (MPI_Offset)(sizeof(boundCondParam) + grid->getNumberOfBlocks() * sizeof(size_t));
    size_t next_file_offset = 0;
    if (size > 1) 
    {
        if (rank == 0) 
        {
            next_file_offset = write_offset + allBytesCount;
            MPI_Send(&next_file_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } 
        else 
        {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_file_offset = write_offset + allBytesCount;
            if (rank < size - 1)
                MPI_Send(&next_file_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_write_at(file_handler, 0, &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);

    MPI_Offset write_offsetIndex;

    for (int nb = 0; nb < blocksCount; nb++) 
    {
        write_offsetIndex = (MPI_Offset)(sizeof(boundCondParam) + bcAddArray[nb].globalID * sizeof(size_t));
        MPI_File_write_at(file_handler, write_offsetIndex, &write_offset, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);

        MPI_File_write_at(file_handler, write_offset, &bcAddArray[nb], 1, boundCondTypeAdd, MPI_STATUS_IGNORE);
        if (bcVector[nb].size() > 0)
            MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAddMigration)), &bcVector[nb][0],
                              bcAddArray[nb].boundCond_count, boundCondType, MPI_STATUS_IGNORE);

        if (bcindexmatrixVector[nb].size() > 0)
            MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAddMigration) + bcAddArray[nb].boundCond_count * sizeof(BoundaryCondition)),
                              &bcindexmatrixVector[nb][0], 1, bcindexmatrixType, MPI_STATUS_IGNORE);

        if (indexContainerVector[nb].size() > 0)
            MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAddMigration) + bcAddArray[nb].boundCond_count * sizeof(BoundaryCondition) +
                              boundCondParamStr.bcindexmatrixCount * sizeof(int)), &indexContainerVector[nb][0], bcAddArray[nb].indexContainer_count, MPI_INT,
                              MPI_STATUS_IGNORE);

        write_offset += bytesCount[nb];
    }

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::writeBoundaryConds time: " << finish - start << " s");
    }

    delete[] bcAddArray;
    delete[] bytesCount;
    delete[] bcVector;
    delete[] bcindexmatrixVector;
    delete[] indexContainerVector;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIOMigrationSimulationObserver::restart(int step)
{
    if (comm->isRoot())
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver restart step: " << step);
    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - start");

    readBlocks(step);

    grid->accept(metisVisitor);

    readDataSet(step);
    readBoundaryConds(step);

    grid->setTimeStep(step);

    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - end");
}

void MPIIOMigrationSimulationObserver::readBlocks(int step) { MPIIOSimulationObserver::readBlocks(step); }

void MPIIOMigrationSimulationObserver::readDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    bool multiPhase1 = false;
    bool multiPhase2 = false;
    size_t blocksCount = 0; // quantity of the blocks, that belong to this process
    dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;

    // read from the grid the blocks, that belong to this process
    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetMigration *dataSetArray = new DataSetMigration[blocksCount];

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetF.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);

    size_t doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
        dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
        dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];
    std::vector<real> doubleValuesArrayF(size_t(blocksCount * doubleCountInBlock)); // real-values in all blocks  Fdistributions
    std::vector<real> doubleValuesArrayH1; // real-values in all blocks  H1distributions
    std::vector<real> doubleValuesArrayH2; // real-values in all blocks  H2distributions

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    size_t ic = 0;
    MPI_Offset read_offset;
    size_t sizeofOneDataSet = size_t(sizeof(DataSetMigration) + doubleCountInBlock * sizeof(real));

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            read_offset = (MPI_Offset)(3 * sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
            MPI_File_read_at(file_handler, read_offset, &dataSetArray[ic], 1, dataSetType, MPI_STATUS_IGNORE);
            MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetMigration)),
                             &doubleValuesArrayF[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
            ic++;
        }
    }

    MPI_File_close(&file_handler);

    //----------------------------------------- H1 ----------------------------------------------------
    ic = 0;
    MPI_Offset fsize;
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);
    MPI_File_get_size(file_handler, &fsize);
    if (fsize > 0)
    {
        multiPhase1 = true;
        doubleValuesArrayH1.resize(blocksCount * doubleCountInBlock);

        sizeofOneDataSet = size_t(doubleCountInBlock * sizeof(real));

        for (int level = minInitLevel; level <= maxInitLevel; level++)
        {
            for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
            {
                read_offset = (MPI_Offset)(block->getGlobalID() * sizeofOneDataSet);
                MPI_File_read_at(file_handler, read_offset, &doubleValuesArrayH1[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
                ic++;
            }
        }

    }
    MPI_File_close(&file_handler);
    //----------------------------------------- H2 ----------------------------------------------------
    ic = 0;
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH2.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_get_size(file_handler, &fsize);
    if (fsize > 0)
    {
        multiPhase2 = true;
        doubleValuesArrayH2.resize(blocksCount * doubleCountInBlock);

        sizeofOneDataSet = size_t(doubleCountInBlock * sizeof(real));

        for (int level = minInitLevel; level <= maxInitLevel; level++)
        {
            for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
            {
                read_offset = (MPI_Offset)(block->getGlobalID() * sizeofOneDataSet);
                MPI_File_read_at(file_handler, read_offset, &doubleValuesArrayH2[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
                ic++;
            }
        }

    }
    MPI_File_close(&file_handler);

    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readDataSet time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readDataSet start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    std::vector<real> vectorsOfValuesF1, vectorsOfValuesF2, vectorsOfValuesF3;
    std::vector<real> vectorsOfValuesH11, vectorsOfValuesH12, vectorsOfValuesH13;
    std::vector<real> vectorsOfValuesH21, vectorsOfValuesH22, vectorsOfValuesH23;

    size_t vectorSize1 = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3];
    size_t vectorSize2 = dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3];
    size_t vectorSize3 = dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

    for (std::size_t n = 0; n < blocksCount; n++) 
    {
        vectorsOfValuesF1.assign(doubleValuesArrayF.data() + index, doubleValuesArrayF.data() + index + vectorSize1);
        if(multiPhase1)
            vectorsOfValuesH11.assign(doubleValuesArrayH1.data() + index, doubleValuesArrayH1.data() + index + vectorSize1);
        if (multiPhase2)
            vectorsOfValuesH21.assign(doubleValuesArrayH2.data() + index, doubleValuesArrayH2.data() + index + vectorSize1);
        index += vectorSize1;

        vectorsOfValuesF2.assign(doubleValuesArrayF.data() + index, doubleValuesArrayF.data() + index + vectorSize2);
        if (multiPhase1)
            vectorsOfValuesH12.assign(doubleValuesArrayH1.data() + index, doubleValuesArrayH1.data() + index + vectorSize2);
        if (multiPhase2)
            vectorsOfValuesH22.assign(doubleValuesArrayH2.data() + index, doubleValuesArrayH2.data() + index + vectorSize2);
        index += vectorSize2;

        vectorsOfValuesF3.assign(doubleValuesArrayF.data() + index, doubleValuesArrayF.data() + index + vectorSize3);
        if (multiPhase1)
            vectorsOfValuesH13.assign(doubleValuesArrayH1.data() + index, doubleValuesArrayH1.data() + index + vectorSize3);
        if (multiPhase2)
            vectorsOfValuesH23.assign(doubleValuesArrayH2.data() + index, doubleValuesArrayH2.data() + index + vectorSize3);
        index += vectorSize3;
 
        SPtr<DistributionArray3D> mFdistributions(new EsoSplit());
        dynamicPointerCast<EsoSplit>(mFdistributions)->setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesF1, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
        dynamicPointerCast<EsoSplit>(mFdistributions)->setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesF2, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
        dynamicPointerCast<EsoSplit>(mFdistributions)->setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
            vectorsOfValuesF3, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));
        
        //----------------------------------------- H1 ----------------------------------------------------
       SPtr<DistributionArray3D> mH1distributions(new EsoSplit());
       if (multiPhase1)
        {
            dynamicPointerCast<EsoSplit>(mH1distributions)->setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesH11, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
            dynamicPointerCast<EsoSplit>(mH1distributions)->setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesH12, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
            dynamicPointerCast<EsoSplit>(mH1distributions)->setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                vectorsOfValuesH13, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

            dynamicPointerCast<EsoSplit>(mH1distributions)->setNX1(dataSetParamStr1.nx1);
            dynamicPointerCast<EsoSplit>(mH1distributions)->setNX2(dataSetParamStr1.nx2);
            dynamicPointerCast<EsoSplit>(mH1distributions)->setNX3(dataSetParamStr1.nx3);
         }

        SPtr<DistributionArray3D> mH2distributions(new EsoSplit());
        if (multiPhase2)
        {
            dynamicPointerCast<EsoSplit>(mH2distributions)->setLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesH21, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
            dynamicPointerCast<EsoSplit>(mH2distributions)->setNonLocalDistributions(CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<real, IndexerX4X3X2X1>(vectorsOfValuesH22, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
            dynamicPointerCast<EsoSplit>(mH2distributions)->setZeroDistributions(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                    vectorsOfValuesH23, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

            dynamicPointerCast<EsoSplit>(mH2distributions)->setNX1(dataSetParamStr1.nx1);
            dynamicPointerCast<EsoSplit>(mH2distributions)->setNX2(dataSetParamStr1.nx2);
            dynamicPointerCast<EsoSplit>(mH2distributions)->setNX3(dataSetParamStr1.nx3);
        }

        dynamicPointerCast<EsoSplit>(mFdistributions)->setNX1(dataSetParamStr1.nx1);
        dynamicPointerCast<EsoSplit>(mFdistributions)->setNX2(dataSetParamStr1.nx2);
        dynamicPointerCast<EsoSplit>(mFdistributions)->setNX3(dataSetParamStr1.nx3);

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetArray[n].globalID);
        this->lbmKernel->setBlock(block);
        this->lbmKernel->setNX(std::array<int, 3>{ {dataSetParamStr1.nx1, dataSetParamStr1.nx2, dataSetParamStr1.nx3}});
        UbTupleInt3 blockNX = grid->getBlockNX();
        this->lbmKernel->setNX(std::array<int, 3>{ { val<1>(blockNX), val<2>(blockNX), val<3>(blockNX) } });
        SPtr<LBMKernel> kernel = this->lbmKernel->clone();
        kernel->setGhostLayerWidth(dataSetArray[n].ghostLayerWidth);
        kernel->setCollisionFactor(dataSetArray[n].collFactor);
        kernel->setDeltaT(dataSetArray[n].deltaT);
        kernel->setCompressible(dataSetArray[n].compressible);
        kernel->setWithForcing(dataSetArray[n].withForcing);
        kernel->setCollisionFactorMultiphase(dataSetArray[n].collFactorL, dataSetArray[n].collFactorG);
        kernel->setDensityRatio(dataSetArray[n].densityRatio);

         SPtr<DataSet3D> dataSetPtr = SPtr<DataSet3D>(new DataSet3D());
        dataSetPtr->setFdistributions(mFdistributions);
        if (multiPhase1)
            dataSetPtr->setHdistributions(mH1distributions);
        if (multiPhase2)
            dataSetPtr->setH2distributions(mH2distributions);
        kernel->setDataSet(dataSetPtr);
        block->setKernel(kernel);
    }

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readDataSet end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetArray;

    //-------------------------------------------------------------

    DSArraysPresence arrPresence;
    MPI_File file_handler1;
    std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler1);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename1);
    MPI_File_read_at(file_handler1, (MPI_Offset)0, &arrPresence, 1, arrayPresenceType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler1);

    if (arrPresence.isAverageDensityArrayPresent)
        readArray(step, AverageDensity, std::string("/cpAverageDensityArray.bin"));

    if (arrPresence.isAverageVelocityArrayPresent)
        readArray(step, AverageVelocity, std::string("/cpAverageVelocityArray.bin"));

    if (arrPresence.isAverageFluktuationsArrayPresent)
        readArray(step, AverageFluktuations, std::string("/cpAverageFluktuationsArray.bin"));

    if (arrPresence.isAverageTripleArrayPresent)
        readArray(step, AverageTriple, std::string("/cpAverageTripleArray.bin"));

    if (arrPresence.isShearStressValArrayPresent)
        readArray(step, ShearStressVal, std::string("/cpShearStressValArray.bin"));

    if (arrPresence.isRelaxationFactorPresent)
        readArray(step, RelaxationFactor, std::string("/cpRelaxationFactor.bin"));
 
    if (arrPresence.isPhaseField1Present)
        readArray(step, PhaseField1, std::string("/cpPhaseField1.bin"));

    if (arrPresence.isPhaseField2Present)
        readArray(step, PhaseField2, std::string("/cpPhaseField2.bin"));

    if (arrPresence.isPressureFieldPresent)
        readArray(step, PressureField, std::string("/cpPressureField.bin"));

}

void MPIIOMigrationSimulationObserver::readArray(int step, Arrays arrType, std::string fname)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readArray start fname = " << fname);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    size_t blocksCount = 0;
    dataSetParam dataSetParamStr;
    memset(&dataSetParamStr, 0, sizeof(dataSetParam));

    // read from the grid the blocks, that belong to this process
    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

    DataSetSmallMigration *dataSetSmallArray = new DataSetSmallMigration[blocksCount];
    size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<real> doubleValuesArray(blocksCount * doubleCountInBlock); // real-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    size_t ic = 0;
    MPI_Offset read_offset;
    size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(real));

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
            MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
            MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)),
                             &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
            ic++;
        }
    }

    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readArray readArray: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    //----------------------------- restore data ---------------------------------
    size_t index = 0;
    size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<real> vectorsOfValues;
    SPtr<CbArray4D<real, IndexerX4X3X2X1>> ___4DArray;
    SPtr<CbArray3D<real, IndexerX3X2X1>> ___3DArray;

    for (std::size_t n = 0; n < blocksCount; n++)
    {
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);

        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill arrays
        switch (arrType) 
        {
            case AverageDensity:
                ___4DArray = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<real, IndexerX4X3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                block->getKernel()->getDataSet()->setAverageDensity(___4DArray);
                break;
            case AverageVelocity:
                ___4DArray = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<real, IndexerX4X3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                block->getKernel()->getDataSet()->setAverageVelocity(___4DArray);
                break;
            case AverageFluktuations:
                ___4DArray = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<real, IndexerX4X3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                block->getKernel()->getDataSet()->setAverageFluctuations(___4DArray);
                break;
            case AverageTriple:
                ___4DArray = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<real, IndexerX4X3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                block->getKernel()->getDataSet()->setAverageTriplecorrelations(___4DArray);
                break;
            case ShearStressVal:
                ___4DArray = CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<real, IndexerX4X3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                block->getKernel()->getDataSet()->setShearStressValues(___4DArray);
                break;
            case RelaxationFactor:
                ___3DArray = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                block->getKernel()->getDataSet()->setRelaxationFactor(___3DArray);
                break;
            case PhaseField1:
                ___3DArray = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                block->getKernel()->getDataSet()->setPhaseField(___3DArray);
                break;
            case PhaseField2:
                ___3DArray = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                block->getKernel()->getDataSet()->setPhaseField2(___3DArray);
                break;
            case PressureField:
                ___3DArray = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(
                    vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                block->getKernel()->getDataSet()->setPressureField(___3DArray);
                break;
            default:
                UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationSimulationObserver::readArray : array type does not exist!"));
                break;
        }
    }

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIOMigrationSimulationObserver::readBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    
    real start {0.};
    real finish {0.};
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    int blocksCount = 0; // quantity of the blocks, that belong to this process

    // read from the grid the blocks, that belong to this process
    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    BCAddMigration *bcAddArray     = new BCAddMigration[blocksCount];
    BoundaryCondition *nullBouCond = new BoundaryCondition();
    memset(nullBouCond, 0, sizeof(BoundaryCondition));
    BoundaryCondition *bcArray;
    int *intArray1;
    int *intArray2;
    std::vector<SPtr<BoundaryConditions>> bcVector;
    std::vector<int> bcindexmatrixV;
    std::vector<int> indexContainerV;

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readBoundaryConds start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    MPI_File_read_at(file_handler, (MPI_Offset)0, &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);
    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    int ic = 0;
    MPI_Offset read_offset1, read_offset2;
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
        {
            read_offset1 = (MPI_Offset)(sizeof(boundCondParam) + block->getGlobalID() * sizeof(size_t));

            MPI_File_read_at(file_handler, read_offset1, &read_offset2, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
            MPI_File_read_at(file_handler, read_offset2, &bcAddArray[ic], 1, boundCondTypeAdd, MPI_STATUS_IGNORE);

            bcArray   = new BoundaryCondition[bcAddArray[ic].boundCond_count];
            intArray1 = new int[boundCondParamStr.bcindexmatrixCount];
            intArray2 = new int[bcAddArray[ic].indexContainer_count];

            if (bcAddArray[ic].boundCond_count > 0) 
            {
                MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAddMigration)), &bcArray[0],
                                 bcAddArray[ic].boundCond_count, boundCondType, MPI_STATUS_IGNORE);
            }
            MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAddMigration) + bcAddArray[ic].boundCond_count * sizeof(BoundaryCondition)),
                             &intArray1[0], 1, bcindexmatrixType, MPI_STATUS_IGNORE);
            if (bcAddArray[ic].indexContainer_count > 0) 
            {
                MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAddMigration) + bcAddArray[ic].boundCond_count * sizeof(BoundaryCondition) +
                                 boundCondParamStr.bcindexmatrixCount * sizeof(int)), &intArray2[0], bcAddArray[ic].indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);
            }

            bcindexmatrixV.resize(0);
            indexContainerV.resize(0);
            bcVector.resize(0);

            for (int ibc = 0; ibc < bcAddArray[ic].boundCond_count; ibc++) 
            {
                SPtr<BoundaryConditions> bc;
                if (memcmp(&bcArray[ibc], nullBouCond, sizeof(BoundaryCondition)) == 0)
                    bc = SPtr<BoundaryConditions>();
                else 
                {
                    bc                         = SPtr<BoundaryConditions>(new BoundaryConditions);
                    bc->noslipBoundaryFlags    = bcArray[ibc].noslipBoundaryFlags;
                    bc->slipBoundaryFlags      = bcArray[ibc].slipBoundaryFlags;
                    bc->densityBoundaryFlags   = bcArray[ibc].densityBoundaryFlags;
                    bc->velocityBoundaryFlags  = bcArray[ibc].velocityBoundaryFlags;
                    bc->wallModelBoundaryFlags = bcArray[ibc].wallModelBoundaryFlags;
                    bc->bcVelocityX1           = bcArray[ibc].bcVelocityX1;
                    bc->bcVelocityX2           = bcArray[ibc].bcVelocityX2;
                    bc->bcVelocityX3           = bcArray[ibc].bcVelocityX3;
                    bc->bcDensity              = bcArray[ibc].bcDensity;
                    bc->bcPhaseField           = bcArray[ibc].bcPhaseField;

                    bc->nx1 = bcArray[ibc].nx1;
                    bc->nx2 = bcArray[ibc].nx2;
                    bc->nx3 = bcArray[ibc].nx3;
                    for (int iq = 0; iq < 26; iq++)
                        bc->setQ(bcArray[ibc].q[iq], iq);
                    bc->setBCStrategyKey(bcArray[ibc].bcStrategyKey);
                }

                bcVector.push_back(bc);
            }

            for (int b1 = 0; b1 < boundCondParamStr.bcindexmatrixCount; b1++)
                bcindexmatrixV.push_back(intArray1[b1]);

            for (int b2 = 0; b2 < bcAddArray[ic].indexContainer_count; b2++)
                indexContainerV.push_back(intArray2[b2]);

            CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);
            SPtr<Block3D> block1 = grid->getBlock(bcAddArray[ic].globalID);

            SPtr<BCSet> bcProc = bcSet->clone(block1->getKernel());
            SPtr<BCArray3D> bcArr(new BCArray3D());
            bcArr->bcindexmatrix  = bcim;
            bcArr->bcvector       = bcVector;
            bcArr->indexContainer = indexContainerV;
            bcProc->setBCArray(bcArr);

            block1->getKernel()->setBCSet(bcProc);

            delete[] bcArray;
            delete[] intArray1;

            ic++;
        }
    }

    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    delete nullBouCond;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationSimulationObserver::readBoundaryConds end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationSimulationObserver::setLBMKernel(SPtr<LBMKernel> kernel) { this->lbmKernel = kernel; }
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationSimulationObserver::setBCSet(SPtr<BCSet> bcSet) { this->bcSet = bcSet; }
