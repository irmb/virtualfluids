#include "MPIIORestartCoProcessor.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Communicator.h"
#include "CoordinateTransformation3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include "LBMKernel.h"
#include "UbFileInputASCII.h"
#include "UbFileOutputASCII.h"
#include "UbScheduler.h"
#include "WbWriter.h"
#include <MemoryUtil.h>
#include <UbSystem.h>

//! BLOCK_SIZE defines the quantity of the BoundaryCondition-structures written as one block to the file
//! To avoid overflow in the parameter \a count of the function MPI_File_write_at
//! structures BoundaryCondition are being written in blocks containing each of them BLOCK_SIZE structures
#define BLOCK_SIZE 1024

using namespace MPIIODataStructures;

MPIIORestartCoProcessor::MPIIORestartCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                 SPtr<Communicator> comm)
    : MPIIOCoProcessor(grid, s, path, comm)
{
    memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));

    //-------------------------   define MPI types  ---------------------------------

    MPI_Datatype typesDataSet[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
    int blocksDataSet[3]         = { 2, 5, 2 };
    MPI_Aint offsetsDatatSet[3], lbDataSet, extentDataSet;

    offsetsDatatSet[0] = 0;
    MPI_Type_get_extent(MPI_DOUBLE, &lbDataSet, &extentDataSet);
    offsetsDatatSet[1] = blocksDataSet[0] * extentDataSet;

    MPI_Type_get_extent(MPI_INT, &lbDataSet, &extentDataSet);
    offsetsDatatSet[2] = offsetsDatatSet[1] + blocksDataSet[1] * extentDataSet;

    MPI_Type_create_struct(3, blocksDataSet, offsetsDatatSet, typesDataSet, &dataSetType);
    MPI_Type_commit(&dataSetType);

    //-----------------------------------------------------------------------

    MPI_Type_contiguous(4, MPI_INT, &dataSetSmallType);
    MPI_Type_commit(&dataSetSmallType);

    //-----------------------------------------------------------------------

    MPI_Type_contiguous(4, MPI_INT, &boundCondParamType);
    MPI_Type_commit(&boundCondParamType);

    //---------------------------------------

    MPI_Type_contiguous(BLOCK_SIZE, boundCondType, &boundCondType1000);
    MPI_Type_commit(&boundCondType1000);

    //---------------------------------------

    MPI_Type_contiguous(6, MPI_INT, &boundCondTypeAdd);
    MPI_Type_commit(&boundCondTypeAdd);
}
//////////////////////////////////////////////////////////////////////////
MPIIORestartCoProcessor::~MPIIORestartCoProcessor()
{
    MPI_Type_free(&dataSetType);
    MPI_Type_free(&dataSetSmallType);
    MPI_Type_free(&boundCondParamType);
    MPI_Type_free(&boundCondType1000);
    MPI_Type_free(&boundCondTypeAdd);
}

//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::process(double step)
{
    if (scheduler->isDue(step)) {
        if (comm->isRoot())
            UBLOG(logINFO, "MPIIORestartCoProcessor save step: " << step);
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
//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::clearAllFiles(int step)
{
    MPI_File file_handler;
    MPI_Info info       = MPI_INFO_NULL;
    MPI_Offset new_size = 0;

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

    MPIIOCoProcessor::clearAllFiles(step);

    std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc10 =
        MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc10 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename10);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::writeBlocks(int step) { MPIIOCoProcessor::writeBlocks(step); }

void MPIIORestartCoProcessor::writeDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;
    DataSetRestart *dataSetArray = new DataSetRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    DSArraysPresence arrPresence;
    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<D3Q27EsoTwist3DSplittedVector> D3Q27EsoTwist3DSplittedVectorPtr;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetArray[ic].x2              = block->getX2();
            dataSetArray[ic].x3              = block->getX3();
            dataSetArray[ic].level           = block->getLevel();
            dataSetArray[ic].ghostLayerWidth = block->getKernel()->getGhostLayerWidth();
            dataSetArray[ic].collFactor      = block->getKernel()->getCollisionFactor();
            dataSetArray[ic].deltaT          = block->getKernel()->getDeltaT();
            dataSetArray[ic].compressible    = block->getKernel()->getCompressible();
            dataSetArray[ic].withForcing     = block->getKernel()->getWithForcing();

            D3Q27EsoTwist3DSplittedVectorPtr = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(
                block->getKernel()->getDataSet()->getFdistributions());
            localDistributions    = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
            nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
            zeroDistributions     = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();

            if (firstBlock) // when first (any) valid block...
            {
                if (localDistributions) {
                    dataSetParamStr1.nx[0] = static_cast<int>(localDistributions->getNX1());
                    dataSetParamStr1.nx[1] = static_cast<int>(localDistributions->getNX2());
                    dataSetParamStr1.nx[2] = static_cast<int>(localDistributions->getNX3());
                    dataSetParamStr1.nx[3] = static_cast<int>(localDistributions->getNX4());
                }
                if (nonLocalDistributions) {
                    dataSetParamStr2.nx[0] = static_cast<int>(nonLocalDistributions->getNX1());
                    dataSetParamStr2.nx[1] = static_cast<int>(nonLocalDistributions->getNX2());
                    dataSetParamStr2.nx[2] = static_cast<int>(nonLocalDistributions->getNX3());
                    dataSetParamStr2.nx[3] = static_cast<int>(nonLocalDistributions->getNX4());
                }
                if (zeroDistributions) {
                    dataSetParamStr3.nx[0] = static_cast<int>(zeroDistributions->getNX1());
                    dataSetParamStr3.nx[1] = static_cast<int>(zeroDistributions->getNX2());
                    dataSetParamStr3.nx[2] = static_cast<int>(zeroDistributions->getNX3());
                    dataSetParamStr3.nx[3] = 1;
                }

                // ... than save some parameters that are equal in all dataSets
                dataSetParamStr1.nx1 = dataSetParamStr2.nx1 = dataSetParamStr3.nx1 =
                    static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX1());
                dataSetParamStr1.nx2 = dataSetParamStr2.nx2 = dataSetParamStr3.nx2 =
                    static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX2());
                dataSetParamStr1.nx3 = dataSetParamStr2.nx3 = dataSetParamStr3.nx3 =
                    static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX3());

                doubleCountInBlock =
                    dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
                    dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
                    dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> averageDensityArray =
                    block->getKernel()->getDataSet()->getAverageDensity();
                if (averageDensityArray)
                    arrPresence.isAverageDensityArrayPresent = true;
                else
                    arrPresence.isAverageDensityArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageVelocityArray3DPtr =
                    block->getKernel()->getDataSet()->getAverageVelocity();
                if (AverageVelocityArray3DPtr)
                    arrPresence.isAverageVelocityArrayPresent = true;
                else
                    arrPresence.isAverageVelocityArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageFluctArray3DPtr =
                    block->getKernel()->getDataSet()->getAverageFluctuations();
                if (AverageFluctArray3DPtr)
                    arrPresence.isAverageFluktuationsArrayPresent = true;
                else
                    arrPresence.isAverageFluktuationsArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageTripleArray3DPtr =
                    block->getKernel()->getDataSet()->getAverageTriplecorrelations();
                if (AverageTripleArray3DPtr)
                    arrPresence.isAverageTripleArrayPresent = true;
                else
                    arrPresence.isAverageTripleArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ShearStressValArray3DPtr =
                    block->getKernel()->getDataSet()->getShearStressValues();
                if (ShearStressValArray3DPtr)
                    arrPresence.isShearStressValArrayPresent = true;
                else
                    arrPresence.isShearStressValArrayPresent = false;

                SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> relaxationFactor3DPtr =
                    block->getKernel()->getDataSet()->getRelaxationFactor();
                if (relaxationFactor3DPtr)
                    arrPresence.isRelaxationFactorPresent = true;
                else
                    arrPresence.isRelaxationFactorPresent = false;

                firstBlock = false;
            }

            if (localDistributions && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) &&
                (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), localDistributions->getDataVector().begin(),
                                         localDistributions->getDataVector().end());
            if (nonLocalDistributions && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) &&
                (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), nonLocalDistributions->getDataVector().begin(),
                                         nonLocalDistributions->getDataVector().end());
            if (zeroDistributions && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) &&
                (dataSetParamStr3.nx[2] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), zeroDistributions->getDataVector().begin(),
                                         zeroDistributions->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + 3 * sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + 3 * sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), &dataSetParamStr2, 1,
                      dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + 2 * sizeof(dataSetParam)), &dataSetParamStr3, 1,
                      dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + 3 * sizeof(dataSetParam)), dataSetArray, blocksCount,
                      dataSetType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + 3 * sizeof(dataSetParam) + blocksCount * sizeof(DataSetRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);

    MPI_Type_free(&dataSetDoubleType);

    delete[] dataSetArray;

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet time: " << finish - start << " s");
    }

    MPI_File file_handler1;
    std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler1);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename1);
    MPI_File_write_at(file_handler1, (MPI_Offset)0, &arrPresence, 1, arrayPresenceType, MPI_STATUS_IGNORE);
    MPI_File_sync(file_handler1);
    MPI_File_close(&file_handler1);

    if (arrPresence.isAverageDensityArrayPresent)
        writeAverageDensityArray(step);

    if (arrPresence.isAverageVelocityArrayPresent)
        writeAverageVelocityArray(step);

    if (arrPresence.isAverageFluktuationsArrayPresent)
        writeAverageFluktuationsArray(step);

    if (arrPresence.isAverageTripleArrayPresent)
        writeAverageTripleArray(step);

    if (arrPresence.isShearStressValArrayPresent)
        writeShearStressValArray(step);

    if (arrPresence.isRelaxationFactorPresent)
        writeRelaxationFactor(step);
}

void MPIIORestartCoProcessor::writeAverageDensityArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values of the AverageDensityArray in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageDensityArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> averageDensityArray;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            averageDensityArray = block->getKernel()->getDataSet()->getAverageDensity();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(averageDensityArray->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(averageDensityArray->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(averageDensityArray->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(averageDensityArray->getNX4());
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) &&
                (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), averageDensityArray->getDataVector().begin(),
                                         averageDensityArray->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageDensityArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageDensityArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeAverageVelocityArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageVelocityArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageVelocityArray3DPtr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) &&
                (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), AverageVelocityArray3DPtr->getDataVector().begin(),
                                         AverageVelocityArray3DPtr->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageVelocityArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);

    MPI_Type_free(&dataSetDoubleType);
    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageVelocityArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeAverageFluktuationsArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageFluktuationsArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageFluctArray3DPtr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) &&
                (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), AverageFluctArray3DPtr->getDataVector().begin(),
                                         AverageFluctArray3DPtr->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageFluktuationsArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename =
        path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageFluktuationsArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeAverageTripleArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageTripleArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageTripleArray3DPtr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) &&
                (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), AverageTripleArray3DPtr->getDataVector().begin(),
                                         AverageTripleArray3DPtr->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageTripleArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeAverageTripleArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeShearStressValArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeShearStressValArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ShearStressValArray3DPtr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) &&
                (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), ShearStressValArray3DPtr->getDataVector().begin(),
                                         ShearStressValArray3DPtr->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeShearStressValArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeShearStressValArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeRelaxationFactor(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeRelaxationFactor start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> RelaxationFactor3DPtr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            dataSetSmallArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            dataSetSmallArray[ic].x2    = block->getX2();
            dataSetSmallArray[ic].x3    = block->getX3();
            dataSetSmallArray[ic].level = block->getLevel();

            RelaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();

            if (firstBlock) // when first (any) valid block...
            {
                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(RelaxationFactor3DPtr->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(RelaxationFactor3DPtr->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(RelaxationFactor3DPtr->getNX3());
                dataSetParamStr.nx[3] = 1;
                doubleCountInBlock =
                    dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), RelaxationFactor3DPtr->getDataVector().begin(),
                                         RelaxationFactor3DPtr->getDataVector().end());

            ic++;
        }
    }

    // register new MPI-types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeRelaxationFactor start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + sizeof(dataSetParam) +
                                blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    // each process writes data identifying blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                      dataSetSmallType, MPI_STATUS_IGNORE);
    // each process writes the dataSet arrays
    if (doubleValuesArray.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                          &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeRelaxationFactor time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::writeBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    int blocksCount          = 0; // quantity of blocks in the grid, max 2147483648 blocks!
    size_t count_boundCond   = 0; // how many BoundaryConditions in all blocks
    int count_indexContainer = 0; // how many indexContainer-values in all blocks
    size_t byteCount         = 0; // how many bytes writes this process in the file

    std::vector<SPtr<Block3D>> blocksVector[25];
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, rank, blocksVector[level]);
        blocksCount += static_cast<int>(blocksVector[level].size());
    }

    BCAddRestart *bcAddArray = new BCAddRestart[blocksCount];
    std::vector<BoundaryCondition> bcVector;
    std::vector<int> bcindexmatrixV;
    std::vector<int> indexContainerV;
    bool bcindexmatrixCountNotInit = true;
    int ic                         = 0;
    SPtr<BCArray3D> bcArr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blocksVector[level]) // all the blocks of the current level
        {
            bcArr = block->getKernel()->getBCProcessor()->getBCArray();

            bcAddArray[ic].x1 =
                block->getX1(); // coordinates of the block needed to find it while regenerating the grid
            bcAddArray[ic].x2                   = block->getX2();
            bcAddArray[ic].x3                   = block->getX3();
            bcAddArray[ic].level                = block->getLevel();
            bcAddArray[ic].boundCond_count      = 0; // how many BoundaryConditions in this block
            bcAddArray[ic].indexContainer_count = 0; // how many indexContainer-values in this block

            for (std::size_t bc = 0; bc < bcArr->getBCVectorSize(); bc++) {
                BoundaryCondition *bouCond = new BoundaryCondition();
                if (bcArr->bcvector[bc] == NULL) {
                    memset(bouCond, 0, sizeof(BoundaryCondition));
                } else {
                    bouCond->noslipBoundaryFlags    = bcArr->bcvector[bc]->getNoSlipBoundary();
                    bouCond->slipBoundaryFlags      = bcArr->bcvector[bc]->getSlipBoundary();
                    bouCond->velocityBoundaryFlags  = bcArr->bcvector[bc]->getVelocityBoundary();
                    bouCond->densityBoundaryFlags   = bcArr->bcvector[bc]->getDensityBoundary();
                    bouCond->wallModelBoundaryFlags = bcArr->bcvector[bc]->getWallModelBoundary();
                    bouCond->bcVelocityX1           = (float)bcArr->bcvector[bc]->getBoundaryVelocityX1();
                    bouCond->bcVelocityX2           = (float)bcArr->bcvector[bc]->getBoundaryVelocityX2();
                    bouCond->bcVelocityX3           = (float)bcArr->bcvector[bc]->getBoundaryVelocityX3();
                    bouCond->bcDensity              = (float)bcArr->bcvector[bc]->getBoundaryDensity();
                    bouCond->bcLodiDensity          = (float)bcArr->bcvector[bc]->getDensityLodiDensity();
                    bouCond->bcLodiVelocityX1       = (float)bcArr->bcvector[bc]->getDensityLodiVelocityX1();
                    bouCond->bcLodiVelocityX2       = (float)bcArr->bcvector[bc]->getDensityLodiVelocityX2();
                    bouCond->bcLodiVelocityX3       = (float)bcArr->bcvector[bc]->getDensityLodiVelocityX3();
                    bouCond->bcLodiLentgh           = (float)bcArr->bcvector[bc]->getDensityLodiLength();
                    bouCond->nx1                    = (float)bcArr->bcvector[bc]->nx1;
                    bouCond->nx2                    = (float)bcArr->bcvector[bc]->nx2;
                    bouCond->nx3                    = (float)bcArr->bcvector[bc]->nx3;
                    for (int iq = 0; iq < 26; iq++)
                        bouCond->q[iq] = (float)bcArr->bcvector[bc]->getQ(iq);
                    bouCond->algorithmType = bcArr->bcvector[bc]->getBcAlgorithmType();
                }

                bcVector.push_back(*bouCond);
                bcAddArray[ic].boundCond_count++;
                count_boundCond++;
            }

            // the quantity of elements in the bcindexmatrix array (CbArray3D<int, IndexerX3X2X1>) in bcArray(BCArray3D)
            // is always equal, this will be the size of the "write-read-block" in MPI_write_.../MPI_read-functions when
            // writing/reading BoundConds
            if (bcindexmatrixCountNotInit) {
                boundCondParamStr.nx1                = static_cast<int>(bcArr->bcindexmatrix.getNX1());
                boundCondParamStr.nx2                = static_cast<int>(bcArr->bcindexmatrix.getNX2());
                boundCondParamStr.nx3                = static_cast<int>(bcArr->bcindexmatrix.getNX3());
                boundCondParamStr.bcindexmatrixCount = static_cast<int>(bcArr->bcindexmatrix.getDataVector().size());
                bcindexmatrixCountNotInit            = false;
            }
            bcindexmatrixV.insert(bcindexmatrixV.end(), bcArr->bcindexmatrix.getDataVector().begin(),
                                  bcArr->bcindexmatrix.getDataVector().end());

            indexContainerV.insert(indexContainerV.end(), bcArr->indexContainer.begin(), bcArr->indexContainer.end());
            bcAddArray[ic].indexContainer_count = static_cast<int>(bcArr->indexContainer.size());
            count_indexContainer += bcAddArray[ic].indexContainer_count;

            ic++;
        }
    }

    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    // how many "big blocks" of BLOCK_SIZE size can by formed
    int bcBlockCount = (int)(count_boundCond / BLOCK_SIZE);
    if (bcBlockCount * BLOCK_SIZE < (int)count_boundCond)
        bcBlockCount += 1;
    for (int i = (int)count_boundCond; i < bcBlockCount * BLOCK_SIZE; i++) {
        BoundaryCondition *bouCond = new BoundaryCondition();
        memset(bouCond, 0, sizeof(BoundaryCondition));
        bcVector.push_back(*bouCond);
    }

    byteCount = bcBlockCount * BLOCK_SIZE * sizeof(BoundaryCondition) + blocksCount * sizeof(BCAddRestart) +
                sizeof(int) * (blocksCount * boundCondParamStr.bcindexmatrixCount + count_indexContainer);

    // write to the file
    // all processes calculate their offsets (quantity of bytes that the process is going to write)
    // and notify the next process (with the rank = rank + 1)
    MPI_Offset write_offset  = (MPI_Offset)(size * (3 * sizeof(int) + sizeof(boundCondParam)));
    size_t next_write_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_write_offset = write_offset + byteCount;
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_write_offset = write_offset + byteCount;
            if (rank < size - 1)
                MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "40");
    MPI_Info_set(info, "striping_unit", "4M");
#endif

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset1 = (MPI_Offset)(rank * (3 * sizeof(int) + sizeof(boundCondParam)));

    // each process writes the quantity of it's blocks
    MPI_File_write_at(file_handler, write_offset1, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // each process writes the quantity of "big blocks" of BLOCK_SIZE of boundary conditions
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1 + sizeof(int)), &bcBlockCount, 1, MPI_INT,
                      MPI_STATUS_IGNORE);
    // each process writes the quantity of indexContainer elements in all blocks
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1 + 2 * sizeof(int)), &count_indexContainer, 1, MPI_INT,
                      MPI_STATUS_IGNORE);
    // each process writes the quantity of bcindexmatrix elements in every block
    MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1 + 3 * sizeof(int)), &boundCondParamStr, 1,
                      boundCondParamType, MPI_STATUS_IGNORE);

    // each process writes data identifying the blocks
    MPI_File_write_at(file_handler, write_offset, bcAddArray, blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
    // each process writes boundary conditions
    if (bcVector.size() > 0)
        MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + blocksCount * sizeof(BCAddRestart)), &bcVector[0],
                          bcBlockCount, boundCondType1000, MPI_STATUS_IGNORE);
    // each process writes bcindexmatrix values
    if (bcindexmatrixV.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + blocksCount * sizeof(BCAddRestart) +
                                       bcBlockCount * BLOCK_SIZE * sizeof(BoundaryCondition)),
                          &bcindexmatrixV[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
    // each process writes indexContainer values
    if (indexContainerV.size() > 0)
        MPI_File_write_at(file_handler,
                          (MPI_Offset)(write_offset + blocksCount * sizeof(BCAddRestart) +
                                       bcBlockCount * BLOCK_SIZE * sizeof(BoundaryCondition) +
                                       blocksCount * boundCondParamStr.bcindexmatrixCount * sizeof(int)),
                          &indexContainerV[0], count_indexContainer, MPI_INT, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds time: " << finish - start << " s");
    }

    delete[] bcAddArray;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIORestartCoProcessor::restart(int step)
{
    if (comm->isRoot())
        UBLOG(logINFO, "MPIIORestartCoProcessor restart step: " << step);
    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - start");

    readBlocks(step);
    readDataSet(step);
    readBoundaryConds(step);

    grid->setTimeStep(step);

    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - end");
}

void MPIIORestartCoProcessor::readBlocks(int step) { MPIIOCoProcessor::readBlocks(step); }

void MPIIORestartCoProcessor::readDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;

    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, read_offset, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), &dataSetParamStr2, 1,
                     dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + 2 * sizeof(dataSetParam)), &dataSetParamStr3, 1,
                     dataSetParamType, MPI_STATUS_IGNORE);

    DataSetRestart *dataSetArray = new DataSetRestart[blocksCount];
    double doubleCountInBlock =
        dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
        dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
        dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];
    std::vector<double> doubleValuesArray(size_t(blocksCount * doubleCountInBlock)); // double-values in all blocks

    //   define MPI_types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + 3 * sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetRestart) + size_t(doubleCountInBlock) * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + 3 * sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetRestart) + size_t(doubleCountInBlock) * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + 3 * sizeof(dataSetParam)), dataSetArray, blocksCount,
                     dataSetType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler,
                     (MPI_Offset)(read_offset + 3 * sizeof(dataSetParam) + blocksCount * sizeof(DataSetRestart)),
                     &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0, vectorSize = 0;
    std::vector<LBMReal> vectorsOfValues1, vectorsOfValues2, vectorsOfValues3;

    for (int n = 0; n < blocksCount; n++) {
        vectorSize = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3];
        vectorsOfValues1.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + vectorSize);
        index += vectorSize;

        vectorSize = dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3];
        vectorsOfValues2.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + vectorSize);
        index += vectorSize;

        vectorSize = dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];
        vectorsOfValues3.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + vectorSize);
        index += vectorSize;

        SPtr<DistributionArray3D> mFdistributions(new D3Q27EsoTwist3DSplittedVector());

        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)
            ->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues1, dataSetParamStr1.nx[0],
                                                        dataSetParamStr1.nx[1], dataSetParamStr1.nx[2],
                                                        dataSetParamStr1.nx[3])));
        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)
            ->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues2, dataSetParamStr2.nx[0],
                                                        dataSetParamStr2.nx[1], dataSetParamStr2.nx[2],
                                                        dataSetParamStr2.nx[3])));
        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)
            ->setZeroDistributions(
                CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                    vectorsOfValues3, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(dataSetParamStr1.nx1);
        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(dataSetParamStr1.nx2);
        dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(dataSetParamStr1.nx3);

        // find the nesessary block and fill it
        SPtr<Block3D> block =
            grid->getBlock(dataSetArray[n].x1, dataSetArray[n].x2, dataSetArray[n].x3, dataSetArray[n].level);
        this->lbmKernel->setBlock(block);
        SPtr<LBMKernel> kernel = this->lbmKernel->clone();
        kernel->setGhostLayerWidth(dataSetArray[n].ghostLayerWidth);
        kernel->setCollisionFactor(dataSetArray[n].collFactor);
        kernel->setDeltaT(dataSetArray[n].deltaT);
        kernel->setCompressible(dataSetArray[n].compressible);
        kernel->setWithForcing(dataSetArray[n].withForcing);
        SPtr<DataSet3D> dataSetPtr = SPtr<DataSet3D>(new DataSet3D());
        dataSetPtr->setFdistributions(mFdistributions);
        kernel->setDataSet(dataSetPtr);
        block->setKernel(kernel);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
        readAverageDensityArray(step);

    if (arrPresence.isAverageVelocityArrayPresent)
        readAverageVelocityArray(step);

    if (arrPresence.isAverageFluktuationsArrayPresent)
        readAverageFluktuationsArray(step);

    if (arrPresence.isAverageTripleArrayPresent)
        readAverageTripleArray(step);

    if (arrPresence.isShearStressValArrayPresent)
        readShearStressValArray(step);

    if (arrPresence.isRelaxationFactorPresent)
        readRelaxationFactor(step);
}

void MPIIORestartCoProcessor::readAverageDensityArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageDensityArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    memset(&dataSetParamStr, 0, sizeof(dataSetParam));

    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageDensityArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageDensityArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill mAverageDensity arrays
        SPtr<AverageValuesArray3D> mAverageDensity;
        mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1],
                                                    dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setAverageDensity(mAverageDensity);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageDensityArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readAverageVelocityArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageVelocityArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageVelocityArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageVelocityArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill mAverageVelocity array
        SPtr<AverageValuesArray3D> mAverageVelocity;
        mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1],
                                                    dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setAverageVelocity(mAverageVelocity);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageVelocityArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readAverageFluktuationsArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageFluktuationsArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename =
        path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageFluktuationsArray time: " << finish - start << " s");
        UBLOG(logINFO,
              "MPIIORestartCoProcessor::readAverageFluktuationsArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill AverageFluktuations array
        SPtr<AverageValuesArray3D> mAverageFluktuations;
        mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1],
                                                    dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setAverageFluctuations(mAverageFluktuations);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageFluktuationsArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readAverageTripleArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageTripleArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageTripleArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageTripleArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill AverageTriplecorrelations array
        SPtr<AverageValuesArray3D> mAverageTriplecorrelations;
        mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1],
                                                    dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setAverageTriplecorrelations(mAverageTriplecorrelations);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readAverageTripleArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readShearStressValArray(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readShearStressValArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readShearStressValArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readShearStressValArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill ShearStressValuesArray array
        SPtr<ShearStressValuesArray3D> mShearStressValues;
        mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
            new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1],
                                                    dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setShearStressValues(mShearStressValues);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readShearStressValArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readRelaxationFactor(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readRelaxationFactor start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // read count of blocks
    int blocksCount = 0;
    dataSetParam dataSetParamStr;
    MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(size * sizeof(int)), &dataSetParamStr, 1, dataSetParamType,
                     MPI_STATUS_IGNORE);

    DataSetSmallRestart *dataSetSmallArray = new DataSetSmallRestart[blocksCount];
    int doubleCountInBlock =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    // calculate the read offset
    MPI_Offset read_offset  = (MPI_Offset)(size * sizeof(int));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + sizeof(dataSetParam) +
                               blocksCount * (sizeof(DataSetSmallRestart) + doubleCountInBlock * sizeof(double));
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallArray, blocksCount,
                     dataSetSmallType, MPI_STATUS_IGNORE);
    if (doubleCountInBlock > 0)
        MPI_File_read_at(file_handler,
                         (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRestart)),
                         &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readRelaxationFactor time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readRelaxationFactor start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    size_t nextVectorSize =
        dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<LBMReal> vectorsOfValues;
    for (int n = 0; n < blocksCount; n++) {
        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
        index += nextVectorSize;

        // fill RelaxationFactor array
        SPtr<RelaxationFactorArray3D> mRelaxationFactor;
        mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].x1, dataSetSmallArray[n].x2, dataSetSmallArray[n].x3,
                                             dataSetSmallArray[n].level);
        block->getKernel()->getDataSet()->setRelaxationFactor(mRelaxationFactor);
    }

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readRelaxationFactor end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

void MPIIORestartCoProcessor::readBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    int blocksCount         = 0;
    int dataCount1000       = 0;
    int dataCount2          = 0;
    MPI_Offset read_offset1 = (MPI_Offset)(rank * (3 * sizeof(int) + sizeof(boundCondParam)));

    // read count of blocks
    MPI_File_read_at(file_handler, read_offset1, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    // read count of big BoundaryCondition blocks
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1 + sizeof(int)), &dataCount1000, 1, MPI_INT,
                     MPI_STATUS_IGNORE);
    // read count of indexContainer values in all blocks
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1 + 2 * sizeof(int)), &dataCount2, 1, MPI_INT,
                     MPI_STATUS_IGNORE);
    // read count of bcindexmatrix values in every block
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1 + 3 * sizeof(int)), &boundCondParamStr, 1,
                     boundCondParamType, MPI_STATUS_IGNORE);

    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    size_t dataCount               = dataCount1000 * BLOCK_SIZE;
    BCAddRestart *bcAddArray       = new BCAddRestart[blocksCount];
    BoundaryCondition *bcArray     = new BoundaryCondition[dataCount];
    BoundaryCondition *nullBouCond = new BoundaryCondition();
    memset(nullBouCond, 0, sizeof(BoundaryCondition));
    int *intArray1 = new int[blocksCount * boundCondParamStr.bcindexmatrixCount];
    int *intArray2 = new int[dataCount2];

    MPI_Offset read_offset  = (MPI_Offset)(size * (3 * sizeof(int) + sizeof(boundCondParam)));
    size_t next_read_offset = 0;

    if (size > 1) {
        if (rank == 0) {
            next_read_offset = read_offset + blocksCount * sizeof(BCAddRestart) +
                               dataCount * sizeof(BoundaryCondition) +
                               (blocksCount * boundCondParamStr.bcindexmatrixCount + dataCount2) * sizeof(int);
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            next_read_offset = read_offset + blocksCount * sizeof(BCAddRestart) +
                               dataCount * sizeof(BoundaryCondition) +
                               (blocksCount * boundCondParamStr.bcindexmatrixCount + dataCount2) * sizeof(int);
            if (rank < size - 1)
                MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_File_read_at(file_handler, read_offset, bcAddArray, blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + blocksCount * sizeof(BCAddRestart)), &bcArray[0],
                     dataCount1000, boundCondType1000, MPI_STATUS_IGNORE);
    MPI_File_read_at(
        file_handler,
        (MPI_Offset)(read_offset + blocksCount * sizeof(BCAddRestart) + dataCount * sizeof(BoundaryCondition)),
        &intArray1[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
    MPI_File_read_at(file_handler,
                     (MPI_Offset)(read_offset + blocksCount * sizeof(BCAddRestart) +
                                  dataCount * sizeof(BoundaryCondition) +
                                  blocksCount * boundCondParamStr.bcindexmatrixCount * sizeof(int)),
                     &intArray2[0], dataCount2, MPI_INT, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    if (comm->isRoot()) {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    int index = 0, index1 = 0, index2 = 0;
    std::vector<SPtr<BoundaryConditions>> bcVector;
    std::vector<int> bcindexmatrixV;
    std::vector<int> indexContainerV;

    for (int n = 0; n < blocksCount; n++) {
        bcVector.resize(0);
        bcindexmatrixV.resize(0);
        indexContainerV.resize(0);

        for (int ibc = 0; ibc < bcAddArray[n].boundCond_count; ibc++) {
            SPtr<BoundaryConditions> bc;
            if (memcmp(&bcArray[index], nullBouCond, sizeof(BoundaryCondition)) == 0)
                bc = SPtr<BoundaryConditions>();
            else {
                bc                         = SPtr<BoundaryConditions>(new BoundaryConditions);
                bc->noslipBoundaryFlags    = bcArray[index].noslipBoundaryFlags;
                bc->slipBoundaryFlags      = bcArray[index].slipBoundaryFlags;
                bc->densityBoundaryFlags   = bcArray[index].densityBoundaryFlags;
                bc->velocityBoundaryFlags  = bcArray[index].velocityBoundaryFlags;
                bc->wallModelBoundaryFlags = bcArray[index].wallModelBoundaryFlags;
                bc->bcVelocityX1           = bcArray[index].bcVelocityX1;
                bc->bcVelocityX2           = bcArray[index].bcVelocityX2;
                bc->bcVelocityX3           = bcArray[index].bcVelocityX3;
                bc->bcDensity              = bcArray[index].bcDensity;
                bc->bcLodiDensity          = bcArray[index].bcLodiDensity;
                bc->bcLodiVelocityX1       = bcArray[index].bcLodiVelocityX1;
                bc->bcLodiVelocityX2       = bcArray[index].bcLodiVelocityX2;
                bc->bcLodiVelocityX3       = bcArray[index].bcLodiVelocityX3;
                bc->bcLodiLentgh           = bcArray[index].bcLodiLentgh;

                bc->nx1 = bcArray[index].nx1;
                bc->nx2 = bcArray[index].nx2;
                bc->nx3 = bcArray[index].nx3;
                for (int iq = 0; iq < 26; iq++)
                    bc->setQ(bcArray[index].q[iq], iq);
                bc->setBcAlgorithmType(bcArray[index].algorithmType);
            }

            bcVector.push_back(bc);
            index++;
        }

        for (int b1 = 0; b1 < boundCondParamStr.bcindexmatrixCount; b1++)
            bcindexmatrixV.push_back(intArray1[index1++]);

        for (int b2 = 0; b2 < bcAddArray[n].indexContainer_count; b2++)
            indexContainerV.push_back(intArray2[index2++]);

        CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2,
                                           boundCondParamStr.nx3);

        SPtr<Block3D> block = grid->getBlock(bcAddArray[n].x1, bcAddArray[n].x2, bcAddArray[n].x3, bcAddArray[n].level);
        SPtr<BCProcessor> bcProc = bcProcessor->clone(block->getKernel());
        SPtr<BCArray3D> bcArr(new BCArray3D());
        bcArr->bcindexmatrix  = bcim;
        bcArr->bcvector       = bcVector;
        bcArr->indexContainer = indexContainerV;
        bcProc->setBCArray(bcArr);

        block->getKernel()->setBCProcessor(bcProc);
    }

    delete nullBouCond;
    delete[] bcArray;
    delete[] bcAddArray;
    delete[] intArray1;
    delete[] intArray2;

    if (comm->isRoot()) {
        UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "
                           << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::setLBMKernel(SPtr<LBMKernel> kernel) { this->lbmKernel = kernel; }
//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::setBCProcessor(SPtr<BCProcessor> bcProcessor) { this->bcProcessor = bcProcessor; }
