#include "MPIIOMigrationBECoProcessor.h"
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
#include "LBMKernel.h"
#include "MetisPartitioningGridVisitor.h"
#include "PointerDefinitions.h"
#include "RenumberGridVisitor.h"
#include "UbFileInputASCII.h"
#include "UbFileOutputASCII.h"
#include "UbScheduler.h"
#include "WbWriter.h"
#include <MemoryUtil.h>
#include <UbSystem.h>

using namespace MPIIODataStructures;

#define MESSAGE_TAG 80
#define SEND_BLOCK_SIZE 100000

MPIIOMigrationBECoProcessor::MPIIOMigrationBECoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, SPtr<Communicator> comm)
    : MPIIOCoProcessor(grid, s, path, comm), nue(-999.999), nuL(-999.999), nuG(-999.999), densityRatio(-999.999)
{
    memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));

    //-------------------------   define MPI types  ---------------------------------

    MPI_Type_contiguous(SEND_BLOCK_SIZE, MPI_INT, &sendBlockIntType);
    MPI_Type_commit(&sendBlockIntType);
}

//////////////////////////////////////////////////////////////////////////
MPIIOMigrationBECoProcessor::~MPIIOMigrationBECoProcessor() { MPI_Type_free(&sendBlockIntType); }

void MPIIOMigrationBECoProcessor::process(double step)
{
    if (scheduler->isDue(step)) {
        if (comm->isRoot())
            UBLOG(logINFO, "MPIIOMigrationBECoProcessor save step: " << step);
        if (comm->isRoot())
            UBLOG(logINFO, "Save check point - start");
        clearAllFiles((int)step);

        writeBlocks((int)step);
        writeDataSet((int)step);
        writeBoundaryConds((int)step);

        writeCpTimeStep((int)step);

        if (comm->isRoot())
            UBLOG(logINFO, "Save check point - end");
    }
}

void MPIIOMigrationBECoProcessor::clearAllFiles(int step)
{
    MPI_File file_handler;
    MPI_Info info       = MPI_INFO_NULL;
    MPI_Offset new_size = 0;

    MPIIOCoProcessor::clearAllFiles(step);

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

    std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
    int rc10 =
        MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc10 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename10);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);

    std::string filename11 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC2.bin";
    int rc11 =
        MPI_File_open(MPI_COMM_WORLD, filename11.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc11 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename11);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);
}

void MPIIOMigrationBECoProcessor::writeBlocks(int step)
{
    grid->deleteBlockIDs();
    RenumberGridVisitor renumber(comm);
    grid->accept(renumber);

    MPIIOCoProcessor::writeBlocks(step);
}

void MPIIOMigrationBECoProcessor::writeDataSet(int step)
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
    int firstGlobalID;
    std::vector<double> doubleValuesArrayF; // double-values (arrays of f's) in all blocks  Fdistribution
    std::vector<double> doubleValuesArrayH1; // double-values (arrays of f's) in all blocks  H1distribution
    // std::vector<double> doubleValuesArrayH2; // double-values (arrays of f's) in all blocks  H2distribution

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool multiPhase = false;
    DSArraysPresence arrPresence;
    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<D3Q27EsoTwist3DSplittedVector> D3Q27EsoTwist3DSplittedVectorPtrF = 0, D3Q27EsoTwist3DSplittedVectorPtrH1 = 0, D3Q27EsoTwist3DSplittedVectorPtrH2 = 0;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF = 0, localDistributionsH1 = 0, localDistributionsH2 = 0;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF = 0, nonLocalDistributionsH1 = 0, nonLocalDistributionsH2 = 0;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsF = 0, zeroDistributionsH1 = 0, zeroDistributionsH2 = 0;
    
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
            D3Q27EsoTwist3DSplittedVectorPtrF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
            localDistributionsF    = D3Q27EsoTwist3DSplittedVectorPtrF->getLocalDistributions();
            nonLocalDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getNonLocalDistributions();
            zeroDistributionsF     = D3Q27EsoTwist3DSplittedVectorPtrF->getZeroDistributions();
 
            D3Q27EsoTwist3DSplittedVectorPtrH1 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getHdistributions());
            if (D3Q27EsoTwist3DSplittedVectorPtrH1 != 0)
            {
                multiPhase = true;
                localDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getLocalDistributions();
                nonLocalDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getNonLocalDistributions();
                zeroDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getZeroDistributions();
            }

            /*D3Q27EsoTwist3DSplittedVectorPtrH2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getH2distributions());
            if (D3Q27EsoTwist3DSplittedVectorPtrH2 != 0)
            {
                localDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getLocalDistributions();
                nonLocalDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getNonLocalDistributions();
                zeroDistributionsH2 = D3Q27EsoTwist3DSplittedVectorPtrH2->getZeroDistributions();
            }*/


            if (firstBlock) // && block->getKernel()) // when first (any) valid block...
            {
                firstGlobalID = block->getGlobalID(); // id of the block needed to find it while regenerating the grid

                if (localDistributionsF) {
                    dataSetParamStr1.nx[0] = static_cast<int>(localDistributionsF->getNX1());
                    dataSetParamStr1.nx[1] = static_cast<int>(localDistributionsF->getNX2());
                    dataSetParamStr1.nx[2] = static_cast<int>(localDistributionsF->getNX3());
                    dataSetParamStr1.nx[3] = static_cast<int>(localDistributionsF->getNX4());
                }

                if (nonLocalDistributionsF) {
                    dataSetParamStr2.nx[0] = static_cast<int>(nonLocalDistributionsF->getNX1());
                    dataSetParamStr2.nx[1] = static_cast<int>(nonLocalDistributionsF->getNX2());
                    dataSetParamStr2.nx[2] = static_cast<int>(nonLocalDistributionsF->getNX3());
                    dataSetParamStr2.nx[3] = static_cast<int>(nonLocalDistributionsF->getNX4());
                }
                if (zeroDistributionsF) {
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

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> averageDensityArray = block->getKernel()->getDataSet()->getAverageDensity();
                if (averageDensityArray)
                    arrPresence.isAverageDensityArrayPresent = true;
                else
                    arrPresence.isAverageDensityArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
                if (AverageVelocityArray3DPtr)
                    arrPresence.isAverageVelocityArrayPresent = true;
                else
                    arrPresence.isAverageVelocityArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
                if (AverageFluctArray3DPtr)
                    arrPresence.isAverageFluktuationsArrayPresent = true;
                else
                    arrPresence.isAverageFluktuationsArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
                if (AverageTripleArray3DPtr)
                    arrPresence.isAverageTripleArrayPresent = true;
                else
                    arrPresence.isAverageTripleArrayPresent = false;

                SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
                if (ShearStressValArray3DPtr)
                    arrPresence.isShearStressValArrayPresent = true;
                else
                    arrPresence.isShearStressValArrayPresent = false;

                SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
                if (relaxationFactor3DPtr)
                    arrPresence.isRelaxationFactorPresent = true;
                else
                    arrPresence.isRelaxationFactorPresent = false;

                SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> phaseField3DPtr1 = block->getKernel()->getDataSet()->getPhaseField();
                if (phaseField3DPtr1)
                    arrPresence.isPhaseField1Present = true;
                else
                    arrPresence.isPhaseField1Present = false;

                SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> phaseField3DPtr2 = block->getKernel()->getDataSet()->getPhaseField2();
                if (phaseField3DPtr2)
                    arrPresence.isPhaseField2Present = true;
                else
                    arrPresence.isPhaseField2Present = false;


                firstBlock = false;
            }

            if (localDistributionsF && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), localDistributionsF->getDataVector().begin(), localDistributionsF->getDataVector().end());
            if (nonLocalDistributionsF && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), nonLocalDistributionsF->getDataVector().begin(), nonLocalDistributionsF->getDataVector().end());
            if (zeroDistributionsF && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                doubleValuesArrayF.insert(doubleValuesArrayF.end(), zeroDistributionsF->getDataVector().begin(), zeroDistributionsF->getDataVector().end());

            if (multiPhase)
            {
                if (localDistributionsH1 && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), localDistributionsH1->getDataVector().begin(), localDistributionsH1->getDataVector().end());
                if (nonLocalDistributionsH1 && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), nonLocalDistributionsH1->getDataVector().begin(), nonLocalDistributionsH1->getDataVector().end());
                if (zeroDistributionsH1 && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                    doubleValuesArrayH1.insert(doubleValuesArrayH1.end(), zeroDistributionsH1->getDataVector().begin(), zeroDistributionsH1->getDataVector().end());
            }

            /*if (D3Q27EsoTwist3DSplittedVectorPtrH2 != 0)
            {
                if (localDistributionsH2 && (dataSetParamStr1.nx[0] > 0) && (dataSetParamStr1.nx[1] > 0) && (dataSetParamStr1.nx[2] > 0) && (dataSetParamStr1.nx[3] > 0))
                doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), localDistributionsH2->getDataVector().begin(), localDistributionsH2->getDataVector().end());
                if (nonLocalDistributionsH2 && (dataSetParamStr2.nx[0] > 0) && (dataSetParamStr2.nx[1] > 0) && (dataSetParamStr2.nx[2] > 0) && (dataSetParamStr2.nx[3] > 0))
                doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), nonLocalDistributionsH2->getDataVector().begin(), nonLocalDistributionsH2->getDataVector().end());
                if (zeroDistributionsH2 && (dataSetParamStr3.nx[0] > 0) && (dataSetParamStr3.nx[1] > 0) && (dataSetParamStr3.nx[2] > 0))
                doubleValuesArrayH2.insert(doubleValuesArrayH2.end(), zeroDistributionsH2->getDataVector().begin(), zeroDistributionsH2->getDataVector().end());
            }*/

            ic++;
        }
    }

    MPI_Type_contiguous(doubleCountInBlock , MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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

    // write to the file
    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetF.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset = (MPI_Offset)(3 * sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

    MPI_File_write_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, write_offset, &doubleValuesArrayF[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);

    //-------------------------------- H1 ------------------------------------------------
    if (multiPhase)
    {
        filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
        rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
        if (rc != MPI_SUCCESS)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);

        write_offset = (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
        MPI_File_write_at(file_handler, write_offset, &doubleValuesArrayH1[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

        MPI_File_sync(file_handler);
        MPI_File_close(&file_handler);
    }

    //-------------------------------- H2 --------------------------------------------------
    /*if (D3Q27EsoTwist3DSplittedVectorPtr2 != 0)
    {
        filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH2.bin";
        rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
        if (rc != MPI_SUCCESS)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);

        write_offset = (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
        MPI_File_write_at(file_handler, write_offset, &doubleValuesArrayH2[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

        MPI_File_sync(file_handler);
        MPI_File_close(&file_handler);
    }    */

    //--------------------------------

    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet time: " << finish - start << " s");
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
        write4DArray(step, AverageDensity, std::string("/cpAverageDensityArray.bin"));
    // writeAverageDensityArray(step);

    if (arrPresence.isAverageVelocityArrayPresent)
        write4DArray(step, AverageVelocity, std::string("/cpAverageVelocityArray.bin"));
    // writeAverageVelocityArray(step);

    if (arrPresence.isAverageFluktuationsArrayPresent)
        write4DArray(step, AverageFluktuations, std::string("/cpAverageFluktuationsArray.bin"));
    // writeAverageFluktuationsArray(step);

    if (arrPresence.isAverageTripleArrayPresent)
        write4DArray(step, AverageTriple, std::string("/cpAverageTripleArray.bin"));
    // writeAverageTripleArray(step);

    if (arrPresence.isShearStressValArrayPresent)
        write4DArray(step, ShearStressVal, std::string("/cpShearStressValArray.bin"));
    // writeShearStressValArray(step);

    if (arrPresence.isRelaxationFactorPresent)
        write3DArray(step, RelaxationFactor, std::string("/cpRelaxationFactor.bin"));
    // writeRelaxationFactor(step);

    if (arrPresence.isPhaseField1Present)
        write3DArray(step, PhaseField1, std::string("/cpPhaseField1.bin"));

    if (arrPresence.isPhaseField2Present)
        write3DArray(step, PhaseField2, std::string("/cpPhaseField2.bin"));
    }

void MPIIOMigrationBECoProcessor::write4DArray(int step, Arrays arrayType, std::string fname)
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

    int firstGlobalID;
    std::vector<double> doubleValuesArray; // double-values of the data array in all blocks
    dataSetParam dataSetParamStr;
    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ___Array;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write4DArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
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
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationBECoProcessor::write4DArray : 4D array type does not exist!"));
                    break;
            }

            if (firstBlock) // when first (any) valid block...
            {
                firstGlobalID = block->getGlobalID();

                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(___Array->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(___Array->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(___Array->getNX3());
                dataSetParamStr.nx[3] = static_cast<int>(___Array->getNX4());
                doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) && (dataSetParamStr.nx[3] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), ___Array->getDataVector().begin(), ___Array->getDataVector().end());

            ic++;
        }
    }

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write4DArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // register new MPI-type depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write4DArray time: " << finish - start << " s");
    }
}

void MPIIOMigrationBECoProcessor::write3DArray(int step, Arrays arrayType, std::string fname)
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

    int firstGlobalID;
    std::vector<double> doubleValuesArray; // double-values of the data array in all blocks
    dataSetParam dataSetParamStr;
    bool firstBlock        = true;
    int doubleCountInBlock = 0;
    int ic                 = 0;
    SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> ___Array;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write3DArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
        {
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
                default:
                    UB_THROW(UbException(UB_EXARGS,
                    "MPIIOMigrationBECoProcessor::write3DArray : 3D array type does not exist!"));
                    break;
            }

            if (firstBlock) // when first (any) valid block...
            {
                firstGlobalID = block->getGlobalID();

                dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
                dataSetParamStr.nx[0] = static_cast<int>(___Array->getNX1());
                dataSetParamStr.nx[1] = static_cast<int>(___Array->getNX2());
                dataSetParamStr.nx[2] = static_cast<int>(___Array->getNX3());
                dataSetParamStr.nx[3] = 1;
                doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

                firstBlock = false;
            }

            if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0))
                doubleValuesArray.insert(doubleValuesArray.end(), ___Array->getDataVector().begin(), ___Array->getDataVector().end());

            ic++;
        }
    }

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write3DArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    // register new MPI-type depending on the block-specific information
    MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::write3DArray time: " << finish - start << " s");
    }
}

void MPIIOMigrationBECoProcessor::writeBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot())
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBoundaryConds start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    int blocksCount          = 0; // quantity of blocks, that belong to this process
    size_t allBytesCount     = 0; // quantity of bytes, that one process writes to the file
    size_t count_boundCond   = 0; // how many BoundaryConditions in all blocks
    int count_indexContainer = 0; // how many indexContainer-values in all blocks

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
    std::vector<int> *indexContainerVector   = new std::vector<int>[blocksCount];
    std::vector<int> bcindexmatrixVector;

    bool bcindexmatrixCountNotInit = true;
    int ic                         = 0;
    SPtr<BCArray3D> bcArr;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) // all the blocks of the current level
        {
            bcArr = block->getKernel()->getBCProcessor()->getBCArray();

            bcAddArray[ic].globalID = block->getGlobalID();                // id of the block needed to find it while regenerating the grid
            bcAddArray[ic].boundCond_count      = 0; // how many BoundaryConditions in this block
            bcAddArray[ic].indexContainer_count = 0; // how many indexContainer-values in this block
            bytesCount[ic]                      = sizeof(BCAddMigration);
            bcVector[ic].resize(0);
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

                bcVector[ic].push_back(*bouCond);
                bcAddArray[ic].boundCond_count++;
                count_boundCond++;
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

            bcindexmatrixVector.insert(bcindexmatrixVector.end(), bcArr->bcindexmatrix.getDataVector().begin(), bcArr->bcindexmatrix.getDataVector().end());

            indexContainerVector[ic].insert(indexContainerVector[ic].begin(), bcArr->indexContainer.begin(), bcArr->indexContainer.end());
            bcAddArray[ic].indexContainer_count = static_cast<int>(bcArr->indexContainer.size());
            count_indexContainer += bcAddArray[ic].indexContainer_count;
            bytesCount[ic] += bcAddArray[ic].indexContainer_count * sizeof(int);

            allBytesCount += bytesCount[ic];

            ic++;
        }
    }

    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    MPI_Info info = MPI_INFO_NULL;
    // MPI_Info_create (&info);
    // MPI_Info_set(info,"romio_cb_write","enable");
    // MPI_Info_set(info,"cb_buffer_size","4194304");
    // MPI_Info_set(info,"striping_unit","4194304");

    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_Offset write_offset = (MPI_Offset)(sizeof(int)) + (MPI_Offset)(bcAddArray[0].globalID) * (MPI_Offset)(boundCondParamStr.bcindexmatrixCount) * (MPI_Offset)(sizeof(int));

    MPI_File_write_at(file_handler, 0, &boundCondParamStr.bcindexmatrixCount, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, write_offset, &bcindexmatrixVector[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC2.bin";
    rc       = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_write_at(file_handler, 0, &boundCondParamStr, 4, MPI_INT, MPI_STATUS_IGNORE);

    write_offset = (MPI_Offset)(sizeof(boundCondParam)) + (MPI_Offset)(grid->getNumberOfBlocks()) * (MPI_Offset)(sizeof(size_t));
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

    MPI_Offset write_offsetIndex;

    for (int nb = 0; nb < blocksCount; nb++) 
    {
        write_offsetIndex = (MPI_Offset)(sizeof(boundCondParam)) + (MPI_Offset)(bcAddArray[nb].globalID) * (MPI_Offset)(sizeof(size_t));
        MPI_File_write_at(file_handler, write_offsetIndex, &write_offset, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);

        MPI_File_write_at(file_handler, write_offset, &bcAddArray[nb], 3, MPI_INT, MPI_STATUS_IGNORE);
        if (bcVector[nb].size() > 0)
            MPI_File_write_at(file_handler, write_offset + (MPI_Offset)(sizeof(BCAddMigration)), &bcVector[nb][0], bcAddArray[nb].boundCond_count, boundCondType, MPI_STATUS_IGNORE);

        if (indexContainerVector[nb].size() > 0)
            MPI_File_write_at(file_handler, write_offset + (MPI_Offset)(sizeof(BCAddMigration)) + (MPI_Offset)(bcAddArray[nb].boundCond_count) * (MPI_Offset)(sizeof(BoundaryCondition)),
                &indexContainerVector[nb][0], bcAddArray[nb].indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);

        write_offset += bytesCount[nb];
    }

    MPI_File_sync(file_handler);
    MPI_File_close(&file_handler);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBoundaryConds time: " << finish - start << " s");
    }

    delete[] bcAddArray;
    delete[] bytesCount;
    delete[] bcVector;
    delete[] indexContainerVector;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIOMigrationBECoProcessor::restart(int step)
{
    if (comm->isRoot())
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor restart step: " << step);
    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - start");

    readBlocks(step);
    SPtr<Grid3DVisitor> newMetisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
    grid->accept(newMetisVisitor);

    readDataSet(step);
    readBoundaryConds(step);

    grid->setTimeStep(step);
    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - end");
}

void MPIIOMigrationBECoProcessor::readBlocks(int step) { MPIIOCoProcessor::readBlocks(step); }

void MPIIOMigrationBECoProcessor::blocksExchange(int tagN, int ind1, int ind2, int doubleCountInBlock, std::vector<double> &pV, std::vector<double> *rawDataReceive)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int indexB = ind1;
    int indexE = ind2;
    //   int myBlocksCount = indexE - indexB;

    int *blocksCounterSend = new int[size];
    int *blocksCounterRec  = new int[size];

    std::vector<double> *rawDataSend = new std::vector<double>[size];
    for (int r = 0; r < size; r++) 
    {
        rawDataSend[r].resize(0);
        blocksCounterSend[r] = 0;
        blocksCounterRec[r]  = 0;
    }

    SPtr<Block3D> tempBlock;
    int tempRank;

    for (int ind = 0; ind < indexE - indexB; ind++)
    {
        tempBlock = grid->getBlock(indexB + int(ind));
        if (!tempBlock)
            throw UbException(UB_EXARGS, "MPIIOMigrationBECoProcessor::blocksExchange -- null block pointer!!!");

        tempRank = tempBlock->getRank();

        if (tempRank == rank) // no need to send data, the process already has it
        {
            blocksCounterRec[tempRank]++;
            rawDataReceive[tempRank].push_back(double(indexB + ind));
            rawDataReceive[tempRank].insert(rawDataReceive[tempRank].end(), pV.begin() + ind * size_t(doubleCountInBlock),
                                            pV.begin() + ind * size_t(doubleCountInBlock) + size_t(doubleCountInBlock));
        } else // we must send data to other processes
        {
            blocksCounterSend[tempRank]++;
            rawDataSend[tempRank].push_back(double(indexB + ind));
            rawDataSend[tempRank].insert(rawDataSend[tempRank].end(), pV.begin() + ind * size_t(doubleCountInBlock),
                                         pV.begin() + ind * size_t(doubleCountInBlock) + size_t(doubleCountInBlock));
        }
    }

    MPI_Request *requests = new MPI_Request[size * 2]; // send + receive
    int requestCount      = 0;

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            MPI_Irecv(&blocksCounterRec[r], 1, MPI_INT, r, tagN, MPI_COMM_WORLD, &requests[requestCount]);
            requestCount++;
        }
    }

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            MPI_Isend(&blocksCounterSend[r], 1, MPI_INT, r, tagN, MPI_COMM_WORLD, &requests[requestCount]);
            requestCount++;
        }
    }

    MPI_Waitall(requestCount, &requests[0], MPI_STATUSES_IGNORE);

    MPI_Type_contiguous(doubleCountInBlock + 1, MPI_DOUBLE, &sendBlockDoubleType);
    MPI_Type_commit(&sendBlockDoubleType);

    for (int r = 0; r < size; r++) 
    {
        if (r != rank)
            rawDataReceive[r].resize(size_t(blocksCounterRec[r]) * size_t(doubleCountInBlock + 1));
    }

    requestCount         = 0;
    int sendRecCount     = 0;
    size_t sendRecOffset = 0;
    const int maxQuant   = 400;
    int restQuant;

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            sendRecCount = int(blocksCounterRec[r] / maxQuant);
            if (sendRecCount * maxQuant < blocksCounterRec[r])
                sendRecCount++;
            requests = (MPI_Request *)realloc(requests, (requestCount + sendRecCount) * sizeof(MPI_Request));

            for (int sc = 0; sc < sendRecCount; sc++)
            {
                restQuant     = (sc < sendRecCount - 1) ? maxQuant : blocksCounterRec[r] - sc * maxQuant;
                sendRecOffset = size_t(sc) * size_t(maxQuant) * size_t((doubleCountInBlock + 1));
                MPI_Irecv(&rawDataReceive[r][sendRecOffset], restQuant, sendBlockDoubleType, r, tagN, MPI_COMM_WORLD, &requests[requestCount]);
                requestCount++;
            }
        }
    }

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            sendRecCount = int(blocksCounterSend[r] / maxQuant);
            if (sendRecCount * maxQuant < blocksCounterSend[r])
                sendRecCount++;
            requests = (MPI_Request *)realloc(requests, (requestCount + sendRecCount) * sizeof(MPI_Request));

            for (int sc = 0; sc < sendRecCount; sc++) 
            {
                restQuant     = (sc < sendRecCount - 1) ? maxQuant : blocksCounterSend[r] - sc * maxQuant;
                sendRecOffset = size_t(sc) * size_t(maxQuant) * size_t((doubleCountInBlock + 1));
                MPI_Isend(&rawDataSend[r][sendRecOffset], restQuant, sendBlockDoubleType, r, tagN, MPI_COMM_WORLD,  &requests[requestCount]);
                requestCount++;
            }
        }
    }

    MPI_Waitall(requestCount, &requests[0], MPI_STATUSES_IGNORE);

    delete[] blocksCounterSend;
    delete[] blocksCounterRec;
    delete[] rawDataSend;
    delete[] requests;
}

void MPIIOMigrationBECoProcessor::readDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (!lbmKernel)
        UB_THROW(UbException(UB_EXARGS, "lbmKernel does not exist!"));
    if (!bcProcessor)
        UB_THROW(UbException(UB_EXARGS, "bcProcessor does not exist!"));
    if (nue == -999.999)
        UB_THROW(UbException(UB_EXARGS, "nue is not initialised!"));
    if (nuL == -999.999 )
        UB_THROW(UbException(UB_EXARGS, "nuL is not initialised!"));
    if (nuG == -999.999)
        UB_THROW(UbException(UB_EXARGS, "nuG is not initialised!"));
    if (densityRatio == -999.999)
        UB_THROW(UbException(UB_EXARGS, "densityRatio is not initialised!"));

    if (comm->isRoot())
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool multiPhase = false;
    dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;

    int blocksCountAll   = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
    int blocksPerProcess = blocksCountAll / size;     // how many blocks has each process

    size_t myBlocksCount;
    if (rank < (size - 1))
        myBlocksCount = blocksPerProcess;
    else
        myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

    int indexB = rank * blocksPerProcess;     // the first "my" block
    int indexE = indexB + int(myBlocksCount); // the latest "my" block

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

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
    std::vector<double> doubleValuesArrayF(size_t(myBlocksCount * doubleCountInBlock)); // double-values in all blocks  Fdistributions
    std::vector<double> doubleValuesArrayH1; // double-values in all blocks  H1distributions
    //std::vector<double> doubleValuesArrayH2; // double-values in all blocks  H2distributions

    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    //--------------------------------- F ---------------------------------------------------------
    MPI_Offset read_offset = (MPI_Offset)(3 * sizeof(dataSetParam)) + (MPI_Offset)(indexB * doubleCountInBlock * sizeof(double));
    MPI_File_read_at(file_handler, read_offset, &doubleValuesArrayF[0], int(myBlocksCount), dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handler);

    //--------------------------------- H1 ---------------------------------------------------------
    MPI_Offset fsize;
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);
    MPI_File_get_size(file_handler, &fsize);
    if (fsize > 0)
    {
        multiPhase = true;
        doubleValuesArrayH1.resize(myBlocksCount * doubleCountInBlock);

        read_offset = (MPI_Offset)(indexB * doubleCountInBlock * sizeof(double)) ;
        MPI_File_read_at(file_handler, read_offset, &doubleValuesArrayH1[0], int(myBlocksCount), dataSetDoubleType, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&file_handler);

    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start of exchange of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    std::vector<double>* rawDataReceiveF = new std::vector<double>[size];
    for (int r = 0; r < size; r++)
        rawDataReceiveF[r].resize(0);
    blocksExchange(MESSAGE_TAG, indexB, indexE, int(doubleCountInBlock), doubleValuesArrayF, rawDataReceiveF);

    std::vector<double>* rawDataReceiveH1 = new std::vector<double>[size];
    for (int r = 0; r < size; r++)
        rawDataReceiveH1[r].resize(0);
    blocksExchange(MESSAGE_TAG, indexB, indexE, int(doubleCountInBlock), doubleValuesArrayH1, rawDataReceiveH1);

    /*    std::vector<double>* rawDataReceiveH2 = new std::vector<double>[size];
        for (int r = 0; r < size; r++)
            rawDataReceiveH2[r].resize(0);
        blocksExchange(MESSAGE_TAG, indexB, indexE, int(doubleCountInBlock), doubleValuesArrayH2, rawDataReceiveH2);*/

    if (comm->isRoot())
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    
    //-------------------------------------- restore blocks ---------------------------------
    int blockID;
    std::vector<double> vectorsOfValuesF1, vectorsOfValuesF2, vectorsOfValuesF3;
    std::vector<double> vectorsOfValuesH11, vectorsOfValuesH12, vectorsOfValuesH13;
    //std::vector<double> vectorsOfValuesH21, vectorsOfValuesH22, vectorsOfValuesH23;

    size_t vectorSize1 = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3];
    size_t vectorSize2 = dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3];
    size_t vectorSize3 = dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

    size_t index;
    for (int r = 0; r < size; r++) 
    {
        index = 0;
        for (int ii = 0; ii < int(rawDataReceiveF[r].size() / doubleCountInBlock); ii++) 
        {
            blockID = (int)(rawDataReceiveF[r][index]);
            index += 1;

            vectorsOfValuesF1.assign(rawDataReceiveF[r].data() + index, rawDataReceiveF[r].data() + index + vectorSize1);
            if(multiPhase)
                vectorsOfValuesH11.assign(rawDataReceiveH1[r].data() + index, rawDataReceiveH1[r].data() + index + vectorSize1);
            //vectorsOfValuesH21.assign(rawDataReceiveH2[r].data() + index, rawDataReceiveH2[r].data() + index + vectorSize1);
            index += vectorSize1;

            vectorsOfValuesF2.assign(rawDataReceiveF[r].data() + index, rawDataReceiveF[r].data() + index + vectorSize2);
            if (multiPhase)
                vectorsOfValuesH12.assign(rawDataReceiveH1[r].data() + index, rawDataReceiveH1[r].data() + index + vectorSize1);
            //vectorsOfValuesH22.assign(rawDataReceiveH2[r].data() + index, rawDataReceiveH2[r].data() + index + vectorSize1);
            index += vectorSize2;

            vectorsOfValuesF3.assign(rawDataReceiveF[r].data() + index, rawDataReceiveF[r].data() + index + vectorSize3);
            if (multiPhase)
                vectorsOfValuesH13.assign(rawDataReceiveH1[r].data() + index, rawDataReceiveH1[r].data() + index + vectorSize1);
                //vectorsOfValuesH23.assign(rawDataReceiveH2[r].data() + index, rawDataReceiveH2[r].data() + index + vectorSize1);
            index += vectorSize3;

            SPtr<DistributionArray3D> mFdistributions(new D3Q27EsoTwist3DSplittedVector());
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesF1, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesF2, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                        vectorsOfValuesF3, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(dataSetParamStr1.nx1);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(dataSetParamStr1.nx2);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(dataSetParamStr1.nx3);

            SPtr<DistributionArray3D> mH1distributions(new D3Q27EsoTwist3DSplittedVector());
            if (multiPhase)
            {
                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH11, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH12, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                    vectorsOfValuesH13, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setNX1(dataSetParamStr1.nx1);
                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setNX2(dataSetParamStr1.nx2);
                dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH1distributions)->setNX3(dataSetParamStr1.nx3);
            }

            /*SPtr<DistributionArray3D> mH2distributions(new D3Q27EsoTwist3DSplittedVector());
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH21, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH22, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                    vectorsOfValuesH23, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX1(dataSetParamStr1.nx1);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX2(dataSetParamStr1.nx2);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX3(dataSetParamStr1.nx3);*/


            // find the nesessary block and fill it
            SPtr<Block3D> block = grid->getBlock(blockID);
            this->lbmKernel->setBlock(block);
            SPtr<LBMKernel> kernel = this->lbmKernel->clone();
            LBMReal collFactor = LBMSystem::calcCollisionFactor(this->nue, block->getLevel());
            LBMReal collFactorL = LBMSystem::calcCollisionFactor(this->nuL, block->getLevel());
            LBMReal collFactorG = LBMSystem::calcCollisionFactor(this->nuG, block->getLevel());
            kernel->setCollisionFactor(collFactor);
            kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
            kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
            kernel->setCollisionFactorMultiphase(collFactorL, collFactorG);
            kernel->setDensityRatio(this->densityRatio);
            SPtr<DataSet3D> dataSetPtr = SPtr<DataSet3D>(new DataSet3D());
            dataSetPtr->setFdistributions(mFdistributions);
            if (multiPhase)
                dataSetPtr->setHdistributions(mH1distributions);
//            dataSetPtr->setHdistributions(mH2distributions);
            kernel->setDataSet(dataSetPtr);
            block->setKernel(kernel);
        }
    }
    //if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    //-------------------------------------------------------------

    DSArraysPresence arrPresence;
    MPI_File file_handler1;
    std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler1);
    if (rc != MPI_SUCCESS)
        return; // throw UbException(UB_EXARGS, "couldn't open file " + filename1);

    MPI_File_read_at(file_handler1, (MPI_Offset)0, &arrPresence, 1, arrayPresenceType, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handler1);

    if (arrPresence.isAverageDensityArrayPresent)
        readArray(step, AverageDensity, std::string("/cpAverageDensityArray.bin"));
    // readAverageDensityArray(step);

    if (arrPresence.isAverageVelocityArrayPresent)
        readArray(step, AverageVelocity, std::string("/cpAverageVelocityArray.bin"));
    //   readAverageVelocityArray(step);

    if (arrPresence.isAverageFluktuationsArrayPresent)
        readArray(step, AverageFluktuations, std::string("/cpAverageFluktuationsArray.bin"));
    //   readAverageFluktuationsArray(step);

    if (arrPresence.isAverageTripleArrayPresent)
        readArray(step, AverageTriple, std::string("/cpAverageTripleArray.bin"));
    //  readAverageTripleArray(step);

    if (arrPresence.isShearStressValArrayPresent)
        readArray(step, ShearStressVal, std::string("/cpShearStressValArray.bin"));
    //   readShearStressValArray(step);

    if (arrPresence.isRelaxationFactorPresent)
        readArray(step, RelaxationFactor, std::string("/cpRelaxationFactor.bin"));
    //   readRelaxationFactor(step);
 
    if (arrPresence.isPhaseField1Present)
        readArray(step, PhaseField1, std::string("/cpPhaseField1.bin"));

    if (arrPresence.isPhaseField2Present)
        readArray(step, PhaseField2, std::string("/cpPhaseField2.bin"));

    delete[] rawDataReceiveF;
//    delete[] rawDataReceiveH1;
//    delete[] rawDataReceiveH2;
}

void MPIIOMigrationBECoProcessor::readArray(int step, Arrays arrType, std::string fname)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    dataSetParam dataSetParamStr;
    memset(&dataSetParamStr, 0, sizeof(dataSetParam));

    int blocksCountAll   = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
    int blocksPerProcess = blocksCountAll / size;     // how many blocks has each process

    size_t myBlocksCount;
    if (rank < (size - 1))
        myBlocksCount = blocksPerProcess;
    else
        myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

    int indexB = rank * blocksPerProcess;     // the first "my" block
    int indexE = indexB + int(myBlocksCount); // the latest "my" block

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

    size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
    MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], int(myBlocksCount), dataSetDoubleType, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handler);
    MPI_Type_free(&dataSetDoubleType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray start of exchange of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    std::vector<double> *rawDataReceive = new std::vector<double>[size];
    for (int r = 0; r < size; r++)
        rawDataReceive[r].resize(0);

    blocksExchange(MESSAGE_TAG + int(arrType), indexB, indexE, int(doubleCountInBlock), doubleValuesArray, rawDataReceive);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray end of exchange of data, rank = " << rank);
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    //----------------------------- restore data ---------------------------------
    int blockID;
    std::vector<double> vectorsOfValues;
    size_t index;
    size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

    for (int r = 0; r < size; r++) 
    {
        index = 0;
        for (int ii = 0; ii < int(rawDataReceive[r].size() / doubleCountInBlock); ii++) 
        {
            blockID = (int)(rawDataReceive[r][index]);
            SPtr<Block3D> block = grid->getBlock(blockID);
            index += 1;

            vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
            index += nextVectorSize;

            // fill arrays
            SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ___4DArray;
            SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> ___3DArray;

            switch (arrType) 
            {
                case AverageDensity:
                    ___4DArray = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(
                            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                    block->getKernel()->getDataSet()->setAverageDensity(___4DArray);
                    break;
                case AverageVelocity:
                    ___4DArray = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(
                            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                    block->getKernel()->getDataSet()->setAverageVelocity(___4DArray);
                    break;
                case AverageFluktuations:
                    ___4DArray = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(
                            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                    block->getKernel()->getDataSet()->setAverageFluctuations(___4DArray);
                    break;
                case AverageTriple:
                    ___4DArray = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(
                            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                    block->getKernel()->getDataSet()->setAverageTriplecorrelations(___4DArray);
                    break;
                case ShearStressVal:
                    ___4DArray = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(
                            vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));
                    block->getKernel()->getDataSet()->setShearStressValues(___4DArray);
                    break;
                case RelaxationFactor:
                    ___3DArray = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                        vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                    block->getKernel()->getDataSet()->setRelaxationFactor(___3DArray);
                    break;
                case PhaseField1:
                    ___3DArray = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                        vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                    block->getKernel()->getDataSet()->setPhaseField(___3DArray);
                    break;
                case PhaseField2:
                    ___3DArray = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                        vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));
                    block->getKernel()->getDataSet()->setPhaseField2(___3DArray);
                    break;
                default:
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationBECoProcessor::readArray : array type does not exist!"));
                    break;
            } 
        }
    }

    delete[] rawDataReceive;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}

void MPIIOMigrationBECoProcessor::readBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
    if (comm->isRoot())
        start = MPI_Wtime();

    int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
    size_t myBlocksCount;
    int blocksPerProcess = blocksCountAll / size; // how many blocks has each process

    if (rank < (size - 1))
        myBlocksCount = blocksPerProcess;
    else
        myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

    int indexB = rank * blocksPerProcess;     // the first "my" block
    int indexE = indexB + int(myBlocksCount); // the latest "my" block

    std::vector<int> bcindexmatrixVAll;

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    int sizeOfBIM;
    MPI_File_read_at(file_handler, (MPI_Offset)0, &sizeOfBIM, 1, MPI_INT, MPI_STATUS_IGNORE);
    bcindexmatrixVAll.resize(myBlocksCount * sizeOfBIM);

    MPI_Type_contiguous(sizeOfBIM, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    MPI_Offset read_offset = (MPI_Offset)(sizeof(int)) + (MPI_Offset)(indexB) * (MPI_Offset)(sizeOfBIM) * (MPI_Offset)(sizeof(int));
    MPI_File_read_at(file_handler, read_offset, &bcindexmatrixVAll[0], int(myBlocksCount), bcindexmatrixType, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds start of exchange of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    std::vector<int> *rawDataReceive = new std::vector<int>[size];
    std::vector<int> *rawDataSend    = new std::vector<int>[size];
    for (int r = 0; r < size; r++) 
    {
        rawDataReceive[r].resize(0);
        rawDataSend[r].resize(0);
        rawDataReceive[r].push_back(0);
        rawDataSend[r].push_back(0);
    }

    SPtr<Block3D> tempBlock;
    int tempRank;
    for (int ind = 0; ind < indexE - indexB; ind++) 
    {
        tempBlock = grid->getBlock(indexB + ind);
        tempRank  = tempBlock->getRank();

        if (tempRank == rank) // no need to send data, the process already has it
        {
            rawDataReceive[tempRank][0]++;
            rawDataReceive[tempRank].push_back(indexB + ind);
            rawDataReceive[tempRank].insert(rawDataReceive[tempRank].end(), bcindexmatrixVAll.begin() + ind * sizeOfBIM, bcindexmatrixVAll.begin() + ind * sizeOfBIM + sizeOfBIM);
        } else // we must send data to other processes
        {
            rawDataSend[tempRank][0]++;
            rawDataSend[tempRank].push_back(indexB + ind);
            rawDataSend[tempRank].insert(rawDataSend[tempRank].end(), bcindexmatrixVAll.begin() + ind * sizeOfBIM, bcindexmatrixVAll.begin() + ind * sizeOfBIM + sizeOfBIM);
        }
    }

    MPI_Request *requests = new MPI_Request[size * 2]; // send + receive
    int requestCount      = 0;
    MPI_Status status;
    int quant;
    int intBlockCount;
    int rds;

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            rds = int(rawDataSend[r].size());
            intBlockCount = (int)(rds / SEND_BLOCK_SIZE);
            if (intBlockCount * SEND_BLOCK_SIZE < rds)
                intBlockCount += 1;

            for (int i = rds; i < intBlockCount * SEND_BLOCK_SIZE; i++)
                rawDataSend[r].push_back(0);

            MPI_Isend(&rawDataSend[r][0], intBlockCount, sendBlockIntType, r, MESSAGE_TAG + 7, MPI_COMM_WORLD, &requests[requestCount]);
            // MPI_Isend(&rawDataSend[r][0], rawDataSend[r].size(), MPI_INT, r, MESSAGE_TAG + 7, MPI_COMM_WORLD,
            // &requests[requestCount]);
            requestCount++;
        }
    }

    for (int r = 0; r < size; r++) 
    {
        if (r != rank) 
        {
            MPI_Probe(r, MESSAGE_TAG + 7, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, sendBlockIntType, &quant);
            rawDataReceive[r].resize(quant * SEND_BLOCK_SIZE);
            MPI_Irecv(&rawDataReceive[r][0], quant, sendBlockIntType, r, MESSAGE_TAG + 7, MPI_COMM_WORLD, &requests[requestCount]);
            requestCount++;
        }
    }

    MPI_Waitall(requestCount, &requests[0], MPI_STATUSES_IGNORE);

    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds end of exchange of data, rank = " << rank);
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC2.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    MPI_File_read_at(file_handler, (MPI_Offset)0, &boundCondParamStr, 4, MPI_INT, MPI_STATUS_IGNORE);

    int blockID;
    size_t index;
    MPI_Offset read_offset1, read_offset2;

    BCAddMigration bcAddArray;
    BoundaryCondition *nullBouCond = new BoundaryCondition();
    memset(nullBouCond, 0, sizeof(BoundaryCondition));
    BoundaryCondition *bcArray = nullptr;
    std::vector<SPtr<BoundaryConditions>> bcVector;
    std::vector<int> indexContainerV;
    std::vector<int> bcindexmatrixV;

    for (int r = 0; r < size; r++) 
    {
        index = 1;

        for (int ii = 0; ii < rawDataReceive[r][0]; ii++) 
        {
            blockID = (int)(rawDataReceive[r][index]);
            index += 1;

            bcindexmatrixV.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + sizeOfBIM);
            index += sizeOfBIM;

            read_offset1 = (MPI_Offset)(sizeof(boundCondParam)) + (MPI_Offset)(blockID) * (MPI_Offset)(sizeof(size_t));

            MPI_File_read_at(file_handler, read_offset1, &read_offset2, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
            MPI_File_read_at(file_handler, read_offset2, &bcAddArray, 3, MPI_INT, MPI_STATUS_IGNORE);

            bcArray = new BoundaryCondition[bcAddArray.boundCond_count];
            indexContainerV.resize(bcAddArray.indexContainer_count);

            if (bcAddArray.boundCond_count > 0)
                MPI_File_read_at(file_handler, read_offset2 + (MPI_Offset)(sizeof(BCAddMigration)), &bcArray[0],
                                 bcAddArray.boundCond_count, boundCondType, MPI_STATUS_IGNORE);

            if (bcAddArray.indexContainer_count > 0)
                MPI_File_read_at(file_handler, read_offset2 + (MPI_Offset)(sizeof(BCAddMigration)) +
                                     (MPI_Offset)(bcAddArray.boundCond_count) * (MPI_Offset)(sizeof(BoundaryCondition)),
                                 &indexContainerV[0], bcAddArray.indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);

            bcVector.resize(0);

            for (int ibc = 0; ibc < bcAddArray.boundCond_count; ibc++) 
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
                    bc->bcLodiDensity          = bcArray[ibc].bcLodiDensity;
                    bc->bcLodiVelocityX1       = bcArray[ibc].bcLodiVelocityX1;
                    bc->bcLodiVelocityX2       = bcArray[ibc].bcLodiVelocityX2;
                    bc->bcLodiVelocityX3       = bcArray[ibc].bcLodiVelocityX3;
                    bc->bcLodiLentgh           = bcArray[ibc].bcLodiLentgh;

                    bc->nx1 = bcArray[ibc].nx1;
                    bc->nx2 = bcArray[ibc].nx2;
                    bc->nx3 = bcArray[ibc].nx3;
                    for (int iq = 0; iq < 26; iq++)
                        bc->setQ(bcArray[ibc].q[iq], iq);
                    bc->setBcAlgorithmType(bcArray[ibc].algorithmType);
                }

                bcVector.push_back(bc);
            }

            CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);
            SPtr<Block3D> block1 = grid->getBlock(blockID);

            SPtr<BCProcessor> bcProc = bcProcessor->clone(block1->getKernel());
            SPtr<BCArray3D> bcArr(new BCArray3D());
            bcArr->bcindexmatrix  = bcim;
            bcArr->bcvector       = bcVector;
            bcArr->indexContainer = indexContainerV;
            bcProc->setBCArray(bcArr);

            block1->getKernel()->setBCProcessor(bcProc);
        }
    }

    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    MPI_File_close(&file_handler);

    delete nullBouCond;
    delete[] bcArray;
    delete[] rawDataReceive;
    delete[] rawDataSend;
    delete[] requests;

    if (comm->isRoot()) 
    {
        finish = MPI_Wtime();
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds end of restore of data, rank = " << rank);
        UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setLBMKernel(SPtr<LBMKernel> kernel) { this->lbmKernel = kernel; }
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setBCProcessor(SPtr<BCProcessor> bcProcessor) { this->bcProcessor = bcProcessor; }
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setNu(double nu) { this->nue = nu; }

void MPIIOMigrationBECoProcessor::setNuLG(double cfL, double cfG) { this->nuL = cfL;  this->nuG = cfG; }

void MPIIOMigrationBECoProcessor::setDensityRatio(double dr) { this->densityRatio = dr; }

