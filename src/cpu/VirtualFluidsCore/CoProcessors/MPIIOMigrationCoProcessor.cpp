#include "MPIIOMigrationCoProcessor.h"
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
#include "RenumberBlockVisitor.h"
#include "UbFileInputASCII.h"
#include "UbFileOutputASCII.h"
#include "UbScheduler.h"
#include "WbWriter.h"
#include <MemoryUtil.h>
#include <UbSystem.h>

using namespace MPIIODataStructures;

MPIIOMigrationCoProcessor::MPIIOMigrationCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, SPtr<Communicator> comm)
    : MPIIOCoProcessor(grid, s, path, comm)
{
    memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));

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
MPIIOMigrationCoProcessor::~MPIIOMigrationCoProcessor()
{
    MPI_Type_free(&dataSetType);
    MPI_Type_free(&dataSetSmallType);
    MPI_Type_free(&boundCondParamType);
    MPI_Type_free(&boundCondTypeAdd);
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::process(double step)
{
    if (scheduler->isDue(step)) 
    {
        if (comm->isRoot())
            UBLOG(logINFO, "MPIIOMigrationCoProcessor save step: " << step);
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

void MPIIOMigrationCoProcessor::clearAllFiles(int step)
{
    MPI_File file_handler;
    MPI_Info info       = MPI_INFO_NULL;
    MPI_Offset new_size = 0;

    MPIIOCoProcessor::clearAllFiles(step);

    UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

    std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
    int rc10 = MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc10 != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename10);
    MPI_File_set_size(file_handler, new_size);
    MPI_File_close(&file_handler);
}

void MPIIOMigrationCoProcessor::writeBlocks(int step)
{
    grid->renumberBlockIDs();
    MPIIOCoProcessor::writeBlocks(step);
}

void MPIIOMigrationCoProcessor::writeDataSet(int step)
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
    std::vector<double> doubleValuesArrayF; // double-values (arrays of f's) in all blocks  Fdistribution
    std::vector<double> doubleValuesArrayH1; // double-values (arrays of f's) in all blocks  H1distribution
    std::vector<double> doubleValuesArrayH2; // double-values (arrays of f's) in all blocks  H2distribution

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeDataSet start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool multiPhase1 = false;
    bool multiPhase2 = false;
    DSArraysPresence arrPresence;
    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<D3Q27EsoTwist3DSplittedVector> D3Q27EsoTwist3DSplittedVectorPtrF = 0, D3Q27EsoTwist3DSplittedVectorPtrH1 = 0, D3Q27EsoTwist3DSplittedVectorPtrH2 = 0;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF = 0, localDistributionsH1 = 0, localDistributionsH2 = 0;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF = 0, nonLocalDistributionsH1 = 0, nonLocalDistributionsH2 = 0;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributionsF = 0, zeroDistributionsH1 = 0, zeroDistributionsH2 = 0;

    SPtr<LBMKernel> kernel;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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

            D3Q27EsoTwist3DSplittedVectorPtrF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
            localDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getLocalDistributions();
            nonLocalDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getNonLocalDistributions();
            zeroDistributionsF = D3Q27EsoTwist3DSplittedVectorPtrF->getZeroDistributions();

            D3Q27EsoTwist3DSplittedVectorPtrH1 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getHdistributions());
            if (D3Q27EsoTwist3DSplittedVectorPtrH1 != 0)
            {
                multiPhase1 = true;
                localDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getLocalDistributions();
                nonLocalDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getNonLocalDistributions();
                zeroDistributionsH1 = D3Q27EsoTwist3DSplittedVectorPtrH1->getZeroDistributions();
            }

            D3Q27EsoTwist3DSplittedVectorPtrH2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getH2distributions());
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeDataSet start MPI IO rank = " << rank);
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

    MPI_File_write_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
    MPI_File_write_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
    
    MPI_Offset write_offset;
    size_t sizeofOneDataSet = sizeof(DataSetMigration) + doubleCountInBlock * sizeof(double);

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

        sizeofOneDataSet = doubleCountInBlock * sizeof(double);

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

        sizeofOneDataSet = doubleCountInBlock * sizeof(double);

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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeDataSet time: " << finish - start << " s");
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

    if (arrPresence.isPhaseField1Present)
        write3DArray(step, PhaseField1, std::string("/cpPhaseField1.bin"));

    if (arrPresence.isPhaseField2Present)
        write3DArray(step, PhaseField2, std::string("/cpPhaseField2.bin"));

}

void MPIIOMigrationCoProcessor::write4DArray(int step, Arrays arrayType, std::string fname)
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
    std::vector<double> doubleValuesArray; // double-values of the AverageDensityArray in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<CbArray4D<LBMReal, IndexerX4X3X2X1>> ___Array;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationCoProcessor::write4DArray : 4D array type does not exist!"));
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::write4DArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
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
    size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::write4DArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::write3DArray(int step, Arrays arrayType, std::string fname)
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
    std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
    dataSetParam dataSetParamStr;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::write3DArray start collect data rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    bool firstBlock           = true;
    size_t doubleCountInBlock = 0;
    int ic                    = 0;
    SPtr<CbArray3D<LBMReal, IndexerX3X2X1>> ___Array;

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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
                default:
                    UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationCoProcessor::write3DArray : 3D array type does not exist!"));
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::write3DArray start MPI IO rank = " << rank);
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

    MPI_File file_handler;
    std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + fname;
    int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    // each process writes common parameters of a dataSet
    MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

    size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::write3DArray time: " << finish - start << " s");
    }

    delete[] dataSetSmallArray;
}

/*
void MPIIOMigrationCoProcessor::writeAverageDensityArray(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values of the AverageDensityArray in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray =
block->getKernel()->getDataSet()->getAverageDensity();

         if (firstBlock) // when first (any) valid block...
         {
            //if (averageDensityArray)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(averageDensityArray->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(averageDensityArray->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(averageDensityArray->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(averageDensityArray->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (averageDensityArray && (dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2]
> 0) && (dataSetParamStr.nx[3] > 0)) doubleValuesArray.insert(doubleValuesArray.end(),
averageDensityArray->getDataVector().begin(), averageDensityArray->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::writeAverageVelocityArray(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageVelocityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr =
block->getKernel()->getDataSet()->getAverageVelocity();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageVelocityArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageVelocityArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) &&
(dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0)) doubleValuesArray.insert(doubleValuesArray.end(),
AverageVelocityArray3DPtr->getDataVector().begin(), AverageVelocityArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageVelocityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageVelocityArray time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::writeAverageFluktuationsArray(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageFluktuationsArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr =
block->getKernel()->getDataSet()->getAverageFluctuations();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageFluctArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageFluctArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) &&
(dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0)) doubleValuesArray.insert(doubleValuesArray.end(),
AverageFluctArray3DPtr->getDataVector().begin(), AverageFluctArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageFluktuationsArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;
#ifdef HLRN_LUSTRE
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageFluktuationsArray time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::writeAverageTripleArray(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageTripleArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr =
block->getKernel()->getDataSet()->getAverageTriplecorrelations();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageTripleArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageTripleArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) &&
(dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0)) doubleValuesArray.insert(doubleValuesArray.end(),
AverageTripleArray3DPtr->getDataVector().begin(), AverageTripleArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageTripleArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageTripleArray time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::writeShearStressValArray(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeShearStressValArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr =
block->getKernel()->getDataSet()->getShearStressValues();

         if (firstBlock) // when first (any) valid block...
         {
            //if (ShearStressValArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (ShearStressValArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) &&
(dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0)) doubleValuesArray.insert(doubleValuesArray.end(),
ShearStressValArray3DPtr->getDataVector().begin(), ShearStressValArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeShearStressValArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeShearStressValArray time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::writeRelaxationFactor(int step)
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeRelaxationFactor start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   bool firstBlock = true;
   size_t doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while
regenerating the grid

         SPtr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr =
block->getKernel()->getDataSet()->getRelaxationFactor();

         if (firstBlock) // when first (any) valid block...
         {
            //if (relaxationFactor3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(relaxationFactor3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(relaxationFactor3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(relaxationFactor3DPtr->getNX3());
            dataSetParamStr.nx[3] = 1;
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (relaxationFactor3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) &&
(dataSetParamStr.nx[2]>0)) doubleValuesArray.insert(doubleValuesArray.end(),
relaxationFactor3DPtr->getDataVector().begin(), relaxationFactor3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeRelaxationFactor start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double);

   for (size_t nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[nb
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeRelaxationFactor time: " << finish - start << " s");
   }

   delete[] dataSetSmallArray;
}
*/
void MPIIOMigrationCoProcessor::writeBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBoundaryConds start collect data rank = " << rank);
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
    std::vector<int> *bcindexmatrixVector    = new std::vector<int>[blocksCount];
    std::vector<int> *indexContainerVector   = new std::vector<int>[blocksCount];

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
            bcindexmatrixVector[ic].insert(bcindexmatrixVector[ic].begin(), bcArr->bcindexmatrix.getDataVector().begin(), bcArr->bcindexmatrix.getDataVector().end());
            bytesCount[ic] += boundCondParamStr.bcindexmatrixCount * sizeof(int);

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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: "<< Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBoundaryConds time: " << finish - start << " s");
    }

    delete[] bcAddArray;
    delete[] bytesCount;
    delete[] bcVector;
    delete[] bcindexmatrixVector;
    delete[] indexContainerVector;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIOMigrationCoProcessor::restart(int step)
{
    if (comm->isRoot())
        UBLOG(logINFO, "MPIIOMigrationCoProcessor restart step: " << step);
    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - start");

    readBlocks(step);

    SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
    grid->accept(metisVisitor);

    readDataSet(step);
    readBoundaryConds(step);

    grid->setTimeStep(step);

    if (comm->isRoot())
        UBLOG(logINFO, "Load check point - end");
}

void MPIIOMigrationCoProcessor::readBlocks(int step) { MPIIOCoProcessor::readBlocks(step); }

void MPIIOMigrationCoProcessor::readDataSet(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
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
    std::vector<double> doubleValuesArrayF(size_t(blocksCount * doubleCountInBlock)); // double-values in all blocks  Fdistributions
    std::vector<double> doubleValuesArrayH1; // double-values in all blocks  H1distributions
    std::vector<double> doubleValuesArrayH2; // double-values in all blocks  H2distributions

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    size_t ic = 0;
    MPI_Offset read_offset;
    size_t sizeofOneDataSet = size_t(sizeof(DataSetMigration) + doubleCountInBlock * sizeof(double));

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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
    MPI_Offset fsize;
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH1.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);
    int fs = MPI_File_get_size(file_handler, &fsize);
    if (fsize > 0)
    {
        multiPhase1 = true;
        doubleValuesArrayH1.resize(blocksCount * doubleCountInBlock);

        sizeofOneDataSet = size_t(doubleCountInBlock * sizeof(double));

        for (int level = minInitLevel; level <= maxInitLevel; level++)
        {
            for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
            {
                read_offset = (MPI_Offset)(block->getGlobalID() * sizeofOneDataSet);
                MPI_File_read_at(file_handler, read_offset, &doubleValuesArrayH1[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
                ic++;
            }
        }

    }
    MPI_File_close(&file_handler);

    //----------------------------------------- H2 ----------------------------------------------------
    filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSetH2.bin";
    rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
    if (rc != MPI_SUCCESS)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    fs = MPI_File_get_size(file_handler, &fsize);
    if (fsize > 0)
    {
        multiPhase2 = true;
        doubleValuesArrayH2.resize(blocksCount * doubleCountInBlock);

        sizeofOneDataSet = size_t(doubleCountInBlock * sizeof(double));

        for (int level = minInitLevel; level <= maxInitLevel; level++)
        {
            for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    size_t index = 0;
    std::vector<double> vectorsOfValuesF1, vectorsOfValuesF2, vectorsOfValuesF3;
    std::vector<double> vectorsOfValuesH11, vectorsOfValuesH12, vectorsOfValuesH13;
    std::vector<double> vectorsOfValuesH21, vectorsOfValuesH22, vectorsOfValuesH23;

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
       if (multiPhase1)
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

        SPtr<DistributionArray3D> mH2distributions(new D3Q27EsoTwist3DSplittedVector());
        if (multiPhase2)
        {
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH21, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(
                    new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValuesH22, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(
                    vectorsOfValuesH23, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX1(dataSetParamStr1.nx1);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX2(dataSetParamStr1.nx2);
            dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mH2distributions)->setNX3(dataSetParamStr1.nx3);
        }

        // find the nesessary block and fill it
        SPtr<Block3D> block = grid->getBlock(dataSetArray[n].globalID);
        this->lbmKernel->setBlock(block);
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet end of restore of data, rank = " << rank);
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

}

void MPIIOMigrationCoProcessor::readArray(int step, Arrays arrType, std::string fname)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readArray start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
    double start, finish;
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
    std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

    // define MPI_types depending on the block-specific information
    MPI_Type_contiguous(int(doubleCountInBlock), MPI_DOUBLE, &dataSetDoubleType);
    MPI_Type_commit(&dataSetDoubleType);

    size_t ic = 0;
    MPI_Offset read_offset;
    size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readArray readArray: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readArray start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    //----------------------------- restore data ---------------------------------
    size_t index = 0;
    size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
    std::vector<double> vectorsOfValues;
    for (std::size_t n = 0; n < blocksCount; n++) 
    {
        SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);

        vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
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
                UB_THROW(UbException(UB_EXARGS, "MPIIOMigrationCoProcessor::readArray : array type does not exist!"));
                break;
        }
    }

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readArray end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    delete[] dataSetSmallArray;
}

/*void MPIIOMigrationCoProcessor::readAverageDensityArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill mAverageDensity arrays
      SPtr<AverageValuesArray3D> mAverageDensity;
      //if
((dataSetParamStr.nx[0]==0)&&(dataSetParamStr.nx[1]==0)&&(dataSetParamStr.nx[2]==0)&&(dataSetParamStr.nx[3]==0))
      //   mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,
IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2],
dataSetParamStr.nx[3]));

      //std::cout << "rank=" << rank << ", dataSetArray[n].globalID=" << dataSetSmallArray[n].globalID << std::endl;
      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageDensity(mAverageDensity);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::readAverageVelocityArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill mAverageVelocity array
      SPtr<AverageValuesArray3D> mAverageVelocity;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) &&
(dataSetParamStr.nx[3] == 0))
      //   mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,
IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2],
dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageVelocity(mAverageVelocity);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::readAverageFluktuationsArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray start of restore of data, rank = " <<
rank); UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill AverageFluktuations array
      SPtr<AverageValuesArray3D> mAverageFluktuations;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) &&
(dataSetParamStr.nx[3] == 0))
      //   mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,
IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2],
dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageFluctuations(mAverageFluktuations);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::readAverageTripleArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill AverageTriplecorrelations array
      SPtr<AverageValuesArray3D> mAverageTriplecorrelations;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) &&
(dataSetParamStr.nx[3] == 0))
      //   mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,
IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2],
dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageTriplecorrelations(mAverageTriplecorrelations);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::readShearStressValArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill ShearStressValuesArray array
      SPtr<ShearStressValuesArray3D> mShearStressValues;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) &&
(dataSetParamStr.nx[3] == 0))
      //   mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal,
IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2],
dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setShearStressValues(mShearStressValues);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

void MPIIOMigrationCoProcessor::readRelaxationFactor(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
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

   DataSetSmallMigration* dataSetSmallArray = new DataSetSmallMigration[blocksCount];
   size_t doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all
blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   size_t ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmallMigration) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmallMigration)), &doubleValuesArray[ic
* doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE); ic++;
      }
   }

   MPI_File_close(&file_handler);
   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] *
dataSetParamStr.nx[3]; std::vector<double> vectorsOfValues; for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill RelaxationFactor array
      SPtr<RelaxationFactorArray3D> mRelaxationFactor;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0))
      //   mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr();
      //else
      mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal,
IndexerX3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setRelaxationFactor(mRelaxationFactor);
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() /
1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}
*/

void MPIIOMigrationCoProcessor::readBoundaryConds(int step)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBoundaryConds start MPI IO rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    double start, finish;
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
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBoundaryConds time: " << finish - start << " s");
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBoundaryConds start of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }

    MPI_File_read_at(file_handler, (MPI_Offset)0, &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);
    MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
    MPI_Type_commit(&bcindexmatrixType);

    int ic = 0;
    MPI_Offset read_offset1, read_offset2;
    for (int level = minInitLevel; level <= maxInitLevel; level++) 
    {
        for (SPtr<Block3D> block : blocksVector[level]) //	blocks of the current level
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

            for (int b1 = 0; b1 < boundCondParamStr.bcindexmatrixCount; b1++)
                bcindexmatrixV.push_back(intArray1[b1]);

            for (int b2 = 0; b2 < bcAddArray[ic].indexContainer_count; b2++)
                indexContainerV.push_back(intArray2[b2]);

            CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);
            SPtr<Block3D> block1 = grid->getBlock(bcAddArray[ic].globalID);

            SPtr<BCProcessor> bcProc = bcProcessor->clone(block1->getKernel());
            SPtr<BCArray3D> bcArr(new BCArray3D());
            bcArr->bcindexmatrix  = bcim;
            bcArr->bcvector       = bcVector;
            bcArr->indexContainer = indexContainerV;
            bcProc->setBCArray(bcArr);

            block1->getKernel()->setBCProcessor(bcProc);

            delete bcArray;
            delete intArray1;

            ic++;
        }
    }

    MPI_File_close(&file_handler);
    MPI_Type_free(&bcindexmatrixType);

    delete nullBouCond;

    if (comm->isRoot()) 
    {
        UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBoundaryConds end of restore of data, rank = " << rank);
        UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
    }
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::setLBMKernel(SPtr<LBMKernel> kernel) { this->lbmKernel = kernel; }
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::setBCProcessor(SPtr<BCProcessor> bcProcessor) { this->bcProcessor = bcProcessor; }
