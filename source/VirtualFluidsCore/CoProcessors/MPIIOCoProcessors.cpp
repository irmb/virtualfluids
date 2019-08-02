#include "MPIIOCoProcessors.h"
#include "Block3D.h"
#include "CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "Communicator.h"
#include "UbScheduler.h"
#include "MPIIODataStructures.h"
#include "UbLogger.h"
#include "MemoryUtil.h"

MPIIOCoProcessor::MPIIOCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<Communicator> comm) :
   CoProcessor(grid, s),
   path(path),
   comm(comm)
{
   UbSystem::makeDirectory(path + "/mpi_io_cp");

   //-------------------------   define MPI types  ---------------------------------

   MPI_Datatype typesGP[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksGP[3] = { 34, 6, 5 };
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
   int blocksBlock[2] = { 13, 1 };
   MPI_Aint offsetsBlock[2], lbBlock, extentBlock;

   offsetsBlock[0] = 0;
   MPI_Type_get_extent(MPI_INT, &lbBlock, &extentBlock);
   offsetsBlock[1] = blocksBlock[0] * extentBlock;

   MPI_Type_create_struct(2, blocksBlock, offsetsBlock, typesBlock, &block3dType);
   MPI_Type_commit(&block3dType);
}

MPIIOCoProcessor::~MPIIOCoProcessor()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&block3dType);
}

void MPIIOCoProcessor::readBlocks(int step)
{
   using namespace MPIIODataStructures;

   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
   MPI_File_read_at(file_handler, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   Block3d* block3dArray = new Block3d[blocksCount];

   GridParam* gridParameters = new GridParam;

   // calculate the read offset
   //MPI_Offset read_offset = (MPI_Offset)(size * sizeof(int));
   MPI_Offset read_offset = (MPI_Offset)(sizeof(int));

   // read parameters of the grid
   MPI_File_read_at(file_handler, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   // read all the blocks
   if (comm->isRoot())
      MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);

   MPI_Bcast(block3dArray, blocksCount, block3dType, comm->getRoot(), MPI_COMM_WORLD);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBlocks time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBlocks start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   // clear the grid
   grid->deleteBlocks();

   // restore the grid
   SPtr<CoordinateTransformation3D> trafo(new CoordinateTransformation3D());
   trafo->Tx1 = gridParameters->trafoParams[0];
   trafo->Tx2 = gridParameters->trafoParams[1];
   trafo->Tx3 = gridParameters->trafoParams[2];
   trafo->Sx1 = gridParameters->trafoParams[3];
   trafo->Sx2 = gridParameters->trafoParams[4];
   trafo->Sx3 = gridParameters->trafoParams[5];
   trafo->alpha = gridParameters->trafoParams[6];
   trafo->beta = gridParameters->trafoParams[7];
   trafo->gamma = gridParameters->trafoParams[8];

   trafo->toX1factorX1 = gridParameters->trafoParams[9];
   trafo->toX1factorX2 = gridParameters->trafoParams[10];
   trafo->toX1factorX3 = gridParameters->trafoParams[11];
   trafo->toX1delta = gridParameters->trafoParams[12];
   trafo->toX2factorX1 = gridParameters->trafoParams[13];
   trafo->toX2factorX2 = gridParameters->trafoParams[14];
   trafo->toX2factorX3 = gridParameters->trafoParams[15];
   trafo->toX2delta = gridParameters->trafoParams[16];
   trafo->toX3factorX1 = gridParameters->trafoParams[17];
   trafo->toX3factorX2 = gridParameters->trafoParams[18];
   trafo->toX3factorX3 = gridParameters->trafoParams[19];
   trafo->toX3delta = gridParameters->trafoParams[20];

   trafo->fromX1factorX1 = gridParameters->trafoParams[21];
   trafo->fromX1factorX2 = gridParameters->trafoParams[22];
   trafo->fromX1factorX3 = gridParameters->trafoParams[23];
   trafo->fromX1delta = gridParameters->trafoParams[24];
   trafo->fromX2factorX1 = gridParameters->trafoParams[25];
   trafo->fromX2factorX2 = gridParameters->trafoParams[26];
   trafo->fromX2factorX3 = gridParameters->trafoParams[27];
   trafo->fromX2delta = gridParameters->trafoParams[28];
   trafo->fromX3factorX1 = gridParameters->trafoParams[29];
   trafo->fromX3factorX2 = gridParameters->trafoParams[30];
   trafo->fromX3factorX3 = gridParameters->trafoParams[31];
   trafo->fromX3delta = gridParameters->trafoParams[32];

   trafo->active = gridParameters->active;
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
   for (int n = 0; n < blocksCount; n++)
   {
      SPtr<Block3D> block(new Block3D(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, block3dArray[n].level));
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

}

void MPIIOCoProcessor::clearAllFiles(int step)
{
   MPI_File file_handler;
   MPI_Info info = MPI_INFO_NULL;
   MPI_Offset new_size = 0;

   UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));

   std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rc1 = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler);
   if (rc1 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename1);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename2 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rc2 = MPI_File_open(MPI_COMM_WORLD, filename2.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc2 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename2);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename3 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
   int rc3 = MPI_File_open(MPI_COMM_WORLD, filename3.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc3 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename3);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename4 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   //MPI_File_delete(filename4.c_str(), info);
   int rc4 = MPI_File_open(MPI_COMM_WORLD, filename4.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc4 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename4);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename5 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   //MPI_File_delete(filename5.c_str(), info);
   int rc5 = MPI_File_open(MPI_COMM_WORLD, filename5.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc5 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename5);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename6 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   //MPI_File_delete(filename6.c_str(), info);
   int rc6 = MPI_File_open(MPI_COMM_WORLD, filename6.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc6 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename6);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename7 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   //MPI_File_delete(filename7.c_str(), info);
   int rc7 = MPI_File_open(MPI_COMM_WORLD, filename7.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc7 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename7);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename8 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   //MPI_File_delete(filename8.c_str(), info);
   int rc8 = MPI_File_open(MPI_COMM_WORLD, filename8.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc8 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename8);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename9 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   //MPI_File_delete(filename9.c_str(), info);
   int rc9 = MPI_File_open(MPI_COMM_WORLD, filename9.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc9 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename9);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
   int rc10 = MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc10 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename10);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);

   std::string filename11 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC2.bin";
   int rc11 = MPI_File_open(MPI_COMM_WORLD, filename11.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc11 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename11);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);
}
