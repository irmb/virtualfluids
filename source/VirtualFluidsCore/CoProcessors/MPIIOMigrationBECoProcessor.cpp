#include "MPIIOMigrationBECoProcessor.h"
#include "D3Q27System.h"
#include "LBMKernel.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <UbSystem.h>
#include <MemoryUtil.h>
#include "BoundaryConditions.h"
#include "Block3D.h"
#include "CoordinateTransformation3D.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "BCArray3D.h"
#include "Communicator.h"
#include "WbWriter.h"
#include "UbScheduler.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "MetisPartitioningGridVisitor.h"
#include "PointerDefinitions.h"
#include "RenumberGridVisitor.h"
#include "UbFileOutputASCII.h"
#include "UbFileInputASCII.h"

using namespace MPIIODataStructures;

#define MESSAGE_TAG 80
#define SEND_BLOCK_SIZE 100000

MPIIOMigrationBECoProcessor::MPIIOMigrationBECoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
   const std::string& path,
   SPtr<Communicator> comm) :
   CoProcessor(grid, s),
   path(path),
   comm(comm), 
   nue(-999.999)
{
   UbSystem::makeDirectory(path + "/mpi_io_cp");

   memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));

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

   //-----------------------------------------------------------------------

   MPI_Datatype typesBC[3] = { MPI_LONG_LONG_INT, MPI_FLOAT, MPI_CHAR };
   int blocksBC[3] = { 5, 38, 1 };
   MPI_Aint offsetsBC[3], lbBC, extentBC;

   offsetsBC[0] = 0;
   MPI_Type_get_extent(MPI_LONG_LONG_INT, &lbBC, &extentBC);
   offsetsBC[1] = blocksBC[0] * extentBC;

   MPI_Type_get_extent(MPI_FLOAT, &lbBC, &extentBC);
   offsetsBC[2] = offsetsBC[1] + blocksBC[1] * extentBC;

   MPI_Type_create_struct(3, blocksBC, offsetsBC, typesBC, &boundCondType);
   MPI_Type_commit(&boundCondType);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(7, MPI_INT, &dataSetParamType);
   MPI_Type_commit(&dataSetParamType);

   //---------------------------------------

   MPI_Type_contiguous(6, MPI_CHAR, &arrayPresenceType);
   MPI_Type_commit(&arrayPresenceType);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(SEND_BLOCK_SIZE, MPI_DOUBLE, &sendBlockDoubleType);
   MPI_Type_commit(&sendBlockDoubleType);

   MPI_Type_contiguous(SEND_BLOCK_SIZE, MPI_INT, &sendBlockIntType);
   MPI_Type_commit(&sendBlockIntType);

}

//////////////////////////////////////////////////////////////////////////
MPIIOMigrationBECoProcessor::~MPIIOMigrationBECoProcessor()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&block3dType);
   MPI_Type_free(&boundCondType);
   MPI_Type_free(&dataSetParamType);
   MPI_Type_free(&sendBlockDoubleType);
   MPI_Type_free(&sendBlockIntType);
   MPI_Type_free(&arrayPresenceType);
}

void MPIIOMigrationBECoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      if (comm->isRoot()) UBLOG(logINFO, "MPIIOMigrationBECoProcessor save step: " << step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - start");
      clearAllFiles((int)step);

      writeBlocks((int)step);
      writeDataSet((int)step);
      writeBoundaryConds((int)step);

      writeCpTimeStep((int)step);

      if (comm->isRoot()) UBLOG(logINFO, "Save check point - end");
   }
}

void MPIIOMigrationBECoProcessor::clearAllFiles(int step)
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

void MPIIOMigrationBECoProcessor::writeCpTimeStep(int step)
{
   if (comm->isRoot())
   {
      UbFileOutputASCII f(path + "/mpi_io_cp/cp.txt");
      f.writeInteger(step);
   }
}
//////////////////////////////////////////////////////////////////////////
int MPIIOMigrationBECoProcessor::readCpTimeStep()
{
   UbFileInputASCII f(path + "/mpi_io_cp/cp.txt");
   int step = f.readInteger();
   return step;
}

void MPIIOMigrationBECoProcessor::writeBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

	grid->deleteBlockIDs();
	RenumberGridVisitor renumber(comm);
	grid->accept(renumber);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBlocks start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blocksCount = 0; // quantity of all the blocks in the grid, max 2147483648 blocks!
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   std::vector<SPtr<Block3D>> blocksVector[25]; // max 25 levels
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      grid->getBlocks(level, blocksVector[level]);
      blocksCount += static_cast<int>(blocksVector[level].size());
   }

   GridParam* gridParameters = new GridParam;
   gridParameters->trafoParams[0] = grid->getCoordinateTransformator()->Tx1;
   gridParameters->trafoParams[1] = grid->getCoordinateTransformator()->Tx2;
   gridParameters->trafoParams[2] = grid->getCoordinateTransformator()->Tx3;
   gridParameters->trafoParams[3] = grid->getCoordinateTransformator()->Sx1;
   gridParameters->trafoParams[4] = grid->getCoordinateTransformator()->Sx2;
   gridParameters->trafoParams[5] = grid->getCoordinateTransformator()->Sx3;
   gridParameters->trafoParams[6] = grid->getCoordinateTransformator()->alpha;
   gridParameters->trafoParams[7] = grid->getCoordinateTransformator()->beta;
   gridParameters->trafoParams[8] = grid->getCoordinateTransformator()->gamma;

   gridParameters->trafoParams[9] = grid->getCoordinateTransformator()->toX1factorX1;
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

   gridParameters->active = grid->getCoordinateTransformator()->active;
   gridParameters->transformation = grid->getCoordinateTransformator()->transformation;
   
   gridParameters->deltaX = grid->getDeltaX(minInitLevel);
   UbTupleInt3 blocknx = grid->getBlockNX();
   gridParameters->blockNx1 = val<1>(blocknx);
   gridParameters->blockNx2 = val<2>(blocknx);
   gridParameters->blockNx3 = val<3>(blocknx);
   gridParameters->nx1 = grid->getNX1();
   gridParameters->nx2 = grid->getNX2();
   gridParameters->nx3 = grid->getNX3();
   gridParameters->periodicX1 = grid->isPeriodicX1();
   gridParameters->periodicX2 = grid->isPeriodicX2();
   gridParameters->periodicX3 = grid->isPeriodicX3();

   //----------------------------------------------------------------------

   Block3d* block3dArray = new Block3d[blocksCount];
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	all the blocks of the current level
      {
         // save data describing the block
         block3dArray[ic].x1 = block->getX1();
         block3dArray[ic].x2 = block->getX2();
         block3dArray[ic].x3 = block->getX3();
         block3dArray[ic].bundle = block->getBundle();
         block3dArray[ic].rank = block->getRank();
         block3dArray[ic].lrank = block->getLocalRank();
         block3dArray[ic].part = block->getPart();
         block3dArray[ic].globalID = block->getGlobalID();
         block3dArray[ic].localID = block->getLocalID();
         block3dArray[ic].level = block->getLevel();
         block3dArray[ic].interpolationFlagCF = block->getCollectionOfInterpolationFlagCF();
         block3dArray[ic].interpolationFlagFC = block->getCollectionOfInterpolationFlagFC();
         block3dArray[ic].counter = block->getMaxGlobalID();
         block3dArray[ic].active = block->isActive();

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBlocks start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   // write to the file
   MPI_File file_handler;
   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(i nfo,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

   // if (comm->isRoot())
   // {
   UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);
   // }

   double start, finish;
   //MPI_Offset write_offset = (MPI_Offset)(size * sizeof(int));
   MPI_Offset write_offset = (MPI_Offset)(sizeof(int));

   if (comm->isRoot())
   {
      start = MPI_Wtime();

      // each process writes the quantity of it's blocks
      MPI_File_write_at(file_handler, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      // each process writes parameters of the grid
      MPI_File_write_at(file_handler, write_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
      // each process writes it's blocks
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);
      //MPI_File_sync(file_handler);
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBlocks time: " << finish - start << " s");
   }

   delete[] block3dArray;
   delete gridParameters;
}

void MPIIOMigrationBECoProcessor::writeDataSet(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (comm->isRoot())
	std::cout << "size = "<<size<<std::endl;

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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   DSArraysPresence arrPresence;
   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
         CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
         CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
         CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();

         if (firstBlock)// && block->getKernel()) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

            if (localDistributions)
            {
               dataSetParamStr1.nx[0] = static_cast<int>(localDistributions->getNX1());
               dataSetParamStr1.nx[1] = static_cast<int>(localDistributions->getNX2());
               dataSetParamStr1.nx[2] = static_cast<int>(localDistributions->getNX3());
               dataSetParamStr1.nx[3] = static_cast<int>(localDistributions->getNX4());
            }

            if (nonLocalDistributions)
            {
               dataSetParamStr2.nx[0] = static_cast<int>(nonLocalDistributions->getNX1());
               dataSetParamStr2.nx[1] = static_cast<int>(nonLocalDistributions->getNX2());
               dataSetParamStr2.nx[2] = static_cast<int>(nonLocalDistributions->getNX3());
               dataSetParamStr2.nx[3] = static_cast<int>(nonLocalDistributions->getNX4());
            }
            if (zeroDistributions)
            {
               dataSetParamStr3.nx[0] = static_cast<int>(zeroDistributions->getNX1());
               dataSetParamStr3.nx[1] = static_cast<int>(zeroDistributions->getNX2());
               dataSetParamStr3.nx[2] = static_cast<int>(zeroDistributions->getNX3());
               dataSetParamStr3.nx[3] = 1;
            }

            // ... than save some parameters that are equal in all blocks
            dataSetParamStr1.nx1 = dataSetParamStr2.nx1 = dataSetParamStr3.nx1 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX1());
            dataSetParamStr1.nx2 = dataSetParamStr2.nx2 = dataSetParamStr3.nx2 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX2());
            dataSetParamStr1.nx3 = dataSetParamStr2.nx3 = dataSetParamStr3.nx3 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX3());

            doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
               dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
               dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray = block->getKernel()->getDataSet()->getAverageDensity();
            if (averageDensityArray)
               arrPresence.isAverageDensityArrayPresent = true;
            else
               arrPresence.isAverageDensityArrayPresent = false;

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
            if (AverageVelocityArray3DPtr)
               arrPresence.isAverageVelocityArrayPresent = true;
            else
               arrPresence.isAverageVelocityArrayPresent = false;

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
            if (AverageFluctArray3DPtr)
               arrPresence.isAverageFluktuationsArrayPresent = true;
            else
               arrPresence.isAverageFluktuationsArrayPresent = false;

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
            if (AverageTripleArray3DPtr)
               arrPresence.isAverageTripleArrayPresent = true;
            else
               arrPresence.isAverageTripleArrayPresent = false;

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
            if (ShearStressValArray3DPtr)
               arrPresence.isShearStressValArrayPresent = true;
            else
               arrPresence.isShearStressValArrayPresent = false;

            SPtr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
            if (relaxationFactor3DPtr)
               arrPresence.isRelaxationFactorPresent = true;
            else
               arrPresence.isRelaxationFactorPresent = false;

            firstBlock = false;
         }

         if (localDistributions && (dataSetParamStr1.nx[0]>0) && (dataSetParamStr1.nx[1]>0) && (dataSetParamStr1.nx[2]>0) && (dataSetParamStr1.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), localDistributions->getDataVector().begin(), localDistributions->getDataVector().end());
         if (nonLocalDistributions && (dataSetParamStr2.nx[0]>0) && (dataSetParamStr2.nx[1]>0) && (dataSetParamStr2.nx[2]>0) && (dataSetParamStr2.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), nonLocalDistributions->getDataVector().begin(), nonLocalDistributions->getDataVector().end());
         if (zeroDistributions && (dataSetParamStr3.nx[0]>0) && (dataSetParamStr3.nx[1]>0) && (dataSetParamStr3.nx[2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), zeroDistributions->getDataVector().begin(), zeroDistributions->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN_LUSTRE
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   // write to the file
   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(3 * sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   MPI_File_write_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeDataSet time: " << finish - start << " s");
   }

   MPI_File file_handler1;
   std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
   rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler1);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename1);
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

void MPIIOMigrationBECoProcessor::writeAverageDensityArray(int step)
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
   std::vector<double> doubleValuesArray; // double-values of the AverageDensityArray in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageDensityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray = block->getKernel()->getDataSet()->getAverageDensity();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(averageDensityArray->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(averageDensityArray->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(averageDensityArray->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(averageDensityArray->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) && (dataSetParamStr.nx[3] > 0))
            doubleValuesArray.insert(doubleValuesArray.end(), averageDensityArray->getDataVector().begin(), averageDensityArray->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageDensityArray time: " << finish - start << " s");
   }
}

void MPIIOMigrationBECoProcessor::writeAverageVelocityArray(int step)
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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageVelocityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageVelocityArray3DPtr->getDataVector().begin(), AverageVelocityArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageVelocityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageVelocityArray time: " << finish - start << " s");
   }
}

void MPIIOMigrationBECoProcessor::writeAverageFluktuationsArray(int step)
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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageFluktuationsArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageFluctArray3DPtr->getDataVector().begin(), AverageFluctArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageFluktuationsArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageFluktuationsArray time: " << finish - start << " s");
   }
}

void MPIIOMigrationBECoProcessor::writeAverageTripleArray(int step)
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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageTripleArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageTripleArray3DPtr->getDataVector().begin(), AverageTripleArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageTripleArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeAverageTripleArray time: " << finish - start << " s");
   }
}

void MPIIOMigrationBECoProcessor::writeShearStressValArray(int step)
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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeShearStressValArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), ShearStressValArray3DPtr->getDataVector().begin(), ShearStressValArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeShearStressValArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeShearStressValArray time: " << finish - start << " s");
   }
}

void MPIIOMigrationBECoProcessor::writeRelaxationFactor(int step)
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
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeRelaxationFactor start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         SPtr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();

         if (firstBlock) // when first (any) valid block...
         {
            firstGlobalID = block->getGlobalID();

            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(relaxationFactor3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(relaxationFactor3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(relaxationFactor3DPtr->getNX3());
            dataSetParamStr.nx[3] = 1;
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

            firstBlock = false;
         }

         if ((dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), relaxationFactor3DPtr->getDataVector().begin(), relaxationFactor3DPtr->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeRelaxationFactor start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;

#ifdef HLRN
   MPI_Info_create(&info);
   MPI_Info_set(info, "striping_factor", "40");
   MPI_Info_set(info, "striping_unit", "4M");
#endif

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(firstGlobalID) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeRelaxationFactor time: " << finish - start << " s");
   }
}

//---------------------------------------------------------------------------------

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

   int blocksCount = 0;    // quantity of blocks, that belong to this process
   size_t allBytesCount = 0;  // quantity of bytes, that one process writes to the file
   size_t count_boundCond = 0;	// how many BoundaryConditions in all blocks
   int count_indexContainer = 0;	// how many indexContainer-values in all blocks

   std::vector<SPtr<Block3D>> blocksVector[25];
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      grid->getBlocks(level, rank, blocksVector[level]);
      blocksCount += static_cast<int>(blocksVector[level].size());
   }

   BCAddMigration* bcAddArray = new BCAddMigration[blocksCount];
   size_t* bytesCount = new size_t[blocksCount];  // quantity of bytes, that each block writes to the file
   std::vector<BoundaryCondition>* bcVector = new std::vector<BoundaryCondition>[blocksCount];
   std::vector<int>* indexContainerVector = new std::vector<int>[blocksCount];
   std::vector<int> bcindexmatrixVector;

   bool bcindexmatrixCountNotInit = true;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  // all the blocks of the current level
      {
         SPtr<BCArray3D> bcArr = block->getKernel()->getBCProcessor()->getBCArray();

         bcAddArray[ic].globalID = block->getGlobalID(); // id of the block needed to find it while regenerating the grid
         bcAddArray[ic].boundCond_count = 0;             // how many BoundaryConditions in this block
         bcAddArray[ic].indexContainer_count = 0;        // how many indexContainer-values in this block
         bytesCount[ic] = sizeof(BCAddMigration);
         bcVector[ic].resize(0);
         indexContainerVector[ic].resize(0);

         for (int bc = 0; bc<bcArr->getBCVectorSize(); bc++)
         {
            BoundaryCondition* bouCond = new BoundaryCondition();
            if (bcArr->bcvector[bc] == NULL)
               memset(bouCond, 0, sizeof(BoundaryCondition));
            else
            {
               bouCond->noslipBoundaryFlags = bcArr->bcvector[bc]->getNoSlipBoundary();
               bouCond->slipBoundaryFlags = bcArr->bcvector[bc]->getSlipBoundary();
               bouCond->velocityBoundaryFlags = bcArr->bcvector[bc]->getVelocityBoundary();
               bouCond->densityBoundaryFlags = bcArr->bcvector[bc]->getDensityBoundary();
               bouCond->wallModelBoundaryFlags = bcArr->bcvector[bc]->getWallModelBoundary();
               bouCond->bcVelocityX1 = bcArr->bcvector[bc]->getBoundaryVelocityX1();
               bouCond->bcVelocityX2 = bcArr->bcvector[bc]->getBoundaryVelocityX2();
               bouCond->bcVelocityX3 = bcArr->bcvector[bc]->getBoundaryVelocityX3();
               bouCond->bcDensity = bcArr->bcvector[bc]->getBoundaryDensity();
               bouCond->bcLodiDensity = bcArr->bcvector[bc]->getDensityLodiDensity();
               bouCond->bcLodiVelocityX1 = bcArr->bcvector[bc]->getDensityLodiVelocityX1();
               bouCond->bcLodiVelocityX2 = bcArr->bcvector[bc]->getDensityLodiVelocityX2();
               bouCond->bcLodiVelocityX3 = bcArr->bcvector[bc]->getDensityLodiVelocityX3();
               bouCond->bcLodiLentgh = bcArr->bcvector[bc]->getDensityLodiLength();
               bouCond->nx1 = bcArr->bcvector[bc]->nx1;
               bouCond->nx2 = bcArr->bcvector[bc]->nx2;
               bouCond->nx3 = bcArr->bcvector[bc]->nx3;
               for (int iq = 0; iq<26; iq++)
                  bouCond->q[iq] = bcArr->bcvector[bc]->getQ(iq);
               bouCond->algorithmType = bcArr->bcvector[bc]->getBcAlgorithmType();
            }

            bcVector[ic].push_back(*bouCond);
            bcAddArray[ic].boundCond_count++;
            count_boundCond++;
            bytesCount[ic] += sizeof(BoundaryCondition);
         }

         if (bcindexmatrixCountNotInit)
         {
            boundCondParamStr.nx1 = static_cast<int>(bcArr->bcindexmatrix.getNX1());
            boundCondParamStr.nx2 = static_cast<int>(bcArr->bcindexmatrix.getNX2());
            boundCondParamStr.nx3 = static_cast<int>(bcArr->bcindexmatrix.getNX3());
            boundCondParamStr.bcindexmatrixCount = static_cast<int>(bcArr->bcindexmatrix.getDataVector().size());
            bcindexmatrixCountNotInit = false;
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

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::writeBoundaryConds start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_Offset write_offset = (MPI_Offset)(sizeof(int)) + (MPI_Offset)(bcAddArray[0].globalID) * (MPI_Offset)(boundCondParamStr.bcindexmatrixCount) * (MPI_Offset)(sizeof(int));

   MPI_File_write_at(file_handler, 0, &boundCondParamStr.bcindexmatrixCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, write_offset, &bcindexmatrixVector[0], blocksCount * boundCondParamStr.bcindexmatrixCount, MPI_INT, MPI_STATUS_IGNORE);

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC2.bin";
   rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

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
   if (comm->isRoot()) UBLOG(logINFO, "MPIIOMigrationBECoProcessor restart step: " << step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - start");

   readBlocks(step);
   SPtr<Grid3DVisitor> newMetisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
   grid->accept(newMetisVisitor);

   readDataSet(step);
   readBoundaryConds(step);

   grid->setTimeStep(step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - end");
}

void MPIIOMigrationBECoProcessor::readBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBlocks start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

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
   for (int n = 0; n<blocksCount; n++)
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

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBlocks end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::blocksExchange(int tagN, int ind1, int ind2, int doubleCountInBlock, std::vector<double>& pV, std::vector<double>* rawDataReceive)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   int indexB = ind1;
   int indexE = ind2;
   int myBlocksCount = indexE - indexB;

   std::vector<double>* rawDataSend = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataSend[r].resize(0);
      rawDataSend[r].push_back(0);
   }

   SPtr<Block3D> tempBlock;
   int tempRank;
   for (int ind = indexB - indexB; ind < indexE - indexB; ind++)
   {
      tempBlock = grid->getBlock(indexB + ind);
      if(!tempBlock)  throw UbException(UB_EXARGS,"MPIIOMigrationBECoProcessor::blocksExchange -- null block pointer!!!" );

      tempRank = tempBlock->getRank();

      if (tempRank == rank) // no need to send data, the process already has it
      {
         rawDataReceive[tempRank][0]++;
         rawDataReceive[tempRank].push_back(double(indexB + ind));
         rawDataReceive[tempRank].insert(rawDataReceive[tempRank].end(), pV.begin() + ind * doubleCountInBlock,
            pV.begin() + ind * doubleCountInBlock + doubleCountInBlock);
      }
      else  // we must send data to other processes
      {
         rawDataSend[tempRank][0]++;
         rawDataSend[tempRank].push_back(double(indexB + ind));
         rawDataSend[tempRank].insert(rawDataSend[tempRank].end(), pV.begin() + ind * doubleCountInBlock,
            pV.begin() + ind * doubleCountInBlock + doubleCountInBlock);
      }
   }

   MPI_Request* requests = new MPI_Request[size * 2]; // send + receive
   int requestCount = 0;
   MPI_Status status;
   int quant;
   int doubleBlockCount;
   int rds;

   for (int r = 0; r < size; r++)
   {
      if (r != rank)
      {
		 rds = rawDataSend[r].size();
         doubleBlockCount = (int)(rds / SEND_BLOCK_SIZE);
         if (doubleBlockCount * SEND_BLOCK_SIZE < rds)
            doubleBlockCount += 1;

	     for (int i = rds; i < doubleBlockCount * SEND_BLOCK_SIZE; i++)
	         rawDataSend[r].push_back(0);

         MPI_Isend(&rawDataSend[r][0], doubleBlockCount, sendBlockDoubleType, r, tagN, MPI_COMM_WORLD, &requests[requestCount]);
         requestCount++;
      }
   }

   for (int r = 0; r < size; r++)
   {
      if (r != rank)
      {
         MPI_Probe(r, tagN, MPI_COMM_WORLD, &status);
         MPI_Get_count(&status, sendBlockDoubleType, &quant);
         rawDataReceive[r].resize(quant * SEND_BLOCK_SIZE);
         MPI_Irecv(&rawDataReceive[r][0], quant, sendBlockDoubleType, r, tagN, MPI_COMM_WORLD, &requests[requestCount]);
         requestCount++;
      }
   }

   MPI_Waitall(requestCount, &requests[0], MPI_STATUSES_IGNORE);
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

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

  if (comm->isRoot())
	std::cout << "size = "<<size<<std::endl;

   dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
      dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
      dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks 

   MPI_Offset read_offset = (MPI_Offset)(3 * sizeof(dataSetParam)) + (MPI_Offset)(indexB * doubleCountInBlock * sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   //-------------------------------------- restore blocks ---------------------------------
   int blockID;
   std::vector<double> vectorsOfValues1, vectorsOfValues2, vectorsOfValues3;

   size_t vectorSize1 = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3];
   size_t vectorSize2 = dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3];
   size_t vectorSize3 = dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

   size_t index;
   for (int r = 0; r < size; r++)
   {
      index = 1;
      for (int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues1.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + vectorSize1);
         index += vectorSize1;

         vectorsOfValues2.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + vectorSize2);
         index += vectorSize2;

         vectorsOfValues3.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + vectorSize3);
         index += vectorSize3;

         SPtr<DistributionArray3D> mFdistributions(new D3Q27EsoTwist3DSplittedVector());

         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues1, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues2, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues3, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(dataSetParamStr1.nx1);
         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(dataSetParamStr1.nx2);
         dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(dataSetParamStr1.nx3);

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         this->lbmKernel->setBlock(block);
         SPtr<LBMKernel> kernel = this->lbmKernel->clone();
         LBMReal collFactor = LBMSystem::calcCollisionFactor(this->nue, block->getLevel());
         kernel->setCollisionFactor(collFactor);
         kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
         kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
         SPtr<DataSet3D> dataSetPtr = SPtr<DataSet3D>(new DataSet3D());
         dataSetPtr->setFdistributions(mFdistributions);
         kernel->setDataSet(dataSetPtr);
         block->setKernel(kernel);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readDataSet end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   //-------------------------------------------------------------

   DSArraysPresence arrPresence;
   MPI_File file_handler1;
   std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
   rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler1);
   if (rc != MPI_SUCCESS) return;// throw UbException(UB_EXARGS, "couldn't open file " + filename1);

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

void MPIIOMigrationBECoProcessor::readAverageDensityArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 1, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   //----------------------------- restore data ---------------------------------
   int blockID;
   std::vector<double> vectorsOfValues;
   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill mAverageDensity arrays
         SPtr<AverageValuesArray3D> mAverageDensity;
         mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setAverageDensity(mAverageDensity);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageDensityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::readAverageVelocityArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 2, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blockID;
   std::vector<double> vectorsOfValues;

   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill mAverageVelocity array
         SPtr<AverageValuesArray3D> mAverageVelocity;
         mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setAverageVelocity(mAverageVelocity);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageVelocityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::readAverageFluktuationsArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 3, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blockID;
   std::vector<double> vectorsOfValues;

   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill AverageFluktuations array
         SPtr<AverageValuesArray3D> mAverageFluktuations;
         mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setAverageFluctuations(mAverageFluktuations);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageFluktuationsArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::readAverageTripleArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int myBlocksCount;
   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 4, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blockID;
   std::vector<double> vectorsOfValues;

   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill AverageTriplecorrelations array
         SPtr<AverageValuesArray3D> mAverageTriplecorrelations;
         mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setAverageTriplecorrelations(mAverageTriplecorrelations);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readAverageTripleArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::readShearStressValArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 5, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blockID;
   std::vector<double> vectorsOfValues;

   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill ShearStressValuesArray array
         SPtr<ShearStressValuesArray3D> mShearStressValues;
         mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setShearStressValues(mShearStressValues);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readShearStressValArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

void MPIIOMigrationBECoProcessor::readRelaxationFactor(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process

   int myBlocksCount;
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);

   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(myBlocksCount * doubleCountInBlock); // double-values in all blocks

   MPI_Offset read_offset = (MPI_Offset)(sizeof(dataSetParam)) + (MPI_Offset)(indexB) * (MPI_Offset)(doubleCountInBlock) * (MPI_Offset)(sizeof(double));
   MPI_File_read_at(file_handler, read_offset, &doubleValuesArray[0], myBlocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<double>* rawDataReceive = new std::vector<double>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataReceive[r].push_back(0);
   }

   blocksExchange(MESSAGE_TAG + 6, indexB, indexE, doubleCountInBlock, doubleValuesArray, rawDataReceive);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor end of exchange of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blockID;
   std::vector<double> vectorsOfValues;

   size_t index;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];

   for (int r = 0; r < size; r++)
   {
      index = 1;
      for(int ii = 0; ii < rawDataReceive[r][0]; ii++)
      {
         blockID = (int)(rawDataReceive[r][index]);
         index += 1;

         vectorsOfValues.assign(rawDataReceive[r].data() + index, rawDataReceive[r].data() + index + nextVectorSize);
         index += nextVectorSize;

         // fill RelaxationFactor array
         SPtr<RelaxationFactorArray3D> mRelaxationFactor;
         mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));

         // find the nesessary block and fill it
         SPtr<Block3D> block = grid->getBlock(blockID);
         block->getKernel()->getDataSet()->setRelaxationFactor(mRelaxationFactor);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readRelaxationFactor end of restore of data, rank = " << rank);
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
   if (comm->isRoot()) start = MPI_Wtime();

   int blocksCountAll = grid->getNumberOfBlocks(); // quantity of all blocks in the grid
   int myBlocksCount;
   int blocksPerProcess = blocksCountAll / size;   // how many blocks has each process
   
   if (rank < (size - 1))
      myBlocksCount = blocksPerProcess;
   else
      myBlocksCount = blocksPerProcess + (blocksCountAll - blocksPerProcess * size);
   
   int indexB = rank * blocksPerProcess;  // the first "my" block
   int indexE = indexB + myBlocksCount;   // the latest "my" block
   
   std::vector<int> bcindexmatrixVAll;

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC1.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   int sizeOfBIM;
   MPI_File_read_at(file_handler, (MPI_Offset)0, &sizeOfBIM, 1, MPI_INT, MPI_STATUS_IGNORE);
   bcindexmatrixVAll.resize(myBlocksCount * sizeOfBIM);
   
   MPI_Offset read_offset = (MPI_Offset)(sizeof(int)) + (MPI_Offset)(indexB) * (MPI_Offset)(sizeOfBIM) * (MPI_Offset)(sizeof(int));
   MPI_File_read_at(file_handler, read_offset, &bcindexmatrixVAll[0], myBlocksCount * sizeOfBIM, MPI_INT, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds start of exchange of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   std::vector<int>* rawDataReceive = new std::vector<int>[size];
   std::vector<int>* rawDataSend = new std::vector<int>[size];
   for (int r = 0; r < size; r++)
   {
      rawDataReceive[r].resize(0);
      rawDataSend[r].resize(0);
      rawDataReceive[r].push_back(0);
      rawDataSend[r].push_back(0);
   }

   SPtr<Block3D> tempBlock;
   int tempRank;
   for (int ind = indexB - indexB; ind < indexE - indexB; ind++)
   {
      tempBlock = grid->getBlock(indexB + ind);
      tempRank = tempBlock->getRank();

      if (tempRank == rank) // no need to send data, the process already has it
      {
         rawDataReceive[tempRank][0]++;
         rawDataReceive[tempRank].push_back(indexB + ind);
         rawDataReceive[tempRank].insert(rawDataReceive[tempRank].end(), bcindexmatrixVAll.begin() + ind * sizeOfBIM,
            bcindexmatrixVAll.begin() + ind * sizeOfBIM + sizeOfBIM);
      }
      else  // we must send data to other processes
      {
         rawDataSend[tempRank][0]++;
         rawDataSend[tempRank].push_back(indexB + ind);
         rawDataSend[tempRank].insert(rawDataSend[tempRank].end(), bcindexmatrixVAll.begin() + ind * sizeOfBIM,
            bcindexmatrixVAll.begin() + ind * sizeOfBIM + sizeOfBIM);
      }
   }

   MPI_Request* requests = new MPI_Request[size * 2]; // send + receive
   int requestCount = 0;
   MPI_Status status;
   int quant;
   int intBlockCount;
   int rds;

   for (int r = 0; r < size; r++)
   {
      if (r != rank)
      {
 		 rds = rawDataSend[r].size();
         intBlockCount = (int)(rds / SEND_BLOCK_SIZE);
         if (intBlockCount * SEND_BLOCK_SIZE < rds)
            intBlockCount += 1;

	     for (int i = rds; i < intBlockCount * SEND_BLOCK_SIZE; i++)
	         rawDataSend[r].push_back(0);

         MPI_Isend(&rawDataSend[r][0], intBlockCount, sendBlockIntType, r, MESSAGE_TAG + 7, MPI_COMM_WORLD, &requests[requestCount]);
         //MPI_Isend(&rawDataSend[r][0], rawDataSend[r].size(), MPI_INT, r, MESSAGE_TAG + 7, MPI_COMM_WORLD, &requests[requestCount]);
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
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &boundCondParamStr, 4, MPI_INT, MPI_STATUS_IGNORE);

   int blockID;
   size_t index;
   MPI_Offset read_offset1, read_offset2;

   BCAddMigration bcAddArray;
   BoundaryCondition* nullBouCond = new BoundaryCondition();
   memset(nullBouCond, 0, sizeof(BoundaryCondition));
   BoundaryCondition* bcArray;
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
            MPI_File_read_at(file_handler, read_offset2 + (MPI_Offset)(sizeof(BCAddMigration)), &bcArray[0], bcAddArray.boundCond_count, boundCondType, MPI_STATUS_IGNORE);

         if (bcAddArray.indexContainer_count > 0)
            MPI_File_read_at(file_handler, read_offset2 + (MPI_Offset)(sizeof(BCAddMigration)) + (MPI_Offset)(bcAddArray.boundCond_count) * (MPI_Offset)(sizeof(BoundaryCondition)),
               &indexContainerV[0], bcAddArray.indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);

         bcVector.resize(0);
            
         for (size_t ibc = 0; ibc<bcAddArray.boundCond_count; ibc++)
         {
            SPtr<BoundaryConditions> bc;
            if (memcmp(&bcArray[ibc], nullBouCond, sizeof(BoundaryCondition)) == 0)
               bc = SPtr<BoundaryConditions>();
            else
            {
               bc = SPtr<BoundaryConditions>(new BoundaryConditions);
               bc->noslipBoundaryFlags = bcArray[ibc].noslipBoundaryFlags;
               bc->slipBoundaryFlags = bcArray[ibc].slipBoundaryFlags;
               bc->densityBoundaryFlags = bcArray[ibc].densityBoundaryFlags;
               bc->velocityBoundaryFlags = bcArray[ibc].velocityBoundaryFlags;
               bc->wallModelBoundaryFlags = bcArray[ibc].wallModelBoundaryFlags;
               bc->bcVelocityX1 = bcArray[ibc].bcVelocityX1;
               bc->bcVelocityX2 = bcArray[ibc].bcVelocityX2;
               bc->bcVelocityX3 = bcArray[ibc].bcVelocityX3;
               bc->bcDensity = bcArray[ibc].bcDensity;
               bc->bcLodiDensity = bcArray[ibc].bcLodiDensity;
               bc->bcLodiVelocityX1 = bcArray[ibc].bcLodiVelocityX1;
               bc->bcLodiVelocityX2 = bcArray[ibc].bcLodiVelocityX2;
               bc->bcLodiVelocityX3 = bcArray[ibc].bcLodiVelocityX3;
               bc->bcLodiLentgh = bcArray[ibc].bcLodiLentgh;

               bc->nx1 = bcArray[ibc].nx1;
               bc->nx2 = bcArray[ibc].nx2;
               bc->nx3 = bcArray[ibc].nx3;
               for (int iq = 0; iq<26; iq++)
                  bc->setQ(bcArray[ibc].q[iq], iq);
               bc->setBcAlgorithmType(bcArray[ibc].algorithmType);
            }

            bcVector.push_back(bc);
         }

         CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);
         SPtr<Block3D> block1 = grid->getBlock(blockID);

         SPtr<BCProcessor> bcProc = bcProcessor->clone(block1->getKernel());
         SPtr<BCArray3D> bcArr(new BCArray3D());
         bcArr->bcindexmatrix = bcim;
         bcArr->bcvector = bcVector;
         bcArr->indexContainer = indexContainerV;
         bcProc->setBCArray(bcArr);

         block1->getKernel()->setBCProcessor(bcProc);
      }
   }
  
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   MPI_File_close(&file_handler);

   delete nullBouCond;

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds end of restore of data, rank = " << rank);
      UBLOG(logINFO, "MPIIOMigrationBECoProcessor::readBoundaryConds time: " << finish - start << " s");
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setLBMKernel(SPtr<LBMKernel> kernel)
{
   this->lbmKernel = kernel;
}
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setBCProcessor(SPtr<BCProcessor> bcProcessor)
{
   this->bcProcessor = bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationBECoProcessor::setNu(double nu)
{
   this->nue = nu;
}

