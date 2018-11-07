#include "MPIIOMigrationCoProcessor.h"
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
#include "RenumberBlockVisitor.h"
#include "MetisPartitioningGridVisitor.h"
#include "PointerDefinitions.h"
#include "UbFileOutputASCII.h"
#include "UbFileInputASCII.h"

MPIIOMigrationCoProcessor::MPIIOMigrationCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
   const std::string& path,
   SPtr<Communicator> comm) :
   CoProcessor(grid, s),
   path(path),
   comm(comm)
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

   MPI_Type_contiguous(7, MPI_INT, &dataSetParamType);
   MPI_Type_commit(&dataSetParamType);

   //-----------------------------------------------------------------------

   MPI_Datatype typesDataSet[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksDataSet[3] = { 2, 2, 2 };
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

   //---------------------------------------

   MPI_Type_contiguous(3, MPI_INT, &boundCondTypeAdd);
   MPI_Type_commit(&boundCondTypeAdd);
   //---------------------------------------

   MPI_Type_contiguous(6, MPI_CHAR, &arrayPresenceType);
   MPI_Type_commit(&arrayPresenceType);

}
//////////////////////////////////////////////////////////////////////////
MPIIOMigrationCoProcessor::~MPIIOMigrationCoProcessor()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&dataSetParamType);
   MPI_Type_free(&block3dType);
   MPI_Type_free(&dataSetType);
   MPI_Type_free(&dataSetSmallType);
   MPI_Type_free(&boundCondParamType);
   MPI_Type_free(&boundCondType);
   MPI_Type_free(&boundCondTypeAdd);
   MPI_Type_free(&arrayPresenceType);
}

//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      if (comm->isRoot()) UBLOG(logINFO, "MPIIOMigrationCoProcessor save step: " << step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - start");
      /*if (comm->isRoot())*/ clearAllFiles((int)step);

      writeBlocks((int)step);
      writeDataSet((int)step);
      writeBoundaryConds((int)step);

      writeCpTimeStep(step);

      if (comm->isRoot()) UBLOG(logINFO, "Save check point - end");
   }
}
////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::clearAllFiles(int step)
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
   MPI_File_delete(filename4.c_str(), info);
   //int rc4 = MPI_File_open(MPI_COMM_WORLD, filename4.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc4 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename4);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename5 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   MPI_File_delete(filename5.c_str(), info);
   //int rc5 = MPI_File_open(MPI_COMM_WORLD, filename5.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc5 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename5);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename6 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   MPI_File_delete(filename6.c_str(), info);
   //int rc6 = MPI_File_open(MPI_COMM_WORLD, filename6.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc6 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename6);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename7 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   MPI_File_delete(filename7.c_str(), info);
   //int rc7 = MPI_File_open(MPI_COMM_WORLD, filename7.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc7 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename7);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename8 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   MPI_File_delete(filename8.c_str(), info);
   //int rc8 = MPI_File_open(MPI_COMM_WORLD, filename8.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc8 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename8);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename9 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   MPI_File_delete(filename9.c_str(), info);
   //int rc9 = MPI_File_open(MPI_COMM_WORLD, filename9.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   //if (rc9 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename9);
   //MPI_File_set_size(file_handler, new_size);
   //MPI_File_close(&file_handler);

   std::string filename10 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   int rc10 = MPI_File_open(MPI_COMM_WORLD, filename10.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc10 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename10);
   MPI_File_set_size(file_handler, new_size);
   MPI_File_close(&file_handler);
}
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::writeCpTimeStep(int step)
{
   if (comm->isRoot())
   {
      UbFileOutputASCII f(path + "/mpi_io_cp/cp.txt");
      f.writeInteger(step);
   }
}
//////////////////////////////////////////////////////////////////////////
int MPIIOMigrationCoProcessor::readCpTimeStep()
{
   UbFileInputASCII f(path + "/mpi_io_cp/cp.txt");
   int step = f.readInteger();
   return step;
}
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::writeBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   grid->renumberBlockIDs();

   //grid->deleteBlockIDs();
   //RenumberBlockVisitor renumber;
   //grid->accept(renumber);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBlocks start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   int blocksCount = 0; // quantity of all the blocks in the grid, max 2147483648 blocks!
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   std::vector<SPtr<Block3D>> blocksVector[25]; // max 25 levels
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      //grid->getBlocks(level, rank, blockVector[level]);
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
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBlocks start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   // write to the file
   MPI_File file_handler;
   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

   // if (comm->isRoot())
   // {
   UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);
   // }

   double start, finish;
   MPI_Offset write_offset = (MPI_Offset)(size * sizeof(int));

   if (comm->isRoot())
   {
      start = MPI_Wtime();

      // each process writes the quantity of it's blocks
      MPI_File_write_at(file_handler, 0/*rank*sizeof(int)*/, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
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
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeBlocks time: " << finish - start << " s");
   }

   delete[] block3dArray;
   delete gridParameters;
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
   DataSet* dataSetArray = new DataSet[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeDataSet start collect data rank = " << rank);
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
         dataSetArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid
         dataSetArray[ic].ghostLayerWidth = block->getKernel()->getGhostLayerWidth();
         dataSetArray[ic].collFactor = block->getKernel()->getCollisionFactor();
         dataSetArray[ic].deltaT = block->getKernel()->getDeltaT();
         dataSetArray[ic].compressible = block->getKernel()->getCompressible();
         dataSetArray[ic].withForcing = block->getKernel()->getWithForcing();

         SPtr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
         CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
         CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
         CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();

         if (firstBlock)// && block->getKernel()) // when first (any) valid block...
         {
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

            SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray = block->getKernel()->getDataSet()->getAverageDencity();
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

   // register new MPI-type depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeDataSet start MPI IO rank = " << rank);
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

   MPI_File_write_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSet) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(3 * sizeof(dataSetParam) + dataSetArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetArray[nb], 1, dataSetType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSet)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   std::vector<double> doubleValuesArray; // double-values of the AverageDensityArray in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray = block->getKernel()->getDataSet()->getAverageDencity();

         if (firstBlock) // when first (any) valid block...
         {
            //if (averageDensityArray)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(averageDensityArray->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(averageDensityArray->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(averageDensityArray->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(averageDensityArray->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (averageDensityArray && (dataSetParamStr.nx[0] > 0) && (dataSetParamStr.nx[1] > 0) && (dataSetParamStr.nx[2] > 0) && (dataSetParamStr.nx[3] > 0))
            doubleValuesArray.insert(doubleValuesArray.end(), averageDensityArray->getDataVector().begin(), averageDensityArray->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageVelocityArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageVelocityArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageVelocityArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageVelocityArray3DPtr->getDataVector().begin(), AverageVelocityArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageVelocityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageFluktuationsArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageFluctArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageFluctArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageFluctArray3DPtr->getDataVector().begin(), AverageFluctArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageFluktuationsArray start MPI IO rank = " << rank);
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

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
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
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();

         if (firstBlock) // when first (any) valid block...
         {
            //if (AverageTripleArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (AverageTripleArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageTripleArray3DPtr->getDataVector().begin(), AverageTripleArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeAverageTripleArray start MPI IO rank = " << rank);
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

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeShearStressValArray start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();

         if (firstBlock) // when first (any) valid block...
         {
            //if (ShearStressValArray3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
            dataSetParamStr.nx[3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (ShearStressValArray3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0) && (dataSetParamStr.nx[3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), ShearStressValArray3DPtr->getDataVector().begin(), ShearStressValArray3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeShearStressValArray start MPI IO rank = " << rank);
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

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 
   dataSetParam dataSetParamStr;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeRelaxationFactor start collect data rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   bool firstBlock = true;
   int doubleCountInBlock = 0;
   int ic = 0;
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         dataSetSmallArray[ic].globalID = block->getGlobalID();     // id of the block needed to find it while regenerating the grid

         SPtr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();

         if (firstBlock) // when first (any) valid block...
         {
            //if (relaxationFactor3DPtr)
            //{
            dataSetParamStr.nx1 = dataSetParamStr.nx2 = dataSetParamStr.nx3 = 0;
            dataSetParamStr.nx[0] = static_cast<int>(relaxationFactor3DPtr->getNX1());
            dataSetParamStr.nx[1] = static_cast<int>(relaxationFactor3DPtr->getNX2());
            dataSetParamStr.nx[2] = static_cast<int>(relaxationFactor3DPtr->getNX3());
            dataSetParamStr.nx[3] = 1;
            doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
            //}
            //else
            //   break;

            firstBlock = false;
         }

         if (relaxationFactor3DPtr && (dataSetParamStr.nx[0]>0) && (dataSetParamStr.nx[1]>0) && (dataSetParamStr.nx[2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), relaxationFactor3DPtr->getDataVector().begin(), relaxationFactor3DPtr->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::writeRelaxationFactor start MPI IO rank = " << rank);
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

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

   MPI_Offset write_offset;
   size_t sizeofOneDataSet = sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double);

   for (int nb = 0; nb < blocksCount; nb++)
   {
      write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallArray[nb].globalID * sizeofOneDataSet);
      MPI_File_write_at(file_handler, write_offset, &dataSetSmallArray[nb], 1, dataSetSmallType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(DataSetSmall)), &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
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

   BCAdd* bcAddArray = new BCAdd[blocksCount];
   size_t* bytesCount = new size_t[blocksCount];  // quantity of bytes, that each block writes to the file
   std::vector<BoundaryCondition>* bcVector = new std::vector<BoundaryCondition>[blocksCount];
   std::vector<int>* bcindexmatrixVector = new std::vector<int>[blocksCount];
   std::vector<int>* indexContainerVector = new std::vector<int>[blocksCount];

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
         bytesCount[ic] = sizeof(BCAdd);
         bcVector[ic].resize(0);
         bcindexmatrixVector[ic].resize(0);
         indexContainerVector[ic].resize(0);

         for (int bc = 0; bc<bcArr->getBCVectorSize(); bc++)
         {
            BoundaryCondition* bouCond = new BoundaryCondition();
            if (bcArr->bcvector[bc] == NULL)
            {
               memset(bouCond, 0, sizeof(BoundaryCondition));
            }
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

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
         MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAdd)), &bcVector[nb][0], bcAddArray[nb].boundCond_count, boundCondType, MPI_STATUS_IGNORE);

      if (bcindexmatrixVector[nb].size() > 0)
         MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAdd) + bcAddArray[nb].boundCond_count * sizeof(BoundaryCondition)),
            &bcindexmatrixVector[nb][0], 1, bcindexmatrixType, MPI_STATUS_IGNORE);

      if (indexContainerVector[nb].size() > 0)
         MPI_File_write_at(file_handler, (MPI_Offset)(write_offset + sizeof(BCAdd) + bcAddArray[nb].boundCond_count * sizeof(BoundaryCondition) + boundCondParamStr.bcindexmatrixCount * sizeof(int)),
            &indexContainerVector[nb][0], bcAddArray[nb].indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);

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
   if (comm->isRoot()) UBLOG(logINFO, "MPIIOMigrationCoProcessor restart step: " << step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - start");

   readBlocks(step);

   SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::KWAY));
   grid->accept(metisVisitor);

   readDataSet(step);
   readBoundaryConds(step);

   grid->setTimeStep(step);

   if (comm->isRoot()) UBLOG(logINFO, "Load check point - end");
   //this->reconnect(grid);
}

void MPIIOMigrationCoProcessor::readBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBlocks start MPI IO rank = " << rank);
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
   //MPI_File_read_at(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   Block3d* block3dArray = new Block3d[blocksCount];

   GridParam* gridParameters = new GridParam;

   // calculate the read offset
   MPI_Offset read_offset = (MPI_Offset)(size * sizeof(int));

   // read parameters of the grid
   MPI_File_read_at(file_handler, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   // read all the blocks
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBlocks time: " << finish - start << " s");
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBlocks start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   // clear the grid
   std::vector<SPtr<Block3D>> blocksVector;
   grid->getBlocks(0, blocksVector);
   int del = 0;
   for (SPtr<Block3D> block : blocksVector)
   {
      grid->deleteBlock(block);
      del++;
   }

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
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readBlocks end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
}

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
   if (comm->isRoot()) start = MPI_Wtime();

   int blocksCount = 0; // quantity of the blocks, that belong to this process
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

   DataSet* dataSetArray = new DataSet[blocksCount];

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   MPI_File_read_at(file_handler, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);

   double doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
      dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
      dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks 

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSet) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(3 * sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetArray[ic], 1, dataSetType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSet)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readDataSet start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0, vectorSize = 0;
   std::vector<double> vectorsOfValues1, vectorsOfValues2, vectorsOfValues3;

   for (int n = 0; n < blocksCount; n++)
   {
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

      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues1, dataSetParamStr1.nx[0], dataSetParamStr1.nx[1], dataSetParamStr1.nx[2], dataSetParamStr1.nx[3])));
      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues2, dataSetParamStr2.nx[0], dataSetParamStr2.nx[1], dataSetParamStr2.nx[2], dataSetParamStr2.nx[3])));
      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues3, dataSetParamStr3.nx[0], dataSetParamStr3.nx[1], dataSetParamStr3.nx[2])));

      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(dataSetParamStr1.nx1);
      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(dataSetParamStr1.nx2);
      dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(dataSetParamStr1.nx3);

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetArray[n].globalID);
      //std::cout << "rank="<<rank<<", dataSetArray[n].globalID=" << dataSetArray[n].globalID << std::endl;

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

   MPI_Type_free(&dataSetDoubleType);

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
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename1);
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

void MPIIOMigrationCoProcessor::readAverageDensityArray(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray start MPI IO rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageDensityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill mAverageDensity arrays
      SPtr<AverageValuesArray3D> mAverageDensity;
      //if ((dataSetParamStr.nx[0]==0)&&(dataSetParamStr.nx[1]==0)&&(dataSetParamStr.nx[2]==0)&&(dataSetParamStr.nx[3]==0))
      //   mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

      //std::cout << "rank=" << rank << ", dataSetArray[n].globalID=" << dataSetSmallArray[n].globalID << std::endl;
      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageDencity(mAverageDensity);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageDensityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageVelocityArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill mAverageVelocity array
      SPtr<AverageValuesArray3D> mAverageVelocity;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) && (dataSetParamStr.nx[3] == 0))
      //   mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageVelocity(mAverageVelocity);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageVelocityArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageFluktuationsArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill AverageFluktuations array
      SPtr<AverageValuesArray3D> mAverageFluktuations;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) && (dataSetParamStr.nx[3] == 0))
      //   mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageFluctuations(mAverageFluktuations);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageFluktuationsArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpAverageTripleArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill AverageTriplecorrelations array
      SPtr<AverageValuesArray3D> mAverageTriplecorrelations;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) && (dataSetParamStr.nx[3] == 0))
      //   mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setAverageTriplecorrelations(mAverageTriplecorrelations);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readAverageTripleArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpShearStressValArray.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill ShearStressValuesArray array
      SPtr<ShearStressValuesArray3D> mShearStressValues;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0) && (dataSetParamStr.nx[3] == 0))
      //   mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      //else
      mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2], dataSetParamStr.nx[3]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setShearStressValues(mShearStressValues);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readShearStressValArray end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
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
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpRelaxationFactor.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

   // read count of blocks
   int blocksCount = 0;
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

   DataSetSmall* dataSetSmallArray = new DataSetSmall[blocksCount];
   int doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> doubleValuesArray(blocksCount * doubleCountInBlock); // double-values in all blocks

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   int ic = 0;
   MPI_Offset read_offset;
   size_t sizeofOneDataSet = size_t(sizeof(DataSetSmall) + doubleCountInBlock * sizeof(double));

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset = (MPI_Offset)(sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         MPI_File_read_at(file_handler, read_offset, &dataSetSmallArray[ic], 1, dataSetSmallType, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset + sizeof(DataSetSmall)), &doubleValuesArray[ic * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         ic++;
      }
   }

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor time: " << finish - start << " s");
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor start of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   size_t index = 0;
   size_t nextVectorSize = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
   std::vector<double> vectorsOfValues;
   for (int n = 0; n < blocksCount; n++)
   {
      vectorsOfValues.assign(doubleValuesArray.data() + index, doubleValuesArray.data() + index + nextVectorSize);
      index += nextVectorSize;

      // fill RelaxationFactor array
      SPtr<RelaxationFactorArray3D> mRelaxationFactor;
      //if ((dataSetParamStr.nx[0] == 0) && (dataSetParamStr.nx[1] == 0) && (dataSetParamStr.nx[2] == 0))
      //   mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr();
      //else
      mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues, dataSetParamStr.nx[0], dataSetParamStr.nx[1], dataSetParamStr.nx[2]));

      // find the nesessary block and fill it
      SPtr<Block3D> block = grid->getBlock(dataSetSmallArray[n].globalID);
      block->getKernel()->getDataSet()->setRelaxationFactor(mRelaxationFactor);
   }

   MPI_Type_free(&dataSetDoubleType);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIOMigrationCoProcessor::readRelaxationFactor end of restore of data, rank = " << rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");
   }

   delete[] dataSetSmallArray;
}

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
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename);

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

   BCAdd* bcAddArray = new BCAdd[blocksCount];
   BoundaryCondition* nullBouCond = new BoundaryCondition();
   memset(nullBouCond, 0, sizeof(BoundaryCondition));
   BoundaryCondition* bcArray;
   int* intArray1;
   int* intArray2;
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
      for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
      {
         read_offset1 = (MPI_Offset)(sizeof(boundCondParam) + block->getGlobalID() * sizeof(size_t));

         MPI_File_read_at(file_handler, read_offset1, &read_offset2, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);
         MPI_File_read_at(file_handler, read_offset2, &bcAddArray[ic], 1, boundCondTypeAdd, MPI_STATUS_IGNORE);

         bcArray = new BoundaryCondition[bcAddArray[ic].boundCond_count];
         intArray1 = new int[boundCondParamStr.bcindexmatrixCount];
         intArray2 = new int[bcAddArray[ic].indexContainer_count];

         if (bcAddArray[ic].boundCond_count > 0)
         {
            MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAdd)), &bcArray[0], bcAddArray[ic].boundCond_count, boundCondType, MPI_STATUS_IGNORE);
         }
         MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAdd) + bcAddArray[ic].boundCond_count * sizeof(BoundaryCondition)),
            &intArray1[0], 1, bcindexmatrixType, MPI_STATUS_IGNORE);
         if (bcAddArray[ic].indexContainer_count > 0)
         {
            MPI_File_read_at(file_handler, (MPI_Offset)(read_offset2 + sizeof(BCAdd) + bcAddArray[ic].boundCond_count * sizeof(BoundaryCondition) + boundCondParamStr.bcindexmatrixCount * sizeof(int)),
               &intArray2[0], bcAddArray[ic].indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);
         }

         bcindexmatrixV.resize(0);
         indexContainerV.resize(0);
         bcVector.resize(0);

         for (size_t ibc = 0; ibc<bcAddArray[ic].boundCond_count; ibc++)
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

         for (int b1 = 0; b1 < boundCondParamStr.bcindexmatrixCount; b1++)
            bcindexmatrixV.push_back(intArray1[b1]);

         for (int b2 = 0; b2 < bcAddArray[ic].indexContainer_count; b2++)
            indexContainerV.push_back(intArray2[b2]);

         CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);
         SPtr<Block3D> block1 = grid->getBlock(bcAddArray[ic].globalID);

         SPtr<BCProcessor> bcProc = bcProcessor->clone(block1->getKernel());
         SPtr<BCArray3D> bcArr(new BCArray3D());
         bcArr->bcindexmatrix = bcim;
         bcArr->bcvector = bcVector;
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
void MPIIOMigrationCoProcessor::setLBMKernel(SPtr<LBMKernel> kernel)
{
   this->lbmKernel = kernel;
}
//////////////////////////////////////////////////////////////////////////
void MPIIOMigrationCoProcessor::setBCProcessor(SPtr<BCProcessor> bcProcessor)
{
   this->bcProcessor = bcProcessor;
}

