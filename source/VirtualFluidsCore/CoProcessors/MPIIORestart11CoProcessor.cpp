#include "MPIIORestart11CoProcessor.h"
#include <boost/foreach.hpp>
#include "D3Q27System.h"
//#include "LBMKernel.h"
#include "CompressibleCumulantLBMKernel.h"
#include "IncompressibleCumulantLBMKernel.h"
#include "ThinWallBCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <UbSystem.h>
#include <MemoryUtil.h>
#include "BoundaryConditions.h"
#include "Block3D.h"
#include "CoordinateTransformation3D.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "BCArray3D.h"

//! BLOCK_SIZE defines the quantity of the BoundaryCondition-structures written as one block to the file
//! To avoid overflow in the parameter \a count of the function MPI_File_write_at 
//! structures BoundaryCondition are being written in blocks containing each of them BLOCK_SIZE structures
#define BLOCK_SIZE 1024

MPIIORestart11CoProcessor::MPIIORestart11CoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
   const std::string& path,
   CommunicatorPtr comm) :
   CoProcessor(grid, s),
   path(path),
   comm(comm),
   mpiTypeFreeFlag(false)
{
   UbSystem::makeDirectory(path+"/mpi_io_cp");

   memset(&dataSetParamStr, 0, sizeof(dataSetParamStr));
   memset(&boundCondParamStr, 0, sizeof(boundCondParamStr));

   //-------------------------   define MPI types  ---------------------------------

   MPI_Datatype typesGP[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksGP[3] = { 34, 6, 5 };
   MPI_Aint offsetsGP[3], lbGP, extentGP;

   offsetsGP[0] = 0;
   MPI_Type_get_extent(MPI_DOUBLE, &lbGP, &extentGP);
   offsetsGP[1] = blocksGP[0]*extentGP;

   MPI_Type_get_extent(MPI_INT, &lbGP, &extentGP);
   offsetsGP[2] = offsetsGP[1]+blocksGP[1]*extentGP;

   MPI_Type_create_struct (3, blocksGP, offsetsGP, typesGP, &gridParamType);
   MPI_Type_commit(&gridParamType);

   //-----------------------------------------------------------------------

   MPI_Datatype typesBlock[2] = { MPI_INT, MPI_CHAR };
   int blocksBlock[2] = { 13, 1 };
   MPI_Aint offsetsBlock[2], lbBlock, extentBlock;

   offsetsBlock[0] = 0;
   MPI_Type_get_extent(MPI_INT, &lbBlock, &extentBlock);
   offsetsBlock[1] = blocksBlock[0]*extentBlock;

   MPI_Type_create_struct(2, blocksBlock, offsetsBlock, typesBlock, &block3dType);
   MPI_Type_commit(&block3dType);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(40, MPI_INT, &dataSetParamType);
   MPI_Type_commit(&dataSetParamType);

   //-----------------------------------------------------------------------

   MPI_Datatype typesDataSet[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksDataSet[3] = { 2, 5, 2 };
   MPI_Aint offsetsDatatSet[3], lbDataSet, extentDataSet;

   offsetsDatatSet[0] = 0;
   MPI_Type_get_extent(MPI_DOUBLE, &lbDataSet, &extentDataSet);
   offsetsDatatSet[1] = blocksDataSet[0]*extentDataSet;

   MPI_Type_get_extent(MPI_INT, &lbDataSet, &extentDataSet);
   offsetsDatatSet[2] = offsetsDatatSet[1]+blocksDataSet[1]*extentDataSet;

   MPI_Type_create_struct(3, blocksDataSet, offsetsDatatSet, typesDataSet, &dataSetType);
   MPI_Type_commit(&dataSetType);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(4, MPI_INT, &boundCondParamType);
   MPI_Type_commit(&boundCondParamType);

   //-----------------------------------------------------------------------

   MPI_Datatype typesBC[3] = { MPI_LONG_LONG_INT, MPI_FLOAT, MPI_CHAR };
   int blocksBC[3] = { 5, 38, 1 };
   MPI_Aint offsetsBC[3], lbBC, extentBC;

   offsetsBC[0] = 0;
   MPI_Type_get_extent(MPI_LONG_LONG_INT, &lbBC, &extentBC);
   offsetsBC[1] = blocksBC[0]*extentBC;

   MPI_Type_get_extent(MPI_FLOAT, &lbBC, &extentBC);
   offsetsBC[2] = offsetsBC[1]+blocksBC[1]*extentBC;

   MPI_Type_create_struct(3, blocksBC, offsetsBC, typesBC, &boundCondType);
   MPI_Type_commit(&boundCondType);

   //---------------------------------------

   MPI_Type_contiguous(BLOCK_SIZE, boundCondType, &boundCondType1000);
   MPI_Type_commit(&boundCondType1000);

   //---------------------------------------

   MPI_Type_contiguous(6, MPI_INT, &boundCondTypeAdd);
   MPI_Type_commit(&boundCondTypeAdd);

}
//////////////////////////////////////////////////////////////////////////
MPIIORestart11CoProcessor::~MPIIORestart11CoProcessor()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&block3dType);
   MPI_Type_free(&dataSetParamType);
   MPI_Type_free(&dataSetType);
   MPI_Type_free(&boundCondParamType);
   MPI_Type_free(&boundCondType);
   MPI_Type_free(&boundCondType1000);
   MPI_Type_free(&boundCondTypeAdd);

   if (mpiTypeFreeFlag)
   {
      MPI_Type_free(&dataSetDoubleType);
      MPI_Type_free(&bcindexmatrixType);
   }
}

//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      if (comm->isRoot()) UBLOG(logINFO, "MPIIORestart11CoProcessor save step: "<<step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - start");
      /*if (comm->isRoot())*/ clearAllFiles((int)step);
      writeBlocks((int)step);
      writeDataSet((int)step);
      writeBoundaryConds((int)step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - end");
      
      //readDataSet((int)step);
      //readBoundaryConds((int)step);
   }
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::clearAllFiles(int step)
{
   MPI_File file_handler1, file_handler2, file_handler3;
   MPI_Info info = MPI_INFO_NULL;
   MPI_Offset new_size = 0;

   UbSystem::makeDirectory(path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));
   std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   //MPI_File_delete(filename1.c_str(), info);
   int rc1 = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler1);
   if (rc1 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename1);
   MPI_File_set_size(file_handler1, new_size);
   //MPI_File_sync(file_handler1);
   MPI_File_close(&file_handler1);

   std::string filename2 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   //MPI_File_delete(filename2.c_str(), info);
   int rc2 = MPI_File_open(MPI_COMM_WORLD, filename2.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler2);
   if (rc2 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename2);
   MPI_File_set_size(file_handler2, new_size);
   //MPI_File_sync(file_handler2);
   MPI_File_close(&file_handler2);

   std::string filename3 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   //MPI_File_delete(filename3.c_str(), info);
   int rc3 = MPI_File_open(MPI_COMM_WORLD, filename3.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &file_handler3);
   if (rc3 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename3);
   MPI_File_set_size(file_handler3, new_size);
   //MPI_File_sync(file_handler3);
   MPI_File_close(&file_handler3);
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::writeBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size=1;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBlocks start collect data rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   std::vector<Block3DPtr> blocksVector[25]; // max 25 levels
   for (int level = minInitLevel; level<=maxInitLevel; level++)
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
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      for(Block3DPtr block : blocksVector[level])  //	all the blocks of the current level
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
         block3dArray[ic].interpolationFlagCF = block->getInterpolationFlagCF();
         block3dArray[ic].interpolationFlagFC = block->getInterpolationFlagFC();
         block3dArray[ic].counter = block->getMaxGlobalID();
         block3dArray[ic].active = block->isActive();

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBlocks start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

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
      MPI_File_write_at(file_handler, (MPI_Offset)(rank*sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      // each process writes parameters of the grid
      MPI_File_write_at(file_handler, write_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
      // each process writes it's blocks
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset +sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);
   }

   MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);
 
   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBlocks time: "<<finish-start<<" s");
   }

   delete[] block3dArray;
   delete gridParameters;
}

void MPIIORestart11CoProcessor::writeDataSet(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   int blocksCount = 0; // quantity of blocks in the grid, max 2147483648 blocks!

   std::vector<Block3DPtr> blocksVector[25];
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      grid->getBlocks(level, rank, blocksVector[level]);
      blocksCount += static_cast<int>(blocksVector[level].size());
   }

   DataSet* dataSetArray = new DataSet[blocksCount];
   std::vector<double> doubleValuesArray; // double-values (arrays of f's) in all blocks 

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeDataSet start collect data rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   bool firstBlock = true;
   int ic = 0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      for(Block3DPtr block : blocksVector[level])  //	blocks of the current level
      {
         dataSetArray[ic].x1 = block->getX1();     // coordinates of the block needed to find it while regenerating the grid
         dataSetArray[ic].x2 = block->getX2();
         dataSetArray[ic].x3 = block->getX3();
         dataSetArray[ic].level = block->getLevel();
         //if (block->getKernel())
         //{
         dataSetArray[ic].ghostLayerWidth = block->getKernel()->getGhostLayerWidth();
         dataSetArray[ic].collFactor = block->getKernel()->getCollisionFactor();
         dataSetArray[ic].deltaT = block->getKernel()->getDeltaT();
         dataSetArray[ic].compressible = block->getKernel()->getCompressible();
         dataSetArray[ic].withForcing = block->getKernel()->getWithForcing();
         //}
         //else
         //{
         //   dataSetArray[ic].ghostLayerWidth = 0;
         //   dataSetArray[ic].collFactor = 0.0;
         //   dataSetArray[ic].deltaT = 0.0;
         //   dataSetArray[ic].compressible = false;
         //   dataSetArray[ic].withForcing = false;
         //}
         //std::cout << "ic="<<ic<<"-"<<dataSetArray[ic].x1 << "," << dataSetArray[ic].x2 << "," << dataSetArray[ic].x3 << "," << dataSetArray[ic].level << "," << dataSetArray[ic].ghostLayerWidth;
         //std::cout << dataSetArray[ic].collFactor<<","<<dataSetArray[ic].deltaT<<","<<dataSetArray[ic].compressible<<","<<dataSetArray[ic].withForcing<<std::endl;
         //dataSetArrayGW[ic].x1 = dataSetArray[ic].x1;
         //dataSetArrayGW[ic].x2 = dataSetArray[ic].x2;
         //dataSetArrayGW[ic].x3 = dataSetArray[ic].x3;
         //dataSetArrayGW[ic].level = dataSetArray[ic].level;
         //dataSetArrayGW[ic].ghostLayerWidth = dataSetArray[ic].ghostLayerWidth;
         //dataSetArrayGW[ic].collFactor = dataSetArray[ic].collFactor;
         //dataSetArrayGW[ic].deltaT = dataSetArray[ic].deltaT;
         //dataSetArrayGW[ic].compressible = dataSetArray[ic].compressible;
         //dataSetArrayGW[ic].withForcing = dataSetArray[ic].withForcing;

         if (firstBlock /*&& block->getKernel()*/) // when first (any) valid block...
         {
            std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > averageDensityArray = block->getKernel()->getDataSet()->getAverageDencity();
            if (averageDensityArray)
            {
               dataSetParamStr.nx[0][0] = static_cast<int>(averageDensityArray->getNX1());
               dataSetParamStr.nx[0][1] = static_cast<int>(averageDensityArray->getNX2());
               dataSetParamStr.nx[0][2] = static_cast<int>(averageDensityArray->getNX3());
               dataSetParamStr.nx[0][3] = static_cast<int>(averageDensityArray->getNX4());
            }

            std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
            if (AverageVelocityArray3DPtr)
            {
               dataSetParamStr.nx[1][0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
               dataSetParamStr.nx[1][1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
               dataSetParamStr.nx[1][2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
               dataSetParamStr.nx[1][3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
            }

            std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
            if (AverageFluctArray3DPtr)
            {
               dataSetParamStr.nx[2][0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
               dataSetParamStr.nx[2][1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
               dataSetParamStr.nx[2][2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
               dataSetParamStr.nx[2][3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
            }

            std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
            if (AverageTripleArray3DPtr)
            {
               dataSetParamStr.nx[3][0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
               dataSetParamStr.nx[3][1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
               dataSetParamStr.nx[3][2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
               dataSetParamStr.nx[3][3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
            }

            std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
            if (ShearStressValArray3DPtr)
            {
               dataSetParamStr.nx[4][0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
               dataSetParamStr.nx[4][1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
               dataSetParamStr.nx[4][2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
               dataSetParamStr.nx[4][3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
            }

            std::shared_ptr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
            if (relaxationFactor3DPtr)
            {
               dataSetParamStr.nx[5][0] = static_cast<int>(relaxationFactor3DPtr->getNX1());
               dataSetParamStr.nx[5][1] = static_cast<int>(relaxationFactor3DPtr->getNX2());
               dataSetParamStr.nx[5][2] = static_cast<int>(relaxationFactor3DPtr->getNX3());
               dataSetParamStr.nx[5][3] = 1;
            }

            std::shared_ptr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
            CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
            if (localDistributions)
            {
               dataSetParamStr.nx[6][0] = static_cast<int>(localDistributions->getNX1());
               dataSetParamStr.nx[6][1] = static_cast<int>(localDistributions->getNX2());
               dataSetParamStr.nx[6][2] = static_cast<int>(localDistributions->getNX3());
               dataSetParamStr.nx[6][3] = static_cast<int>(localDistributions->getNX4());
            }

            CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
            if (nonLocalDistributions)
            {
               dataSetParamStr.nx[7][0] = static_cast<int>(nonLocalDistributions->getNX1());
               dataSetParamStr.nx[7][1] = static_cast<int>(nonLocalDistributions->getNX2());
               dataSetParamStr.nx[7][2] = static_cast<int>(nonLocalDistributions->getNX3());
               dataSetParamStr.nx[7][3] = static_cast<int>(nonLocalDistributions->getNX4());
            }

            CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();
            if (zeroDistributions)
            {
               dataSetParamStr.nx[8][0] = static_cast<int>(zeroDistributions->getNX1());
               dataSetParamStr.nx[8][1] = static_cast<int>(zeroDistributions->getNX2());
               dataSetParamStr.nx[8][2] = static_cast<int>(zeroDistributions->getNX3());
               dataSetParamStr.nx[8][3] = 1;
            }

            // ... than save some parameters that are equal in all dataSets
            dataSetParamStr.nx1 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX1());
            dataSetParamStr.nx2 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX2());
            dataSetParamStr.nx3 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX3());

            firstBlock = false;

            // how many elements are in all arrays of DataSet (equal in all blocks)
            int doubleCount = 0, temp;
            for (int i = 0; i<9; i++)   // 9 arrays ( averageValues, averageVelocity, averageFluktuations,
            {                 // averageTriplecorrelations, shearStressValues, relaxationFactor, 3 * fdistributions
               temp = 1;
               for (int ii = 0; ii < 4; ii++)
               {
                  temp *= dataSetParamStr.nx[i][ii];
                  //std::cout << ",dataSetParamStr.nx[" << i << "][" << ii << "]" << "=" << dataSetParamStr.nx[i][ii];
               }
               doubleCount += temp;
            }
            dataSetParamStr.doubleCountInBlock = doubleCount;
         }
         //std::cout << ",doubleCountInBlock="<<dataSetParamStr.doubleCountInBlock<< "," << dataSetParamStr.nx1 << "," << dataSetParamStr.nx2 << "," << dataSetParamStr.nx3 << std::endl;

         std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageValuesArray3DPtr = block->getKernel()->getDataSet()->getAverageDencity();
         if (AverageValuesArray3DPtr&&(dataSetParamStr.nx[0][0]>0)&&(dataSetParamStr.nx[0][1]>0)&&(dataSetParamStr.nx[0][2]>0)&&(dataSetParamStr.nx[0][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageValuesArray3DPtr->getDataVector().begin(), AverageValuesArray3DPtr->getDataVector().end());

         std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
         if (AverageVelocityArray3DPtr&&(dataSetParamStr.nx[1][0]>0)&&(dataSetParamStr.nx[1][1]>0)&&(dataSetParamStr.nx[1][2]>0)&&(dataSetParamStr.nx[1][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageVelocityArray3DPtr->getDataVector().begin(), AverageVelocityArray3DPtr->getDataVector().end());

         std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
         if (AverageFluctArray3DPtr&&(dataSetParamStr.nx[2][0]>0)&&(dataSetParamStr.nx[2][1]>0)&&(dataSetParamStr.nx[2][2]>0)&&(dataSetParamStr.nx[2][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageFluctArray3DPtr->getDataVector().begin(), AverageFluctArray3DPtr->getDataVector().end());

         std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
         if (AverageTripleArray3DPtr&&(dataSetParamStr.nx[3][0]>0)&&(dataSetParamStr.nx[3][1]>0)&&(dataSetParamStr.nx[3][2]>0)&&(dataSetParamStr.nx[3][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageTripleArray3DPtr->getDataVector().begin(), AverageTripleArray3DPtr->getDataVector().end());

         std::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
         if (ShearStressValArray3DPtr&&(dataSetParamStr.nx[4][0]>0)&&(dataSetParamStr.nx[4][1]>0)&&(dataSetParamStr.nx[4][2]>0)&&(dataSetParamStr.nx[4][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), ShearStressValArray3DPtr->getDataVector().begin(), ShearStressValArray3DPtr->getDataVector().end());

         std::shared_ptr< CbArray3D<LBMReal, IndexerX3X2X1> > RelaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
         if (RelaxationFactor3DPtr&&(dataSetParamStr.nx[5][0]>0)&&(dataSetParamStr.nx[5][1]>0)&&(dataSetParamStr.nx[5][2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), RelaxationFactor3DPtr->getDataVector().begin(), RelaxationFactor3DPtr->getDataVector().end());

         std::shared_ptr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
         CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
         if (localDistributions&&(dataSetParamStr.nx[6][0]>0)&&(dataSetParamStr.nx[6][1]>0)&&(dataSetParamStr.nx[6][2]>0)&&(dataSetParamStr.nx[6][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), localDistributions->getDataVector().begin(), localDistributions->getDataVector().end());

         CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
         if (nonLocalDistributions&&(dataSetParamStr.nx[7][0]>0)&&(dataSetParamStr.nx[7][1]>0)&&(dataSetParamStr.nx[7][2]>0)&&(dataSetParamStr.nx[7][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), nonLocalDistributions->getDataVector().begin(), nonLocalDistributions->getDataVector().end());

         CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();
         if (zeroDistributions&&(dataSetParamStr.nx[8][0]>0)&&(dataSetParamStr.nx[8][1]>0)&&(dataSetParamStr.nx[8][2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), zeroDistributions->getDataVector().begin(), zeroDistributions->getDataVector().end());

         ic++;
      }
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(dataSetParamStr.doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);
   mpiTypeFreeFlag = true;

   //doubleValuesArrayGW.assign(doubleValuesArray.begin(), doubleValuesArray.end());
   
   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeDataSet start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   // write to the file
   // all processes calculate their offsets (quantity of bytes that the process is going to write) 
   // and notify the next process (with the rank = rank + 1)
   MPI_Offset write_offset = (MPI_Offset)(size * sizeof(int));
   size_t next_write_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_write_offset = write_offset + sizeof(dataSetParam) + blocksCount * (sizeof(DataSet)+ dataSetParamStr.doubleCountInBlock*sizeof(double));
         MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_write_offset = write_offset + sizeof(dataSetParam) + blocksCount * (sizeof(DataSet)+ dataSetParamStr.doubleCountInBlock*sizeof(double));
         if (rank<size-1)
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
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
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE| MPI_MODE_WRONLY, info, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   //std::cout << "writeDataSet rank=" << rank << ",blocksCount=" << blocksCount;
   //std::cout << ", rank*sizeof(int)=" << (MPI_Offset)(rank * sizeof(int)) << ", write_offset=" << write_offset << std::endl;

   // each process writes the quantity of it's blocks
   MPI_File_write_at(file_handler, (MPI_Offset)(rank*sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // each process writes common parameters of a dataSet
   MPI_File_write_at(file_handler, write_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   // each process writes data identifying blocks
   MPI_File_write_at(file_handler, (MPI_Offset)(write_offset+sizeof(dataSetParam)), dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
   // each process writes the dataSet arrays
   MPI_File_write_at(file_handler, (MPI_Offset)(write_offset+sizeof(dataSetParam)+blocksCount*sizeof(DataSet)), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
   MPI_File_sync(file_handler);
   
   //int blockC;
   //MPI_File_read_at(file_handler, (MPI_Offset)(rank * sizeof(int)), &blockC, 1, MPI_INT, MPI_STATUS_IGNORE);
   //std::cout << "readDataSet rank=" << rank << ", blockC=" << blockC << std::endl;
   
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeDataSet time: "<<finish-start<<" s");
   }

   delete[] dataSetArray;
}

void MPIIORestart11CoProcessor::writeBoundaryConds(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBoundaryConds start collect data rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   int blocksCount = 0;          // quantity of blocks in the grid, max 2147483648 blocks!
   size_t count_boundCond = 0;	// how many BoundaryConditions in all blocks
   int count_indexContainer = 0;	// how many indexContainer-values in all blocks
   size_t byteCount = 0;			// how many bytes writes this process in the file 

   std::vector<Block3DPtr> blocksVector[25];
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      grid->getBlocks(level, rank, blocksVector[level]);
      blocksCount += static_cast<int>(blocksVector[level].size());
   }

   BCAdd* bcAddArray = new BCAdd[blocksCount];
   std::vector<BoundaryCondition> bcVector;
   std::vector<int> bcindexmatrixV;
   std::vector<int> indexContainerV;
   bool bcindexmatrixCountNotInit = true;

   int ic = 0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      for(Block3DPtr block : blocksVector[level])  // all the blocks of the current level
      {
         BCArray3DPtr bcArr = block->getKernel()->getBCProcessor()->getBCArray();

         bcAddArray[ic].x1 = block->getX1(); // coordinates of the block needed to find it while regenerating the grid
         bcAddArray[ic].x2 = block->getX2();
         bcAddArray[ic].x3 = block->getX3();
         bcAddArray[ic].level = block->getLevel();
         bcAddArray[ic].boundCond_count = 0; // how many BoundaryConditions in this block
         bcAddArray[ic].indexContainer_count = 0;  // how many indexContainer-values in this block

         for (int bc = 0; bc<bcArr->getBCVectorSize(); bc++)
         {
            BoundaryCondition* bouCond = new BoundaryCondition();
            if (bcArr->bcvector[bc]==NULL)
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
            //std::cout << "writeBoundaryConds noslipBoundaryFlags="<< bouCond->noslipBoundaryFlags << std::endl;
            bcVector.push_back(*bouCond);
            //bcVectorGW.push_back(*bouCond);
            //if (bcVector[count_boundCond].noslipBoundaryFlags != bcVectorGW[count_boundCond].noslipBoundaryFlags)
            //   std::cout << "bcVector[count_boundCond].noslipBoundaryFlags != bcVectorGW[count_boundCond].noslipBoundaryFlags!!!" << std::endl;
            bcAddArray[ic].boundCond_count++;
            count_boundCond++;
         }

         // the quantity of elements in the bcindexmatrix array (CbArray3D<int, IndexerX3X2X1>) in bcArray(BCArray3D) is always equal,
         // this will be the size of the "write-read-block" in MPI_write_.../MPI_read-functions when writing/reading BoundConds
         if (bcindexmatrixCountNotInit)
         {
            boundCondParamStr.nx1 = static_cast<int>(bcArr->bcindexmatrix.getNX1());
            boundCondParamStr.nx2 = static_cast<int>(bcArr->bcindexmatrix.getNX2());
            boundCondParamStr.nx3 = static_cast<int>(bcArr->bcindexmatrix.getNX3());
            boundCondParamStr.bcindexmatrixCount = static_cast<int>(bcArr->bcindexmatrix.getDataVector().size());
            bcindexmatrixCountNotInit = false;
         }
         bcindexmatrixV.insert(bcindexmatrixV.end(), bcArr->bcindexmatrix.getDataVector().begin(), bcArr->bcindexmatrix.getDataVector().end());

         indexContainerV.insert(indexContainerV.end(), bcArr->indexContainer.begin(), bcArr->indexContainer.end());
         bcAddArray[ic].indexContainer_count = static_cast<int>(bcArr->indexContainer.size());
         count_indexContainer += bcAddArray[ic].indexContainer_count;

         ic++;
      }
   }

   //bcindexmatrixVGW.assign(bcindexmatrixV.begin(), bcindexmatrixV.end());
   //indexContainerVGW.assign(indexContainerV.begin(), indexContainerV.end());
   
   MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
   MPI_Type_commit(&bcindexmatrixType);
   mpiTypeFreeFlag = true;

   //how many "big blocks" of BLOCK_SIZE size can by formed
   int bcBlockCount = (int)(count_boundCond/BLOCK_SIZE);
   if (bcBlockCount * BLOCK_SIZE<count_boundCond)
      bcBlockCount += 1;
   for (int i = (int)count_boundCond; i<bcBlockCount * BLOCK_SIZE; i++)
   {
      BoundaryCondition* bouCond = new BoundaryCondition();
      memset(bouCond, 0, sizeof(BoundaryCondition));
      bcVector.push_back(*bouCond);
   }

   byteCount = bcBlockCount * BLOCK_SIZE * sizeof(BoundaryCondition) + blocksCount * sizeof(BCAdd) + sizeof(int) * (blocksCount * boundCondParamStr.bcindexmatrixCount + count_indexContainer);

   // write to the file
   // all processes calculate their offsets (quantity of bytes that the process is going to write) 
   // and notify the next process (with the rank = rank + 1)
   MPI_Offset write_offset = (MPI_Offset)(size * (3 * sizeof(int) + sizeof(boundCondParam)));
   size_t next_write_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_write_offset = write_offset + byteCount;
         MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&write_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_write_offset = write_offset + byteCount;
         if (rank<size-1)
            MPI_Send(&next_write_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBoundaryConds start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
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
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpBC.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   MPI_Offset write_offset1 = (MPI_Offset)(rank * (3 * sizeof(int) + sizeof(boundCondParam)));

   // each process writes the quantity of it's blocks
   MPI_File_write_at(file_handler, write_offset1, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // each process writes the quantity of "big blocks" of BLOCK_SIZE of boundary conditions
   MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1+sizeof(int)), &bcBlockCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // each process writes the quantity of indexContainer elements in all blocks
   MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1+2*sizeof(int)), &count_indexContainer, 1, MPI_INT, MPI_STATUS_IGNORE);
   // each process writes the quantity of bcindexmatrix elements in every block
   MPI_File_write_at(file_handler, (MPI_Offset)(write_offset1+3*sizeof(int)), &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);

   //std::cout << "rank=" << rank << ",(rank*write_offset1)=" << rank*write_offset1<< ",blocksCount=" << blocksCount ;
   //std::cout << ", " << rank*write_offset1 + sizeof(int) << ",bcBlockCount=" << bcBlockCount;
   //std::cout << ", " << rank*write_offset1 + 2 * sizeof(int) << ",count_indexContainer=" << count_indexContainer;
   //std::cout << ", " << rank*write_offset1 + 3 * sizeof(int) << ",boundCondParamStr=" << boundCondParamStr.bcindexmatrixCount << std::endl;

   // each process writes data identifying the blocks
   MPI_File_write_at(file_handler, write_offset, bcAddArray, blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
   // each process writes boundary conditions
   if (bcVector.size()>0)
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset+blocksCount*sizeof(BCAdd)), &bcVector[0], bcBlockCount, boundCondType1000, MPI_STATUS_IGNORE);
   // each process writes bcindexmatrix values
   if (bcindexmatrixV.size()>0)
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset+blocksCount*sizeof(BCAdd)+bcBlockCount*BLOCK_SIZE*sizeof(BoundaryCondition)), &bcindexmatrixV[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
   // each process writes indexContainer values
   if (indexContainerV.size()>0)
      MPI_File_write_at(file_handler, (MPI_Offset)(write_offset+blocksCount*sizeof(BCAdd)+bcBlockCount*BLOCK_SIZE*sizeof(BoundaryCondition)+blocksCount*boundCondParamStr.bcindexmatrixCount*sizeof(int)), &indexContainerV[0], count_indexContainer, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_sync(file_handler);

   //std::cout <<"rank="<<rank<<",blocksCount="<< blocksCount<<", "<< bcBlockCount<<", "<< count_indexContainer<<", "<< bcindexmatrixCount << std::endl;
   //std::cout <<"rank="<<rank<<",write_offset="<< write_offset <<", "<< write_offset + blocksCount * sizeof(BCAdd) <<", "<< write_offset + blocksCount * sizeof(BCAdd) + bcBlockCount*BLOCK_SIZE * sizeof(BoundaryCondition) <<", "<< write_offset + blocksCount * sizeof(BCAdd) + bcBlockCount*BLOCK_SIZE * sizeof(BoundaryCondition) + blocksCount*bcindexmatrixCount * sizeof(int)<< std::endl;

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::writeBoundaryConds time: "<<finish-start<<" s");
   }

   delete[] bcAddArray;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIORestart11CoProcessor::restart(int step)
{
   if (comm->isRoot()) UBLOG(logINFO, "MPIIORestart11CoProcessor restart step: "<<step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - start");
   readBlocks(step);
   readDataSet(step);
   readBoundaryConds(step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - end");
   this->reconnect(grid);
}

void MPIIORestart11CoProcessor::readBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBlocks start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpBlocks.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   // read count of blocks
   int blocksCount = 0;
   //MPI_File_read_at(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   Block3d* block3dArray = new Block3d[blocksCount];

   // calculate the read offset
   MPI_Offset read_offset = (MPI_Offset)(size * sizeof(int));

   GridParam* gridParameters = new GridParam;

   // read parameters of the grid
   MPI_File_read_at(file_handler, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   // read all the blocks
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBlocks time: "<<finish-start<<" s");
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBlocks start of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   // clear the grid
   std::vector<Block3DPtr> blocksVector;
   grid->getBlocks(0, blocksVector);
   int del = 0;
   for(Block3DPtr block : blocksVector)
   {
      grid->deleteBlock(block);
      del++;
   }

   // restore the grid
   CoordinateTransformation3DPtr trafo(new CoordinateTransformation3D());
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
      Block3DPtr block(new Block3D(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, block3dArray[n].level));
      block->setActive(block3dArray[n].active);
      block->setBundle(block3dArray[n].bundle);
      block->setRank(block3dArray[n].rank);
      block->setLocalRank(block3dArray[n].lrank);
      block->setGlobalID(block3dArray[n].globalID);
      block->setLocalID(block3dArray[n].localID);
      block->setPart(block3dArray[n].part);
      block->setLevel(block3dArray[n].level);
      block->interpolationFlagCF = block3dArray[n].interpolationFlagCF;
      block->interpolationFlagFC = block3dArray[n].interpolationFlagFC;

      grid->addBlock(block);
   }

   delete gridParameters;
   delete[] block3dArray;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBlocks end of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
}

void MPIIORestart11CoProcessor::readDataSet(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readDataSet start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   // read count of blocks
   int blocksCount = 0;
   MPI_File_read_at(file_handler, (MPI_Offset)(rank*sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(size*sizeof(int)), &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);
   //std::cout <<"MPIIORestart11CoProcessor::readDataSet rank=" << rank <<", dataSetParamStr.doubleCountInBlock="<< dataSetParamStr.doubleCountInBlock << std::endl;
   //std::cout << ",dataSetParamStr.nx[6][0]" << "=" << dataSetParamStr.nx[6][0] << "," << dataSetParamStr.nx[6][1] << "," << dataSetParamStr.nx[6][2] << "," << dataSetParamStr.nx[6][3];
   //std::cout << ",doubleCountInBlock=" << dataSetParamStr.doubleCountInBlock << "," << dataSetParamStr.nx1 << "," << dataSetParamStr.nx2 << "," << dataSetParamStr.nx3 << std::endl;

   DataSet* dataSetArray = new DataSet[blocksCount];
   std::vector<double> doubleValuesArray(blocksCount * dataSetParamStr.doubleCountInBlock); // double-values in all blocks 
   
   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(dataSetParamStr.doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);
   mpiTypeFreeFlag = true;
   //std::cout << "MPIIORestart11CoProcessor::readDataSet rank=" << rank << " 123=" << dataSetParamStr.doubleCountInBlock << std::endl;

   // calculate the read offset
   MPI_Offset read_offset = (MPI_Offset)(size * sizeof(int));
   size_t next_read_offset = 0;

   if(size > 1)
   {
   	if(rank == 0)
   	{
   		next_read_offset = read_offset + sizeof(dataSetParam) + blocksCount * (sizeof(DataSet) + dataSetParamStr.doubleCountInBlock * sizeof(double));
   		MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
   	}
   	else
   	{
   		MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_read_offset = read_offset + sizeof(dataSetParam) + blocksCount * (sizeof(DataSet) + dataSetParamStr.doubleCountInBlock * sizeof(double));
   		if(rank < size - 1)
   			MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
   	}
   }

   /*int chunkFlag = 0;

   if (rank == 0)
   {
      MPI_File_read_at(file_handler, read_offset, dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+blocksCount*sizeof(DataSet)), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
       
      for (int i=1; i<size; i+=chunk)
      {
         for (int j=i; j<i+chunk; j++)
         {
            if (j < size)
            {
               MPI_Send(&chunkFlag, 1, MPI_INT, j, 77, MPI_COMM_WORLD);
               //UBLOG(logINFO, "j= "<<j);
            }
         }
         for (int j=i; j<i+chunk; j++)
         {
            if (j < size)
            {
               MPI_Recv(&chunkFlag, 1, MPI_INT, j, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
         }
      }
   }
   else
   {
      MPI_Recv(&chunkFlag, 1, MPI_INT, 0, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handler, read_offset, dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+blocksCount*sizeof(DataSet)), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
      MPI_Send(&chunkFlag, 1, MPI_INT, 0, 77, MPI_COMM_WORLD);
      //UBLOG(logINFO, "read rank= "<<rank);
   }*/

   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+sizeof(dataSetParam)), dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+sizeof(dataSetParam)+blocksCount*sizeof(DataSet)), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
   MPI_File_close(&file_handler);

   //for (int ch = 0; ch < blocksCount; ch++)
   //{
   //   if ((dataSetArrayGW[ch].x1 != dataSetArray[ch].x1) ||
   //      (dataSetArrayGW[ch].x2 != dataSetArray[ch].x2) ||
   //      (dataSetArrayGW[ch].x3 != dataSetArray[ch].x3) ||
   //      (dataSetArrayGW[ch].level != dataSetArray[ch].level) ||
   //      (dataSetArrayGW[ch].ghostLayerWidth != dataSetArray[ch].ghostLayerWidth) ||
   //      (dataSetArrayGW[ch].collFactor != dataSetArray[ch].collFactor) ||
   //      (dataSetArrayGW[ch].deltaT != dataSetArray[ch].deltaT) ||
   //      (dataSetArrayGW[ch].compressible != dataSetArray[ch].compressible) ||
   //      (dataSetArrayGW[ch].withForcing != dataSetArray[ch].withForcing)) 
   //      std::cout << "dataSetArrayGW != rank" << rank << ", !!!!!====="<< std::endl;
   //}
   //std::cout << "doubleValuesArrayGW.size" << doubleValuesArrayGW.size() << ", " << doubleValuesArray.size() << std::endl;
   //for (int vl = 0; vl < doubleValuesArrayGW.size(); vl++)
   //   if(doubleValuesArrayGW[vl] != doubleValuesArray[vl])
   //      std::cout << "doubleValuesArrayGW != rank" << rank << ", !!!!!====="<< std::endl;

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readDataSet time: "<<finish-start<<" s");
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readDataSet start of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   size_t index = 0, nextVectorSize = 0;
   std::vector<double> vectorsOfValues[9];
   for (int n = 0; n<blocksCount; n++)
   {
      for (int b = 0; b<9; b++) // assign approciate vectors for 9 dataSet arrays
      {
         nextVectorSize = dataSetParamStr.nx[b][0]* dataSetParamStr.nx[b][1]* dataSetParamStr.nx[b][2]* dataSetParamStr.nx[b][3];
         vectorsOfValues[b].assign(doubleValuesArray.data()+index, doubleValuesArray.data()+index+nextVectorSize);
         index += nextVectorSize;
      }

      // fill dataSet arrays
      AverageValuesArray3DPtr mAverageDensity;
      if ((dataSetParamStr.nx[0][0]==0)&&(dataSetParamStr.nx[0][1]==0)&&(dataSetParamStr.nx[0][2]==0)&&(dataSetParamStr.nx[0][3]==0))
         mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageDensity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[0], dataSetParamStr.nx[0][0], dataSetParamStr.nx[0][1], dataSetParamStr.nx[0][2], dataSetParamStr.nx[0][3]));

      AverageValuesArray3DPtr mAverageVelocity;
      if ((dataSetParamStr.nx[1][0]==0)&&(dataSetParamStr.nx[1][1]==0)&&(dataSetParamStr.nx[1][2]==0)&&(dataSetParamStr.nx[1][3]==0))
         mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[1], dataSetParamStr.nx[1][0], dataSetParamStr.nx[1][1], dataSetParamStr.nx[1][2], dataSetParamStr.nx[1][3]));

      AverageValuesArray3DPtr mAverageFluktuations;
      if ((dataSetParamStr.nx[2][0]==0)&&(dataSetParamStr.nx[2][1]==0)&&(dataSetParamStr.nx[2][2]==0)&&(dataSetParamStr.nx[2][3]==0))
         mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[2], dataSetParamStr.nx[2][0], dataSetParamStr.nx[2][1], dataSetParamStr.nx[2][2], dataSetParamStr.nx[2][3]));

      AverageValuesArray3DPtr mAverageTriplecorrelations;
      if ((dataSetParamStr.nx[3][0]==0)&&(dataSetParamStr.nx[3][1]==0)&&(dataSetParamStr.nx[3][2]==0)&&(dataSetParamStr.nx[3][3]==0))
         mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[3], dataSetParamStr.nx[3][0], dataSetParamStr.nx[3][1], dataSetParamStr.nx[3][2], dataSetParamStr.nx[3][3]));

      ShearStressValuesArray3DPtr mShearStressValues;
      if ((dataSetParamStr.nx[4][0]==0)&&(dataSetParamStr.nx[4][1]==0)&&(dataSetParamStr.nx[4][2]==0)&&(dataSetParamStr.nx[4][3]==0))
         mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[4], dataSetParamStr.nx[4][0], dataSetParamStr.nx[4][1], dataSetParamStr.nx[4][2], dataSetParamStr.nx[4][3]));

      RelaxationFactorArray3DPtr mRelaxationFactor;
      if ((dataSetParamStr.nx[5][0]==0)&&(dataSetParamStr.nx[5][1]==0)&&(dataSetParamStr.nx[5][2]==0))
         mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr();
      else
         mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues[5], dataSetParamStr.nx[5][0], dataSetParamStr.nx[5][1], dataSetParamStr.nx[5][2]));

      //DistributionArray3DPtr mFdistributions(new D3Q27EsoTwist3DSplittedVector(dataSetParamStr.nx1, dataSetParamStr.nx2, dataSetParamStr.nx3, -999.0));
      DistributionArray3DPtr mFdistributions(new D3Q27EsoTwist3DSplittedVector());

      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[6], dataSetParamStr.nx[6][0], dataSetParamStr.nx[6][1], dataSetParamStr.nx[6][2], dataSetParamStr.nx[6][3])));
      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[7], dataSetParamStr.nx[7][0], dataSetParamStr.nx[7][1], dataSetParamStr.nx[7][2], dataSetParamStr.nx[7][3])));
      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues[8], dataSetParamStr.nx[8][0], dataSetParamStr.nx[8][1], dataSetParamStr.nx[8][2])));

      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(dataSetParamStr.nx1);
      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(dataSetParamStr.nx2);
      std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(dataSetParamStr.nx3);

      DataSet3DPtr dataSetPtr = DataSet3DPtr(new DataSet3D());
      dataSetPtr->setAverageDencity(mAverageDensity);
      dataSetPtr->setAverageVelocity(mAverageVelocity);
      dataSetPtr->setAverageFluctuations(mAverageFluktuations);
      dataSetPtr->setAverageTriplecorrelations(mAverageTriplecorrelations);
      dataSetPtr->setShearStressValues(mShearStressValues);
      dataSetPtr->setRelaxationFactor(mRelaxationFactor);
      dataSetPtr->setFdistributions(mFdistributions);

      // find the nesessary block and fill it
      Block3DPtr block = grid->getBlock(dataSetArray[n].x1, dataSetArray[n].x2, dataSetArray[n].x3, dataSetArray[n].level);
      //LBMKernelPtr kernel(new CompressibleCumulantLBMKernel());
      //LBMKernelPtr kernel(new IncompressibleCumulantLBMKernel());
      LBMKernelPtr kernel = this->lbmKernel->clone();
      kernel->setGhostLayerWidth(dataSetArray[n].ghostLayerWidth);
      kernel->setCollisionFactor(dataSetArray[n].collFactor);
      kernel->setDeltaT(dataSetArray[n].deltaT);
      kernel->setCompressible(dataSetArray[n].compressible);
      kernel->setWithForcing(dataSetArray[n].withForcing);
      kernel->setDataSet(dataSetPtr);
      block->setKernel(kernel);
      //block->getKernel()->setDataSet(dataSetPtr);
   }

   delete[] dataSetArray;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet end of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
}

void MPIIORestart11CoProcessor::readBoundaryConds(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBoundaryConds start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handler;
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpBC.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   int blocksCount = 0;
   int dataCount1000 = 0;
   int dataCount2 = 0;
   MPI_Offset read_offset1 = (MPI_Offset)(rank * (3 * sizeof(int) + sizeof(boundCondParam)));

   // read count of blocks
   MPI_File_read_at(file_handler, read_offset1, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // read count of big BoundaryCondition blocks
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1+sizeof(int)), &dataCount1000, 1, MPI_INT, MPI_STATUS_IGNORE);
   // read count of indexContainer values in all blocks
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1+2*sizeof(int)), &dataCount2, 1, MPI_INT, MPI_STATUS_IGNORE);
   // read count of bcindexmatrix values in every block
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset1+3*sizeof(int)), &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);

   //std::cout << "rank=" << rank << ",(rank*read_offset1)=" << rank*read_offset1 << ",blocksCount=" << blocksCount;
   //std::cout << ", " << rank*read_offset1 + sizeof(int) << ",bcBlockCount=" << dataCount1000;
   //std::cout << ", " << rank*read_offset1 + 2 * sizeof(int) << ",count_indexContainer=" << dataCount2;
   //std::cout << ", " << rank*read_offset1 + 3 * sizeof(int) << ",boundCondParamStr=" << boundCondParamStr.bcindexmatrixCount << std::endl;
   //std::cout << "readrank=" << rank << ",blocksCount=" << blocksCount << ", " << dataCount1000 << ", " << dataCount2 << ", " << boundCondParamStr.bcindexmatrixCount << std::endl;

   MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
   MPI_Type_commit(&bcindexmatrixType);
   mpiTypeFreeFlag = true;

   size_t dataCount = dataCount1000 * BLOCK_SIZE;
   BCAdd* bcAddArray = new BCAdd[blocksCount];
   BoundaryCondition* bcArray = new BoundaryCondition[dataCount];
   BoundaryCondition* nullBouCond = new BoundaryCondition();
   memset(nullBouCond, 0, sizeof(BoundaryCondition));
   int* intArray1 = new int[blocksCount * boundCondParamStr.bcindexmatrixCount];
   int* intArray2 = new int[dataCount2];

   MPI_Offset read_offset = (MPI_Offset)(size * (3 * sizeof(int) + sizeof(boundCondParam)));
   size_t next_read_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_read_offset = read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+(blocksCount * boundCondParamStr.bcindexmatrixCount + dataCount2)*sizeof(int);
         MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_read_offset = read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+(blocksCount * boundCondParamStr.bcindexmatrixCount + dataCount2)*sizeof(int);
         if (rank<size-1)
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
   }
   //std::cout << "readrank=" << rank << ",read_offset=" << read_offset << ", " << read_offset + blocksCount * sizeof(BCAdd) << ", " << read_offset + blocksCount * sizeof(BCAdd) + dataCount * sizeof(BoundaryCondition) << ", " << read_offset + blocksCount * sizeof(BCAdd) + dataCount * sizeof(BoundaryCondition) + blocksCount * bcindexmatrixCount * sizeof(int) << std::endl;

   MPI_File_read_at(file_handler, read_offset, bcAddArray, blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+blocksCount*sizeof(BCAdd)), &bcArray[0], dataCount1000, boundCondType1000, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)), &intArray1[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handler, (MPI_Offset)(read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+blocksCount * boundCondParamStr.bcindexmatrixCount*sizeof(int)), &intArray2[0], dataCount2, MPI_INT, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBoundaryConds time: "<<finish-start<<" s");
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBoundaryConds start of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   int index = 0, index1 = 0, index2 = 0;
   std::vector<BoundaryConditionsPtr> bcVector;
   std::vector<int> bcindexmatrixV;
   std::vector<int> indexContainerV;

   for (size_t n = 0; n<blocksCount; n++)
   {
      bcVector.resize(0);
      bcindexmatrixV.resize(0);
      indexContainerV.resize(0);

      for (size_t ibc = 0; ibc<bcAddArray[n].boundCond_count; ibc++)
      {
         BoundaryConditionsPtr bc;
         if (memcmp(&bcArray[index], nullBouCond, sizeof(BoundaryCondition))==0)
            bc = BoundaryConditionsPtr();
         else
         {
            bc = BoundaryConditionsPtr(new BoundaryConditions);
            bc->noslipBoundaryFlags = bcArray[index].noslipBoundaryFlags;
            bc->slipBoundaryFlags = bcArray[index].slipBoundaryFlags;
            bc->densityBoundaryFlags = bcArray[index].densityBoundaryFlags;
            bc->velocityBoundaryFlags = bcArray[index].velocityBoundaryFlags;
            bc->wallModelBoundaryFlags = bcArray[index].wallModelBoundaryFlags;
            bc->bcVelocityX1 = bcArray[index].bcVelocityX1;
            bc->bcVelocityX2 = bcArray[index].bcVelocityX2;
            bc->bcVelocityX3 = bcArray[index].bcVelocityX3;
            bc->bcDensity = bcArray[index].bcDensity;
            bc->bcLodiDensity = bcArray[index].bcLodiDensity;
            bc->bcLodiVelocityX1 = bcArray[index].bcLodiVelocityX1;
            bc->bcLodiVelocityX2 = bcArray[index].bcLodiVelocityX2;
            bc->bcLodiVelocityX3 = bcArray[index].bcLodiVelocityX3;
            bc->bcLodiLentgh = bcArray[index].bcLodiLentgh;

            bc->nx1 = bcArray[index].nx1;
            bc->nx2 = bcArray[index].nx2;
            bc->nx3 = bcArray[index].nx3;
            for (int iq = 0; iq<26; iq++)
               bc->setQ(bcArray[index].q[iq], iq);
            bc->setBcAlgorithmType(bcArray[index].algorithmType);

            //if (bcVectorGW[index].nx1 != bc->nx1)
            //   std::cout << "readBoundaryConds nx1 !!!!===" << bcVectorGW[index].nx1 << " ---- " << bc->nx1 << std::endl;
            //if (bcVectorGW[index].nx2 != bc->nx2)
            //   std::cout << "readBoundaryConds nx2 !!!!===" << bcVectorGW[index].nx2 << " ---- " << bc->nx2 << std::endl;
            //if (bcVectorGW[index].nx3 != bc->nx3)
            //   std::cout << "readBoundaryConds nx3 !!!!===" << bcVectorGW[index].nx3 << " ---- " << bc->nx3 << std::endl;
            //if (bcVectorGW[index].algorithmType != bc->algorithmType)
            //   std::cout << "readBoundaryConds algorithmType !!!!===" << bcVectorGW[index].algorithmType << " ---- " << bc->algorithmType << std::endl;
            //for (int iq = 0; iq<26; iq++)
            //   if (bcVectorGW[index].q[iq] != bc->q[iq])
            //   std::cout << "readBoundaryConds q !!!!===" /*<< bcVectorGW[index].q << " ---- " << bc->q*/ << std::endl;
            //std::cout << "readBoundaryConds BoundaryConditionsPtr !!!!===" <<std::endl;

         }

         bcVector.push_back(bc);
         index++;
      }

      for (int b1 = 0; b1 < boundCondParamStr.bcindexmatrixCount; b1++)
      {
         //if (bcindexmatrixVGW[index1] != intArray1[index1])
         //   std::cout << "readBoundaryConds bcindexmatrixVGW !!!!===" << std::endl;
         bcindexmatrixV.push_back(intArray1[index1++]);
      }
      for (int b2 = 0; b2 < bcAddArray[n].indexContainer_count; b2++)
      {
         //if (indexContainerVGW[index2] != intArray2[index2])
         //   std::cout << "readBoundaryConds indexContainerVGW !!!!===" << std::endl;
         indexContainerV.push_back(intArray2[index2++]);
      }

      CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, boundCondParamStr.nx1, boundCondParamStr.nx2, boundCondParamStr.nx3);

      Block3DPtr block = grid->getBlock(bcAddArray[n].x1, bcAddArray[n].x2, bcAddArray[n].x3, bcAddArray[n].level);
      //if(!block) std::cout << "readBoundaryConds can't find the block!!!" << std::endl;
      BCProcessorPtr bcProc = bcProcessor->clone(block->getKernel());
      //if(!bcProc) std::cout << "readBoundaryConds can't find the bcProc!!!" << std::endl;
      BCArray3DPtr bcArr(new BCArray3D());
      bcArr->bcindexmatrix = bcim;
      bcArr->bcvector = bcVector;
      bcArr->indexContainer = indexContainerV;
      bcProc->setBCArray(bcArr);
      
      //if (!(block->getKernel())) 
      //   std::cout << "readBoundaryConds kernel=" << block->getKernel() <<" "<< bcAddArray[n].x1 << " " << bcAddArray[n].x2 << " " << bcAddArray[n].x3 << std::endl;
      block->getKernel()->setBCProcessor(bcProc);
   }

   delete nullBouCond;
   delete[] bcArray;
   delete[] bcAddArray;
   delete[] intArray1;
   delete[] intArray2;
   
   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestart11CoProcessor::readBoundaryConds end of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::setChunk(int val)
{
   chunk = val;
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::setLBMKernel(LBMKernelPtr kernel)
{
   this->lbmKernel = kernel;
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestart11CoProcessor::setBCProcessor(BCProcessorPtr bcProcessor)
{
   this->bcProcessor = bcProcessor;
}

