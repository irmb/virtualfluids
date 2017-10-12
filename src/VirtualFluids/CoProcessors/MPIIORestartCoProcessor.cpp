#include "MPIIORestartCoProcessor.h"
#include <boost/foreach.hpp>
#include "D3Q27System.h"
//#include "LBMKernel.h"
#include "CompressibleCumulantLBMKernel.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <UbSystem.h>
#include <MemoryUtil.h>

//! BLOCK_SIZE defines the quantity of the BoundaryCondition-structures written as one block to the file
//! To avoid overflow in the parameter \a count of the function MPI_File_write_at 
//! structures BoundaryCondition are being written in blocks containing each of them BLOCK_SIZE structures
#define BLOCK_SIZE 1024

MPIIORestartCoProcessor::MPIIORestartCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
   const std::string& path,
   CommunicatorPtr comm) :
   CoProcessor(grid, s),
   path(path),
   comm(comm),
   mpiTypeFreeFlag(false)
{
   UbSystem::makeDirectory(path+"/mpi_io_cp");

   memset(&blockParamStr, 0, sizeof(blockParamStr));

   //-------------------------   define MPI types  ---------------------------------

   MPI_Datatype typesGP[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksGP[3] = { 34, 6, 5 };
   MPI_Aint offsetsGP[3], lbGP, extentGP;

   offsetsGP[0] = 0;
   MPI_Type_get_extent(MPI_DOUBLE, &lbGP, &extentGP);
   offsetsGP[1] = blocksGP[0]*extentGP;

   MPI_Type_get_extent(MPI_INT, &lbGP, &extentGP);
   offsetsGP[2] = offsetsGP[1]+blocksGP[1]*extentGP;

   MPI_Type_struct(3, blocksGP, offsetsGP, typesGP, &gridParamType);
   MPI_Type_commit(&gridParamType);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(41, MPI_INT, &blockParamType);
   MPI_Type_commit(&blockParamType);

   //-----------------------------------------------------------------------

   //MPI_Datatype types[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   //int blocks[3] = { 2, 14, 3 };
   //MPI_Aint offsets[3], lb, extent;

   //offsets[0] = 0;
   //MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
   //offsets[1] = blocks[0] * extent;

   //MPI_Type_get_extent(MPI_INT, &lb, &extent);
   //offsets[2] = offsets[1] + blocks[1] * extent;

   //MPI_Type_create_struct(3, blocks, offsets, types, &block3dType);
   //MPI_Type_commit(&block3dType);

   MPI_Datatype typesBlock[2] = { MPI_INT, MPI_CHAR };
   int blocksBlock[2] = { 13, 1 };
   MPI_Aint offsetsBlock[2], lbBlock, extentBlock;

   offsetsBlock[0] = 0;
   MPI_Type_get_extent(MPI_INT, &lbBlock, &extentBlock);
   offsetsBlock[1] = blocksBlock[0]*extentBlock;

   MPI_Type_create_struct(2, blocksBlock, offsetsBlock, typesBlock, &block3dType);
   MPI_Type_commit(&block3dType);

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

   //MPI_Type_contiguous(4, MPI_INT, &dataSetType);
   MPI_Type_commit(&dataSetType);

   //-----------------------------------------------------------------------

   MPI_Datatype typesBC[3] = { MPI_LONG_LONG_INT, MPI_FLOAT, MPI_CHAR };
   int blocksBC[3] = { 5, 38, 1 };
   MPI_Aint offsetsBC[3], lbBC, extentBC;

   offsetsBC[0] = 0;
   MPI_Type_get_extent(MPI_LONG_LONG_INT, &lbBC, &extentBC);
   offsetsBC[1] = blocksBC[0]*extentBC;

   MPI_Type_get_extent(MPI_FLOAT, &lbBC, &extentBC);
   offsetsBC[2] = blocksBC[1]*extentBC;

   //MPI_Type_struct(3, blocksBC, offsetsBC, typesBC, &boundCondType);
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
MPIIORestartCoProcessor::~MPIIORestartCoProcessor()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&blockParamType);
   MPI_Type_free(&block3dType);
   MPI_Type_free(&dataSetType);
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
void MPIIORestartCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      if (comm->isRoot()) UBLOG(logINFO, "MPIIORestartCoProcessor save step: "<<step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - start");
      writeBlocks((int)step);
      writeDataSet((int)step);
      writeBoundaryConds((int)step);
      if (comm->isRoot()) UBLOG(logINFO, "Save check point - end");
   }
}
//////////////////////////////////////////////////////////////////////////
void MPIIORestartCoProcessor::writeBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size=1;
   MPI_File file_handler;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBlocks start collect data rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

  // if (comm->isRoot())
  // {
      UbSystem::makeDirectory(path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step));
      std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpBlocks.bin";
      int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler);
      if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);
  // }

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
   bool firstBlock = true;
   int ic = 0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocksVector[level])  //	all the blocks of the current level
      {
         if (firstBlock && block->getKernel()) // when first (any) valid block...
         {
            boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageValuesArray3DPtr = block->getKernel()->getDataSet()->getAverageValues();
            if (AverageValuesArray3DPtr)
            {
               blockParamStr.nx[0][0] = static_cast<int>(AverageValuesArray3DPtr->getNX1());
               blockParamStr.nx[0][1] = static_cast<int>(AverageValuesArray3DPtr->getNX2());
               blockParamStr.nx[0][2] = static_cast<int>(AverageValuesArray3DPtr->getNX3());
               blockParamStr.nx[0][3] = static_cast<int>(AverageValuesArray3DPtr->getNX4());
            }

            boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
            if (AverageVelocityArray3DPtr)
            {
               blockParamStr.nx[1][0] = static_cast<int>(AverageVelocityArray3DPtr->getNX1());
               blockParamStr.nx[1][1] = static_cast<int>(AverageVelocityArray3DPtr->getNX2());
               blockParamStr.nx[1][2] = static_cast<int>(AverageVelocityArray3DPtr->getNX3());
               blockParamStr.nx[1][3] = static_cast<int>(AverageVelocityArray3DPtr->getNX4());
            }

            boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
            if (AverageFluctArray3DPtr)
            {
               blockParamStr.nx[2][0] = static_cast<int>(AverageFluctArray3DPtr->getNX1());
               blockParamStr.nx[2][1] = static_cast<int>(AverageFluctArray3DPtr->getNX2());
               blockParamStr.nx[2][2] = static_cast<int>(AverageFluctArray3DPtr->getNX3());
               blockParamStr.nx[2][3] = static_cast<int>(AverageFluctArray3DPtr->getNX4());
            }

            boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
            if (AverageTripleArray3DPtr)
            {
               blockParamStr.nx[3][0] = static_cast<int>(AverageTripleArray3DPtr->getNX1());
               blockParamStr.nx[3][1] = static_cast<int>(AverageTripleArray3DPtr->getNX2());
               blockParamStr.nx[3][2] = static_cast<int>(AverageTripleArray3DPtr->getNX3());
               blockParamStr.nx[3][3] = static_cast<int>(AverageTripleArray3DPtr->getNX4());
            }

            boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
            if (ShearStressValArray3DPtr)
            {
               blockParamStr.nx[4][0] = static_cast<int>(ShearStressValArray3DPtr->getNX1());
               blockParamStr.nx[4][1] = static_cast<int>(ShearStressValArray3DPtr->getNX2());
               blockParamStr.nx[4][2] = static_cast<int>(ShearStressValArray3DPtr->getNX3());
               blockParamStr.nx[4][3] = static_cast<int>(ShearStressValArray3DPtr->getNX4());
            }

            boost::shared_ptr< CbArray3D<LBMReal, IndexerX3X2X1> > relaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
            if (relaxationFactor3DPtr)
            {
               blockParamStr.nx[5][0] = static_cast<int>(relaxationFactor3DPtr->getNX1());
               blockParamStr.nx[5][1] = static_cast<int>(relaxationFactor3DPtr->getNX2());
               blockParamStr.nx[5][2] = static_cast<int>(relaxationFactor3DPtr->getNX3());
               blockParamStr.nx[5][3] = 1;
            }

            boost::shared_ptr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
            CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
            if (localDistributions)
            {
               blockParamStr.nx[6][0] = static_cast<int>(localDistributions->getNX1());
               blockParamStr.nx[6][1] = static_cast<int>(localDistributions->getNX2());
               blockParamStr.nx[6][2] = static_cast<int>(localDistributions->getNX3());
               blockParamStr.nx[6][3] = static_cast<int>(localDistributions->getNX4());
            }

            CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
            if (nonLocalDistributions)
            {
               blockParamStr.nx[7][0] = static_cast<int>(nonLocalDistributions->getNX1());
               blockParamStr.nx[7][1] = static_cast<int>(nonLocalDistributions->getNX2());
               blockParamStr.nx[7][2] = static_cast<int>(nonLocalDistributions->getNX3());
               blockParamStr.nx[7][3] = static_cast<int>(nonLocalDistributions->getNX4());
            }

            CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();
            if (zeroDistributions)
            {
               blockParamStr.nx[8][0] = static_cast<int>(zeroDistributions->getNX1());
               blockParamStr.nx[8][1] = static_cast<int>(zeroDistributions->getNX2());
               blockParamStr.nx[8][2] = static_cast<int>(zeroDistributions->getNX3());
               blockParamStr.nx[8][3] = 1;
            }

            // ... than save some parameters that are equal in all blocks
            blockParamStr.nx1 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX1());
            blockParamStr.nx2 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX2());
            blockParamStr.nx3 = static_cast<int>(block->getKernel()->getDataSet()->getFdistributions()->getNX3());

            firstBlock = false;

            // how many elements are in all arrays of DataSet (equal in all blocks)
            int doubleCount = 0, temp;
            for (int i = 0; i<9; i++)   // 9 arrays ( averageValues, averageVelocity, averageFluktuations,
            {                 // averageTriplecorrelations, shearStressValues, relaxationFactor, 3 * fdistributions
               temp = 1;
               for (int ii = 0; ii<4; ii++)
                  temp *= blockParamStr.nx[i][ii];
               doubleCount += temp;
            }
            blockParamStr.doubleCountInBlock = doubleCount;

            // the quantity of elements in the bcindexmatrix array (CbArray3D<int, IndexerX3X2X1>) in bcArray(BCArray3D) is always equal,
            // this will be the size of the "write-read-block" in MPI_write_.../MPI_read-functions when writing/reading BoundConds
            BCArray3D bcArr = block->getKernel()->getBCProcessor()->getBCArray();
            blockParamStr.bcindexmatrix_count = static_cast<int>(bcArr.bcindexmatrix.getDataVector().size());
         }

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
         //if(block->getKernel())
         //{
         //	block3dArray[ic].ghostLayerWidth = block->getKernel()->getGhostLayerWidth();
         //	block3dArray[ic].collFactor = block->getKernel()->getCollisionFactor();
   //         block3dArray[ic].deltaT = block->getKernel()->getDeltaT();
         //	block3dArray[ic].compressible = block->getKernel()->getCompressible();
         //	block3dArray[ic].withForcing = block->getKernel()->getWithForcing();
         //}
         //else 
         //{
         //	block3dArray[ic].ghostLayerWidth = 0;
         //	block3dArray[ic].collFactor = 0.0;
         //	block3dArray[ic].deltaT = 0.0;
         //	block3dArray[ic].compressible = false;
         //	block3dArray[ic].withForcing = false;
         //}

         ic++;
      }
   }

   // write to the file
   // all processes calculate their offsets (quantity of bytes that the process is going to write) 
   // and notify the next process (with the rank = rank + 1)
   size_t view_offset = size*sizeof(int);
   size_t next_view_offset = 0;

   //if (size>1)
   //{
   //   if (rank==0)
   //   {
   //      next_view_offset = view_offset+sizeof(GridParam)+sizeof(BlockParam)+blocksCount*sizeof(Block3d);
   //      MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
   //   }
   //   else
   //   {
   //      MPI_Recv(&view_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   //      next_view_offset = view_offset+sizeof(GridParam)+sizeof(BlockParam)+blocksCount*sizeof(Block3d);
   //      if (rank<size-1)
   //         MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
   //   }
   //}

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBlocks start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   if (comm->isRoot())
   {
      // each process writes the quantity of it's blocks
      MPI_File_write_at(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      // each process writes parameters of the grid
      MPI_File_write_at(file_handler, view_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
      // each process writes common parameters of a block
      MPI_File_write_at(file_handler, view_offset+sizeof(GridParam), &blockParamStr, 1, blockParamType, MPI_STATUS_IGNORE);
      // each process writes it's blocks
      MPI_File_write_at(file_handler, view_offset+sizeof(GridParam)+sizeof(BlockParam), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);
      //MPI_File_sync(file_handler);
   }
   MPI_File_close(&file_handler);
 
   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBlocks time: "<<finish-start<<" s");
   }

   // register new MPI-types depending on the block-specific information
   MPI_Type_contiguous(blockParamStr.doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   MPI_Type_contiguous(blockParamStr.bcindexmatrix_count, MPI_INT, &bcindexmatrixType);
   MPI_Type_commit(&bcindexmatrixType);

   mpiTypeFreeFlag = true;

   delete[] block3dArray;
   delete gridParameters;
}

void MPIIORestartCoProcessor::writeDataSet(int step)
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
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet start collect data rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   int ic = 0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocksVector[level])  //	blocks of the current level
      {
         dataSetArray[ic].x1 = block->getX1();     // coordinates of the block needed to find it while regenerating the grid
         dataSetArray[ic].x2 = block->getX2();
         dataSetArray[ic].x3 = block->getX3();
         dataSetArray[ic].level = block->getLevel();
         if (block->getKernel())
         {
            dataSetArray[ic].ghostLayerWidth = block->getKernel()->getGhostLayerWidth();
            dataSetArray[ic].collFactor = block->getKernel()->getCollisionFactor();
            dataSetArray[ic].deltaT = block->getKernel()->getDeltaT();
            dataSetArray[ic].compressible = block->getKernel()->getCompressible();
            dataSetArray[ic].withForcing = block->getKernel()->getWithForcing();
         }
         else
         {
            dataSetArray[ic].ghostLayerWidth = 0;
            dataSetArray[ic].collFactor = 0.0;
            dataSetArray[ic].deltaT = 0.0;
            dataSetArray[ic].compressible = false;
            dataSetArray[ic].withForcing = false;
         }

         boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageValuesArray3DPtr = block->getKernel()->getDataSet()->getAverageValues();
         if (AverageValuesArray3DPtr&&(blockParamStr.nx[0][0]>0)&&(blockParamStr.nx[0][1]>0)&&(blockParamStr.nx[0][2]>0)&&(blockParamStr.nx[0][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageValuesArray3DPtr->getDataVector().begin(), AverageValuesArray3DPtr->getDataVector().end());

         boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageVelocityArray3DPtr = block->getKernel()->getDataSet()->getAverageVelocity();
         if (AverageVelocityArray3DPtr&&(blockParamStr.nx[1][0]>0)&&(blockParamStr.nx[1][1]>0)&&(blockParamStr.nx[1][2]>0)&&(blockParamStr.nx[1][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageVelocityArray3DPtr->getDataVector().begin(), AverageVelocityArray3DPtr->getDataVector().end());

         boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageFluctArray3DPtr = block->getKernel()->getDataSet()->getAverageFluctuations();
         if (AverageFluctArray3DPtr&&(blockParamStr.nx[2][0]>0)&&(blockParamStr.nx[2][1]>0)&&(blockParamStr.nx[2][2]>0)&&(blockParamStr.nx[2][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageFluctArray3DPtr->getDataVector().begin(), AverageFluctArray3DPtr->getDataVector().end());

         boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > AverageTripleArray3DPtr = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
         if (AverageTripleArray3DPtr&&(blockParamStr.nx[3][0]>0)&&(blockParamStr.nx[3][1]>0)&&(blockParamStr.nx[3][2]>0)&&(blockParamStr.nx[3][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), AverageTripleArray3DPtr->getDataVector().begin(), AverageTripleArray3DPtr->getDataVector().end());

         boost::shared_ptr< CbArray4D<LBMReal, IndexerX4X3X2X1> > ShearStressValArray3DPtr = block->getKernel()->getDataSet()->getShearStressValues();
         if (ShearStressValArray3DPtr&&(blockParamStr.nx[4][0]>0)&&(blockParamStr.nx[4][1]>0)&&(blockParamStr.nx[4][2]>0)&&(blockParamStr.nx[4][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), ShearStressValArray3DPtr->getDataVector().begin(), ShearStressValArray3DPtr->getDataVector().end());

         boost::shared_ptr< CbArray3D<LBMReal, IndexerX3X2X1> > RelaxationFactor3DPtr = block->getKernel()->getDataSet()->getRelaxationFactor();
         if (RelaxationFactor3DPtr&&(blockParamStr.nx[5][0]>0)&&(blockParamStr.nx[5][1]>0)&&(blockParamStr.nx[5][2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), RelaxationFactor3DPtr->getDataVector().begin(), RelaxationFactor3DPtr->getDataVector().end());

         boost::shared_ptr< D3Q27EsoTwist3DSplittedVector > D3Q27EsoTwist3DSplittedVectorPtr = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(block->getKernel()->getDataSet()->getFdistributions());
         CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getLocalDistributions();
         if (localDistributions&&(blockParamStr.nx[6][0]>0)&&(blockParamStr.nx[6][1]>0)&&(blockParamStr.nx[6][2]>0)&&(blockParamStr.nx[6][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), localDistributions->getDataVector().begin(), localDistributions->getDataVector().end());

         CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getNonLocalDistributions();
         if (nonLocalDistributions&&(blockParamStr.nx[7][0]>0)&&(blockParamStr.nx[7][1]>0)&&(blockParamStr.nx[7][2]>0)&&(blockParamStr.nx[7][3]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), nonLocalDistributions->getDataVector().begin(), nonLocalDistributions->getDataVector().end());

         CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions = D3Q27EsoTwist3DSplittedVectorPtr->getZeroDistributions();
         if (zeroDistributions&&(blockParamStr.nx[8][0]>0)&&(blockParamStr.nx[8][1]>0)&&(blockParamStr.nx[8][2]>0))
            doubleValuesArray.insert(doubleValuesArray.end(), zeroDistributions->getDataVector().begin(), zeroDistributions->getDataVector().end());

         ic++;
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   // write to the file
   // all processes calculate their offsets (quantity of bytes that the process is going to write) 
   // and notify the next process (with the rank = rank + 1)
   size_t view_offset = size*sizeof(int);
   size_t next_view_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_view_offset = view_offset+blocksCount * (sizeof(DataSet)+blockParamStr.doubleCountInBlock*sizeof(double));
         MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&view_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_view_offset = view_offset+blocksCount * (sizeof(DataSet)+blockParamStr.doubleCountInBlock*sizeof(double));
         if (rank<size-1)
            MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

   MPI_File file_handler;
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpDataSet.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   // each process writes the quantity of it's blocks
   MPI_File_write_at_all(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // each process writes data identifying blocks
   MPI_File_write_at_all(file_handler, view_offset, dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
   // each process writes the dataSet arrays
   MPI_File_write_at_all(file_handler, view_offset+blocksCount*sizeof(DataSet), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeDataSet time: "<<finish-start<<" s");
   }

   delete[] dataSetArray;
}

void MPIIORestartCoProcessor::writeBoundaryConds(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds start collect data rank = "<<rank);
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

   int ic = 0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blocksVector[level])  // all the blocks of the current level
      {
         BCArray3D bcArr = block->getKernel()->getBCProcessor()->getBCArray();

         bcAddArray[ic].x1 = block->getX1(); // coordinates of the block needed to find it while regenerating the grid
         bcAddArray[ic].x2 = block->getX2();
         bcAddArray[ic].x3 = block->getX3();
         bcAddArray[ic].level = block->getLevel();
         bcAddArray[ic].boundCond_count = 0; // how many BoundaryConditions in this block
         bcAddArray[ic].indexContainer_count = 0;  // how many indexContainer-values in this block

         for (int bc = 0; bc<bcArr.getBCVectorSize(); bc++)
         {
            BoundaryCondition* bouCond = new BoundaryCondition();
            if (bcArr.bcvector[bc]==NULL)
            {
               memset(bouCond, 0, sizeof(BoundaryCondition));
            }
            else
            {
               bouCond->noslipBoundaryFlags = bcArr.bcvector[bc]->getNoSlipBoundary();
               bouCond->slipBoundaryFlags = bcArr.bcvector[bc]->getSlipBoundary();
               bouCond->velocityBoundaryFlags = bcArr.bcvector[bc]->getVelocityBoundary();
               bouCond->densityBoundaryFlags = bcArr.bcvector[bc]->getDensityBoundary();
               bouCond->wallModelBoundaryFlags = bcArr.bcvector[bc]->getWallModelBoundary();
               bouCond->bcVelocityX1 = bcArr.bcvector[bc]->getBoundaryVelocityX1();
               bouCond->bcVelocityX2 = bcArr.bcvector[bc]->getBoundaryVelocityX2();
               bouCond->bcVelocityX3 = bcArr.bcvector[bc]->getBoundaryVelocityX3();
               bouCond->bcDensity = bcArr.bcvector[bc]->getBoundaryDensity();
               bouCond->bcLodiDensity = bcArr.bcvector[bc]->getDensityLodiDensity();
               bouCond->bcLodiVelocityX1 = bcArr.bcvector[bc]->getDensityLodiVelocityX1();
               bouCond->bcLodiVelocityX2 = bcArr.bcvector[bc]->getDensityLodiVelocityX2();
               bouCond->bcLodiVelocityX3 = bcArr.bcvector[bc]->getDensityLodiVelocityX3();
               bouCond->bcLodiLentgh = bcArr.bcvector[bc]->getDensityLodiLength();
               bouCond->nx1 = bcArr.bcvector[bc]->nx1;
               bouCond->nx2 = bcArr.bcvector[bc]->nx2;
               bouCond->nx3 = bcArr.bcvector[bc]->nx3;
               for (int iq = 0; iq<26; iq++)
                  bouCond->q[iq] = bcArr.bcvector[bc]->getQ(iq);
               bouCond->algorithmType = bcArr.bcvector[bc]->getBcAlgorithmType();
            }

            bcVector.push_back(*bouCond);
            bcAddArray[ic].boundCond_count++;
            count_boundCond++;
         }

         bcindexmatrixV.insert(bcindexmatrixV.end(), bcArr.bcindexmatrix.getDataVector().begin(), bcArr.bcindexmatrix.getDataVector().end());

         indexContainerV.insert(indexContainerV.end(), bcArr.indexContainer.begin(), bcArr.indexContainer.end());
      
         bcAddArray[ic].indexContainer_count = static_cast<int>(bcArr.indexContainer.size());
         count_indexContainer += bcAddArray[ic].indexContainer_count;

         ic++;
      }
   }

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

   byteCount = bcBlockCount * BLOCK_SIZE*sizeof(BoundaryCondition)+blocksCount*sizeof(BCAdd)+sizeof(int) * (blocksCount * blockParamStr.bcindexmatrix_count+count_indexContainer);

   // write to the file
   // all processes calculate their offsets (quantity of bytes that the process is going to write) 
   // and notify the next process (with the rank = rank + 1)
   size_t view_offset = 3*size*sizeof(int);
   size_t next_view_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_view_offset = view_offset+byteCount;
         MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&view_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_view_offset = view_offset+byteCount;
         if (rank<size-1)
            MPI_Send(&next_view_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds start MPI IO rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_Info info = MPI_INFO_NULL;
   //MPI_Info_create (&info);
   //MPI_Info_set(info,"romio_cb_write","enable");
   //MPI_Info_set(info,"cb_buffer_size","4194304");
   //MPI_Info_set(info,"striping_unit","4194304");

   MPI_File file_handler;
   std::string filename = path+"/mpi_io_cp/mpi_io_cp_"+UbSystem::toString(step)+"/cpBC.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &file_handler);
   if (rc!=MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file "+filename);

   // each process writes the quantity of it's blocks
   MPI_File_write_at_all(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);	//	blocks quantity
   // each process writes the quantity of "big blocks" of BLOCK_SIZE of boundary conditions
   MPI_File_write_at_all(file_handler, (rank+size)*sizeof(int), &bcBlockCount, 1, MPI_INT, MPI_STATUS_IGNORE); // quantity of BoundConds / BLOCK_SIZE
   // each process writes the quantity of indexContainer elements in all blocks
   MPI_File_write_at_all(file_handler, (rank+2*size)*sizeof(int), &count_indexContainer, 1, MPI_INT, MPI_STATUS_IGNORE); // quantity of indexContainer	

   // each process writes data identifying the blocks
   MPI_File_write_at_all(file_handler, view_offset, bcAddArray, blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
   // each process writes boundary conditions
   if (bcVector.size()>0)
   {
      MPI_File_write_at(file_handler, view_offset+blocksCount*sizeof(BCAdd), &bcVector[0], bcBlockCount, boundCondType1000, MPI_STATUS_IGNORE);
   }
   // each process writes bcindexmatrix values
   MPI_File_write_at_all(file_handler, view_offset+blocksCount*sizeof(BCAdd)+bcBlockCount * BLOCK_SIZE*sizeof(BoundaryCondition), &bcindexmatrixV[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
   // each process writes indexContainer values
   MPI_File_write_at_all(file_handler, view_offset+blocksCount*sizeof(BCAdd)+bcBlockCount * BLOCK_SIZE*sizeof(BoundaryCondition)+blocksCount * blockParamStr.bcindexmatrix_count*sizeof(int), &indexContainerV[0], count_indexContainer, MPI_INT, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::writeBoundaryConds time: "<<finish-start<<" s");
   }

   delete[] bcAddArray;
}

//------------------------------------------- READ -----------------------------------------------
void MPIIORestartCoProcessor::restart(int step)
{
   if (comm->isRoot()) UBLOG(logINFO, "MPIIORestartCoProcessor restart step: "<<step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - start");
   readBlocks(step);
   readDataSet(step);
   readBoundaryConds(step);
   if (comm->isRoot()) UBLOG(logINFO, "Load check point - end");
   this->reconnect(grid);
}

void MPIIORestartCoProcessor::readBlocks(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //MPI_Comm_size(MPI_COMM_WORLD, &size);
   size = 1;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBlocks start MPI IO rank = "<<rank);
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
   size_t read_offset = size*sizeof(int);
   //size_t read_offset = sizeof(int);
   size_t next_read_offset = 0;

   //if (size>1)
   //{
   //   if (rank==0)
   //   {
   //      next_read_offset = read_offset+sizeof(GridParam)+sizeof(BlockParam)+blocksCount*sizeof(Block3d);
   //      MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
   //   }
   //   else
   //   {
   //      MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   //      next_read_offset = read_offset+sizeof(GridParam)+sizeof(BlockParam)+blocksCount*sizeof(Block3d);
   //      if (rank<size-1)
   //         MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
   //   }
   //}

   GridParam* gridParameters = new GridParam;

   // read parameters of the grid
   MPI_File_read_at(file_handler, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   // read parameters of a block
   MPI_File_read_at(file_handler, read_offset+sizeof(GridParam), &blockParamStr, 1, blockParamType, MPI_STATUS_IGNORE);
   // read all the blocks
   MPI_File_read_at(file_handler, read_offset+sizeof(GridParam)+sizeof(BlockParam), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBlocks time: "<<finish-start<<" s");
   }

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBlocks start of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   // clear the grid
   std::vector<Block3DPtr> blocksVector;
   grid->getBlocks(0, blocksVector);
   int del = 0;
   BOOST_FOREACH(Block3DPtr block, blocksVector)
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
      //    LBMKernelPtr kernel(new CompressibleCumulantLBMKernel(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, CompressibleCumulantLBMKernel::NORMAL));
          //kernel->setCollisionFactor(block3dArray[n].collFactor);
          //kernel->setGhostLayerWidth(block3dArray[n].ghostLayerWidth);
          //kernel->setCompressible(block3dArray[n].compressible);
          //kernel->setWithForcing(block3dArray[n].withForcing);
          //kernel->setDeltaT(block3dArray[n].deltaT);

      Block3DPtr block(new Block3D(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, block3dArray[n].level));
      block->setActive(block3dArray[n].active);
      block->setBundle(block3dArray[n].bundle);
      block->setRank(block3dArray[n].rank);
      block->setLocalRank(block3dArray[n].lrank);
      block->setGlobalID(block3dArray[n].globalID);
      block->setLocalID(block3dArray[n].localID);
      block->setPart(block3dArray[n].part),
         block->setLevel(block3dArray[n].level);
      //block->setInterpolationFlagCF(block3dArray[n].interpolationFlagCF);
      //block->setInterpolationFlagFC(block3dArray[n].interpolationFlagFC);
      block->interpolationFlagCF = block3dArray[n].interpolationFlagCF;
      block->interpolationFlagFC = block3dArray[n].interpolationFlagFC;

      //if(rank == block3dArray[n].rank)
      //	block->setKernel(kernel);

      grid->addBlock(block);
   }

   // define MPI_types depending on the block-specific information
   MPI_Type_contiguous(blockParamStr.doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
   MPI_Type_commit(&dataSetDoubleType);

   MPI_Type_contiguous(blockParamStr.bcindexmatrix_count, MPI_INT, &bcindexmatrixType);
   MPI_Type_commit(&bcindexmatrixType);

   mpiTypeFreeFlag = true;

   delete gridParameters;
   delete[] block3dArray;

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBlocks end of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
}

void MPIIORestartCoProcessor::readDataSet(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet start MPI IO rank = "<<rank);
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
   MPI_File_read_at_all(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);

   DataSet* dataSetArray = new DataSet[blocksCount];
   std::vector<double> doubleValuesArray(blocksCount * blockParamStr.doubleCountInBlock); // double-values in all blocks 

   // calculate the read offset
   size_t read_offset = size*sizeof(int);
   size_t next_read_offset = 0;

   if(size > 1)
   {
   	if(rank == 0)
   	{
   		next_read_offset = read_offset + blocksCount * (sizeof(DataSet) + blockParamStr.doubleCountInBlock * sizeof(double));
   		MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
   	}
   	else
   	{
   		MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_read_offset = read_offset + blocksCount * (sizeof(DataSet) + blockParamStr.doubleCountInBlock * sizeof(double));
   		if(rank < size - 1)
   			MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank + 1, 5, MPI_COMM_WORLD);
   	}
   }

   MPI_File_read_at_all(file_handler, read_offset, dataSetArray, blocksCount, dataSetType, MPI_STATUS_IGNORE);
   MPI_File_read_at_all(file_handler, read_offset+blocksCount*sizeof(DataSet), &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);
   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet time: "<<finish-start<<" s");
      UBLOG(logINFO, "MPIIORestartCoProcessor::readDataSet start of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }

   size_t index = 0, nextVectorSize = 0;
   std::vector<double> vectorsOfValues[9];
   for (int n = 0; n<blocksCount; n++)
   {
      for (int b = 0; b<9; b++) // assign approciate vectors for 9 dataSet arrays
      {
         nextVectorSize = blockParamStr.nx[b][0]*blockParamStr.nx[b][1]*blockParamStr.nx[b][2]*blockParamStr.nx[b][3];
         vectorsOfValues[b].assign(doubleValuesArray.data()+index, doubleValuesArray.data()+index+nextVectorSize);
         index += nextVectorSize;
      }

      // fill dataSet arrays
      AverageValuesArray3DPtr mAverageValues;
      if ((blockParamStr.nx[0][0]==0)&&(blockParamStr.nx[0][1]==0)&&(blockParamStr.nx[0][2]==0)&&(blockParamStr.nx[0][3]==0))
         mAverageValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[0], blockParamStr.nx[0][0], blockParamStr.nx[0][1], blockParamStr.nx[0][2], blockParamStr.nx[0][3]));

      AverageValuesArray3DPtr mAverageVelocity;
      if ((blockParamStr.nx[1][0]==0)&&(blockParamStr.nx[1][1]==0)&&(blockParamStr.nx[1][2]==0)&&(blockParamStr.nx[1][3]==0))
         mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageVelocity = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[1], blockParamStr.nx[1][0], blockParamStr.nx[1][1], blockParamStr.nx[1][2], blockParamStr.nx[1][3]));

      AverageValuesArray3DPtr mAverageFluktuations;
      if ((blockParamStr.nx[2][0]==0)&&(blockParamStr.nx[2][1]==0)&&(blockParamStr.nx[2][2]==0)&&(blockParamStr.nx[2][3]==0))
         mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageFluktuations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[2], blockParamStr.nx[2][0], blockParamStr.nx[2][1], blockParamStr.nx[2][2], blockParamStr.nx[2][3]));

      AverageValuesArray3DPtr mAverageTriplecorrelations;
      if ((blockParamStr.nx[3][0]==0)&&(blockParamStr.nx[3][1]==0)&&(blockParamStr.nx[3][2]==0)&&(blockParamStr.nx[3][3]==0))
         mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mAverageTriplecorrelations = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[3], blockParamStr.nx[3][0], blockParamStr.nx[3][1], blockParamStr.nx[3][2], blockParamStr.nx[3][3]));

      ShearStressValuesArray3DPtr mShearStressValues;
      if ((blockParamStr.nx[4][0]==0)&&(blockParamStr.nx[4][1]==0)&&(blockParamStr.nx[4][2]==0)&&(blockParamStr.nx[4][3]==0))
         mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr();
      else
         mShearStressValues = CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[4], blockParamStr.nx[4][0], blockParamStr.nx[4][1], blockParamStr.nx[4][2], blockParamStr.nx[4][3]));

      RelaxationFactorArray3DPtr mRelaxationFactor;
      if ((blockParamStr.nx[5][0]==0)&&(blockParamStr.nx[5][1]==0)&&(blockParamStr.nx[5][2]==0))
         mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr();
      else
         mRelaxationFactor = CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues[5], blockParamStr.nx[5][0], blockParamStr.nx[5][1], blockParamStr.nx[5][2]));

      DistributionArray3DPtr mFdistributions(new D3Q27EsoTwist3DSplittedVector(blockParamStr.nx1, blockParamStr.nx2, blockParamStr.nx3, -999.0));

      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[6], blockParamStr.nx[6][0], blockParamStr.nx[6][1], blockParamStr.nx[6][2], blockParamStr.nx[6][3])));
      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNonLocalDistributions(CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr(new CbArray4D<LBMReal, IndexerX4X3X2X1>(vectorsOfValues[7], blockParamStr.nx[7][0], blockParamStr.nx[7][1], blockParamStr.nx[7][2], blockParamStr.nx[7][3])));
      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setZeroDistributions(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<LBMReal, IndexerX3X2X1>(vectorsOfValues[8], blockParamStr.nx[8][0], blockParamStr.nx[8][1], blockParamStr.nx[8][2])));

      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX1(blockParamStr.nx1);
      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX2(blockParamStr.nx2);
      boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(mFdistributions)->setNX3(blockParamStr.nx3);

      DataSet3DPtr dataSetPtr = DataSet3DPtr(new DataSet3D());
      dataSetPtr->setAverageValues(mAverageValues);
      dataSetPtr->setAverageVelocity(mAverageVelocity);
      dataSetPtr->setAverageFluctuations(mAverageFluktuations);
      dataSetPtr->setAverageTriplecorrelations(mAverageTriplecorrelations);
      dataSetPtr->setShearStressValues(mShearStressValues);
      dataSetPtr->setRelaxationFactor(mRelaxationFactor);
      dataSetPtr->setFdistributions(mFdistributions);

      // find the nesessary block and fill it
      Block3DPtr block = grid->getBlock(dataSetArray[n].x1, dataSetArray[n].x2, dataSetArray[n].x3, dataSetArray[n].level);
      //LBMKernelPtr kernel(new CompressibleCumulantLBMKernel(0, 0, 0, CompressibleCumulantLBMKernel::NORMAL));
      LBMKernelPtr kernel(new CompressibleCumulantLBMKernel());
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

void MPIIORestartCoProcessor::readBoundaryConds(int step)
{
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds start MPI IO rank = "<<rank);
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
   // read count of blocks
   MPI_File_read_at_all(file_handler, rank*sizeof(int), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   // read count of big BoundaryCondition blocks
   MPI_File_read_at_all(file_handler, (rank+size)*sizeof(int), &dataCount1000, 1, MPI_INT, MPI_STATUS_IGNORE);
   // read count of indexContainer values in all blocks
   MPI_File_read_at_all(file_handler, (rank+2*size)*sizeof(int), &dataCount2, 1, MPI_INT, MPI_STATUS_IGNORE);

   size_t dataCount = dataCount1000 * BLOCK_SIZE;
   BCAdd* bcAddArray = new BCAdd[blocksCount];
   BoundaryCondition* bcArray = new BoundaryCondition[dataCount];
   BoundaryCondition* nullBouCond = new BoundaryCondition();
   memset(nullBouCond, 0, sizeof(BoundaryCondition));
   int* intArray1 = new int[blocksCount * blockParamStr.bcindexmatrix_count];
   int* intArray2 = new int[dataCount2];

   size_t read_offset = 3*size*sizeof(int);
   size_t next_read_offset = 0;

   if (size>1)
   {
      if (rank==0)
      {
         next_read_offset = read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+(blocksCount * blockParamStr.bcindexmatrix_count+dataCount2)*sizeof(int);
         MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, 1, 5, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(&read_offset, 1, MPI_LONG_LONG_INT, rank-1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         next_read_offset = read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+(blocksCount * blockParamStr.bcindexmatrix_count+dataCount2)*sizeof(int);
         if (rank<size-1)
            MPI_Send(&next_read_offset, 1, MPI_LONG_LONG_INT, rank+1, 5, MPI_COMM_WORLD);
      }
   }

   MPI_File_read_at_all(file_handler, read_offset, &bcAddArray[0], blocksCount, boundCondTypeAdd, MPI_STATUS_IGNORE);
   MPI_File_read_at_all(file_handler, read_offset+blocksCount*sizeof(BCAdd), &bcArray[0], dataCount1000, boundCondType1000, MPI_STATUS_IGNORE);
   MPI_File_read_at_all(file_handler, read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition), &intArray1[0], blocksCount, bcindexmatrixType, MPI_STATUS_IGNORE);
   MPI_File_read_at_all(file_handler, read_offset+blocksCount*sizeof(BCAdd)+dataCount*sizeof(BoundaryCondition)+blocksCount * blockParamStr.bcindexmatrix_count*sizeof(int), &intArray2[0], dataCount2, MPI_INT, MPI_STATUS_IGNORE);
   //MPI_File_sync(file_handler);

   MPI_File_close(&file_handler);

   if (comm->isRoot())
   {
      finish = MPI_Wtime();
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds time: "<<finish-start<<" s");
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds start of restore of data, rank = "<<rank);
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
         }

         bcVector.push_back(bc);
         index++;
      }

      for (int b1 = 0; b1<blockParamStr.bcindexmatrix_count; b1++)
         bcindexmatrixV.push_back(intArray1[index1++]);

      for (int b2 = 0; b2<bcAddArray[n].indexContainer_count; b2++)
         indexContainerV.push_back(intArray2[index2++]);

      CbArray3D<int, IndexerX3X2X1> bcim(bcindexmatrixV, blockParamStr.nx1, blockParamStr.nx2, blockParamStr.nx3);

      BCProcessorPtr bcProc(new BCProcessor());
      BCArray3D &bcArr = bcProc->getBCArray();
      bcArr.bcindexmatrix = bcim;
      bcArr.bcvector = bcVector;
      bcArr.indexContainer = indexContainerV;

      Block3DPtr block = grid->getBlock(bcAddArray[n].x1, bcAddArray[n].x2, bcAddArray[n].x3, bcAddArray[n].level);
      block->getKernel()->setBCProcessor(bcProc);
   }

   delete nullBouCond;
   delete[] bcArray;
   delete[] bcAddArray;
   delete[] intArray1;
   delete[] intArray2;
   
   if (comm->isRoot())
   {
      UBLOG(logINFO, "MPIIORestartCoProcessor::readBoundaryConds end of restore of data, rank = "<<rank);
      UBLOG(logINFO, "Physical Memory currently used by current process: "<<Utilities::getPhysMemUsedByMe()/1073741824.0<<" GB");
   }
}

