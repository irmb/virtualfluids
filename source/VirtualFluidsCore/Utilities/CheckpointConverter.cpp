#include "CheckpointConverter.h"
#include <MemoryUtil.h>
#include "BoundaryConditions.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "Communicator.h"
#include "CoordinateTransformation3D.h"
#include <stdio.h>

#define BLOCK_SIZE 1024

CheckpointConverter::CheckpointConverter(SPtr<Grid3D> grid, const std::string& path, SPtr<Communicator> comm) :
   grid(grid), path(path), comm(comm)
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

   MPI_Type_contiguous(BLOCK_SIZE, boundCondType, &boundCondType1000);
   MPI_Type_commit(&boundCondType1000);

   //---------------------------------------

   MPI_Type_contiguous(7, MPI_INT, &dataSetParamType);
   MPI_Type_commit(&dataSetParamType);

   //---------------------------------------

   MPI_Datatype typesDataSetRead[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksDataSetRead[3] = { 2, 5, 2 };
   MPI_Aint offsetsDataSetRead[3], lbDataSetRead, extentDataSetRead;

   offsetsDataSetRead[0] = 0;
   MPI_Type_get_extent(MPI_DOUBLE, &lbDataSetRead, &extentDataSetRead);
   offsetsDataSetRead[1] = blocksDataSetRead[0] * extentDataSetRead;

   MPI_Type_get_extent(MPI_INT, &lbDataSetRead, &extentDataSetRead);
   offsetsDataSetRead[2] = offsetsDataSetRead[1] + blocksDataSetRead[1] * extentDataSetRead;

   MPI_Type_create_struct(3, blocksDataSetRead, offsetsDataSetRead, typesDataSetRead, &dataSetTypeRead);
   MPI_Type_commit(&dataSetTypeRead);

   //-----------------------------------------------------------------------

   MPI_Datatype typesDataSetWrite[3] = { MPI_DOUBLE, MPI_INT, MPI_CHAR };
   int blocksDataSetWrite[3] = { 2, 2, 2 };
   MPI_Aint offsetsDataSetWrite[3], lbDataSetWrite, extentDataSetWrite;

   offsetsDataSetWrite[0] = 0;
   MPI_Type_get_extent(MPI_DOUBLE, &lbDataSetWrite, &extentDataSetWrite);
   offsetsDataSetWrite[1] = blocksDataSetWrite[0] * extentDataSetWrite;

   MPI_Type_get_extent(MPI_INT, &lbDataSetWrite, &extentDataSetWrite);
   offsetsDataSetWrite[2] = offsetsDataSetWrite[1] + blocksDataSetWrite[1] * extentDataSetWrite;

   MPI_Type_create_struct(3, blocksDataSetWrite, offsetsDataSetWrite, typesDataSetWrite, &dataSetTypeWrite);
   MPI_Type_commit(&dataSetTypeWrite);


   //-----------------------------------------------------------------------

   MPI_Type_contiguous(4, MPI_INT, &dataSetSmallTypeRead);
   MPI_Type_commit(&dataSetSmallTypeRead);

   MPI_Type_contiguous(1, MPI_INT, &dataSetSmallTypeWrite);
   MPI_Type_commit(&dataSetSmallTypeWrite);

   //-----------------------------------------------------------------------

   MPI_Type_contiguous(4, MPI_INT, &boundCondParamType);
   MPI_Type_commit(&boundCondParamType);

   MPI_Type_contiguous(6, MPI_INT, &boundCondTypeAddRead);
   MPI_Type_commit(&boundCondTypeAddRead);
  
   MPI_Type_contiguous(3, MPI_INT, &boundCondTypeAddWrite);
   MPI_Type_commit(&boundCondTypeAddWrite);

}

//////////////////////////////////////////////////////////////////////////
CheckpointConverter::~CheckpointConverter()
{
   MPI_Type_free(&gridParamType);
   MPI_Type_free(&block3dType);
   MPI_Type_free(&boundCondType);
   MPI_Type_free(&dataSetParamType);
   MPI_Type_free(&dataSetTypeRead);
   MPI_Type_free(&dataSetTypeWrite);
   MPI_Type_free(&boundCondType1000);
   MPI_Type_free(&dataSetSmallTypeRead);
   MPI_Type_free(&dataSetSmallTypeWrite);
   MPI_Type_free(&boundCondParamType);
   MPI_Type_free(&boundCondTypeAddRead);
   MPI_Type_free(&boundCondTypeAddWrite);
}

//------------------------------------------- READ -----------------------------------------------
void CheckpointConverter::convert(int step, int procCount)
{
   UBLOG(logINFO, "UtilConvertor::convert start ");

   convertBlocks(step, procCount);
   convertDataSet(step, procCount);
   convertBC(step, procCount);
   
   UBLOG(logINFO, "UtilConvertor::convert finish ");
}

void CheckpointConverter::convertBlocks(int step, int procCount)
{
   //   UBLOG(logINFO, "Physical Memory currently used by current process: " << Utilities::getPhysMemUsedByMe() / 1073741824.0 << " GB");

   double start, finish;
   start = MPI_Wtime();

   // file to read from
   MPI_File file_handlerR;
   std::string filenameR = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rcR = MPI_File_open(MPI_COMM_WORLD, filenameR.c_str(), MPI_MODE_RDONLY/*MPI_MODE_RDWR*/, MPI_INFO_NULL, &file_handlerR);
   if (rcR != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameR);

   // file to write to
   MPI_File file_handlerW;
   UbSystem::makeDirectory(path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step));
   std::string filenameW = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBlocks.bin";
   int rcW = MPI_File_open(MPI_COMM_WORLD, filenameW.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handlerW);
   if (rcW != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameW);

   // read count of blocks
   int blocksCount = 0;
   MPI_File_read_at(file_handlerR, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   Block3d* block3dArray = new Block3d[blocksCount];
   GridParam* gridParameters = new GridParam;

   // calculate the read offset
   procCount = 1; // readBlocks and writeBlocks in both MPIIORestartCoProcessor and MPIIOMigrationCoProcessor have size == 1!
   MPI_Offset read_offset = (MPI_Offset)(procCount * sizeof(int));

   // read parameters of the grid and blocks
   MPI_File_read_at(file_handlerR, read_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);

   // clear the grid
   std::vector<SPtr<Block3D>> blocksVector[25];
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   if (minInitLevel > -1)
   {
      int maxInitLevel = this->grid->getFinestInitializedLevel();
      for (int level = minInitLevel; level <= maxInitLevel; level++)
      {
         grid->getBlocks(level, blocksVector[level]);
         for (SPtr<Block3D> block : blocksVector[level])  //	blocks of the current level
            grid->deleteBlock(block);
      }
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

   // renumber blocks
   grid->renumberBlockIDs();

   // refresh globalID in all the blocks
   SPtr<Block3D> block;
   for (int n = 0; n < blocksCount; n++)
   {
      block = grid->getBlock(block3dArray[n].x1, block3dArray[n].x2, block3dArray[n].x3, block3dArray[n].level);
      block3dArray[n].globalID = block->getGlobalID();
   }

   // write all data to the file
   MPI_Offset write_offset = read_offset;

   //for (int ic = 0; ic < blocksCount; ic++)
   //{
   //   std::cout << "rank=" << 0 << ",ic=" << ic<< ",write_offset=" << write_offset + sizeof(GridParam);
   //   std::cout << " " << block3dArray[ic].x1 << " " << block3dArray[ic].x2 << " " << block3dArray[ic].x3 << " " << block3dArray[ic].bundle;
   //   std::cout << " " << block3dArray[ic].rank << " " << block3dArray[ic].lrank << " " << block3dArray[ic].part << " " << block3dArray[ic].globalID;
   //   std::cout << " " << block3dArray[ic].localID << " " << block3dArray[ic].level << " " << block3dArray[ic].interpolationFlagCF << " " << block3dArray[ic].interpolationFlagFC;
   //   std::cout << " " << block3dArray[ic].counter << " " << block3dArray[ic].active << std::endl;
   //}
   MPI_File_write_at(file_handlerW, 0, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handlerW, write_offset, gridParameters, 1, gridParamType, MPI_STATUS_IGNORE);
   MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(GridParam)), &block3dArray[0], blocksCount, block3dType, MPI_STATUS_IGNORE);

   MPI_File_close(&file_handlerR);
   MPI_File_close(&file_handlerW);

   finish = MPI_Wtime();
   UBLOG(logINFO, "UtilConvertor::convertBlocks time: " << finish - start << " s");

   delete gridParameters;
   delete[] block3dArray;
}

void CheckpointConverter::convertDataSet(int step, int procCount)
{
   // file to read from
   MPI_File file_handlerR;
   std::string filenameR = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rcR = MPI_File_open(MPI_COMM_WORLD, filenameR.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handlerR);
   if (rcR != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameR);

   // file to write to
   MPI_File file_handlerW;
   std::string filenameW = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   int rcW = MPI_File_open(MPI_COMM_WORLD, filenameW.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY/*MPI_MODE_RDWR*/, MPI_INFO_NULL, &file_handlerW);
   if (rcW != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameW);

   //ßß
   //MPI_File file_handler7;
   //std::string filename7 = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpDataSet.bin";
   //int rc7 = MPI_File_open(MPI_COMM_WORLD, filename7.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler7);
   //if (rc7 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename7);
   //ßß

   double start, finish;
   start = MPI_Wtime();

   int blocksCount = 0;
   dataSetParam dataSetParamStr1, dataSetParamStr2, dataSetParamStr3;
   DataSetRead* dataSetReadArray;
   DataSetWrite* dataSetWriteArray;
   size_t doubleCountInBlock;
   std::vector<double> doubleValuesArray;
   size_t sizeofOneDataSet;

   // calculate the read offset
   MPI_Offset read_offset = (MPI_Offset)(procCount * sizeof(int));
   MPI_Offset write_offset;

   for (int pc = 0; pc < procCount; pc++)
   {
      // read count of blocks and parameters of data arrays
      MPI_File_read_at(file_handlerR, (MPI_Offset)(pc * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, read_offset, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
      doubleCountInBlock = dataSetParamStr1.nx[0] * dataSetParamStr1.nx[1] * dataSetParamStr1.nx[2] * dataSetParamStr1.nx[3] +
         dataSetParamStr2.nx[0] * dataSetParamStr2.nx[1] * dataSetParamStr2.nx[2] * dataSetParamStr2.nx[3] +
         dataSetParamStr3.nx[0] * dataSetParamStr3.nx[1] * dataSetParamStr3.nx[2] * dataSetParamStr3.nx[3];

      MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
      MPI_Type_commit(&dataSetDoubleType);
      //ßß
      //DataSetWrite* dataSetArray123 = new DataSetWrite[blocksCount];
      //std::vector<double> doubleValuesArray123(blocksCount * doubleCountInBlock); // double-values in all blocks 
      //ßß

      dataSetReadArray = new DataSetRead[blocksCount];
      dataSetWriteArray = new DataSetWrite[blocksCount];
      doubleValuesArray.resize(blocksCount * doubleCountInBlock);

      //std::cout << "rank=" << pc << ", read_offset=" << read_offset + 3 * sizeof(dataSetParam) << std::endl;

      // read data
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 3 * sizeof(dataSetParam)), dataSetReadArray, blocksCount, dataSetTypeRead, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 3 * sizeof(dataSetParam) + blocksCount * sizeof(DataSetRead)), 
         &doubleValuesArray[0], blocksCount * doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);
      //MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 3 * sizeof(dataSetParam) + blocksCount * sizeof(DataSetRead)), 
      //   &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

      //std::cout << "rank=" << pc << ", doubleValuesArray[0]=" << doubleValuesArray[0] << ", doubleValuesArray[last]=" << doubleValuesArray[blocksCount * (int)doubleCountInBlock - 1] << std::endl;
      //std::cout << "rank1=" << pc << " " << std::setprecision(15) << doubleValuesArray[3046] << std::endl;
      //std::cout << "rank2=" << pc << " " << std::setprecision(15) << doubleValuesArray[3047] << std::endl;
      //std::cout << "rank3=" << pc << " " << std::setprecision(15) << doubleValuesArray[3048] << std::endl;
      //std::cout << "rank4=" << pc << " " << std::setprecision(15) << doubleValuesArray[3049] << std::endl;
      //std::cout << "rank5=" << pc << " " << std::setprecision(15) << doubleValuesArray[3050] << std::endl;

      // offset to read the data of the next process
      read_offset = read_offset + (MPI_Offset)(3 * sizeof(dataSetParam) + blocksCount * (sizeof(DataSetRead) + doubleCountInBlock * sizeof(double)));

      // write parameters of data arrays
      MPI_File_write_at(file_handlerW, (MPI_Offset)0, &dataSetParamStr1, 1, dataSetParamType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handlerW, (MPI_Offset)(sizeof(dataSetParam)), &dataSetParamStr2, 1, dataSetParamType, MPI_STATUS_IGNORE);
      MPI_File_write_at(file_handlerW, (MPI_Offset)(2 * sizeof(dataSetParam)), &dataSetParamStr3, 1, dataSetParamType, MPI_STATUS_IGNORE);
      
      sizeofOneDataSet = sizeof(DataSetWrite) + doubleCountInBlock * sizeof(double);

      // write blocks and their data arrays
      for (int nb = 0; nb < blocksCount; nb++)
      {
         SPtr<Block3D> block = grid->getBlock(dataSetReadArray[nb].x1, dataSetReadArray[nb].x2, dataSetReadArray[nb].x3, dataSetReadArray[nb].level);
         dataSetWriteArray[nb].globalID = block->getGlobalID();
         //std::cout << "rank=" << pc << ",nb=" << nb << dataSetReadArray[nb].x1 << " " << dataSetReadArray[nb].x2 << " " << dataSetReadArray[nb].x3 << " " << dataSetReadArray[nb].level;
         //std::cout << ", block="<< block->getGlobalID() <<std::endl;
         dataSetWriteArray[nb].ghostLayerWidth = dataSetReadArray[nb].ghostLayerWidth;
         dataSetWriteArray[nb].collFactor = dataSetReadArray[nb].collFactor;
         dataSetWriteArray[nb].deltaT = dataSetReadArray[nb].deltaT;
         dataSetWriteArray[nb].compressible = dataSetReadArray[nb].compressible;
         dataSetWriteArray[nb].withForcing = dataSetReadArray[nb].withForcing;

         write_offset = (MPI_Offset)(3 * sizeof(dataSetParam) + dataSetWriteArray[nb].globalID * sizeofOneDataSet);
         //std::cout << "rank=" << pc << ", block->getGlobalID()=" << block->getGlobalID() << ", doubleCountInBlock=" << doubleCountInBlock << std::endl;
         MPI_File_write_at(file_handlerW, write_offset, &dataSetWriteArray[nb], 1, dataSetTypeWrite, MPI_STATUS_IGNORE);
         MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(DataSetWrite)), &doubleValuesArray[nb * doubleCountInBlock], 
            doubleCountInBlock, MPI_DOUBLE, MPI_STATUS_IGNORE);
         //MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(DataSetWrite)), &doubleValuesArray[nb * doubleCountInBlock], 
         //   1, dataSetDoubleType, MPI_STATUS_IGNORE);
//ßß
         //MPI_Offset read_offset123 = (MPI_Offset)(3 * sizeof(dataSetParam) + block->getGlobalID() * sizeofOneDataSet);
         //MPI_File_read_at(file_handlerW, read_offset123, &dataSetArray123[nb], 1, dataSetTypeWrite, MPI_STATUS_IGNORE);
         //MPI_File_read_at(file_handlerW, (MPI_Offset)(read_offset123 + sizeof(DataSetWrite)), &doubleValuesArray123[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
         //std::cout << ", doubleValuesArray123[0]=" << doubleValuesArray123[0] << ", doubleValuesArray123[last]=" << doubleValuesArray123[doubleCountInBlock - 1] << std::endl;
        
         //for (int as = 0; as < doubleCountInBlock; as++)
         //   if (doubleValuesArray123[nb * doubleCountInBlock + as] != doubleValuesArray[nb * doubleCountInBlock + as])
         //      std::cout << "nb="<<nb<< ":"<<doubleValuesArray123[nb * doubleCountInBlock + as] << "!=" << doubleValuesArray[nb * doubleCountInBlock + as] << std::endl;
//ßß
      }
      //std::cout << "rank11=" << pc << " " << std::setprecision(15) << doubleValuesArray123[3046] << std::endl;
      //std::cout << "rank22=" << pc << " " << std::setprecision(15) << doubleValuesArray123[3047] << std::endl;
      //std::cout << "rank33=" << pc << " " << std::setprecision(15) << doubleValuesArray123[3048] << std::endl;
      //std::cout << "rank44=" << pc << " " << std::setprecision(15) << doubleValuesArray123[3049] << std::endl;
      //std::cout << "rank55=" << pc << " " << std::setprecision(15) << doubleValuesArray123[3050] << std::endl;

      MPI_Type_free(&dataSetDoubleType);
      delete[] dataSetReadArray;
      delete[] dataSetWriteArray;
   }
//ßß
   //MPI_File_close(&file_handler7);
//ßß
   MPI_File_close(&file_handlerR);

   MPI_File_sync(file_handlerW);
   MPI_File_close(&file_handlerW);

   DSArraysPresence arrPresence;
   MPI_File file_handler1;
   std::string filename1 = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
   int rc = MPI_File_open(MPI_COMM_WORLD, filename1.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handler1);
   if (rc != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename1);
   MPI_File_read_at(file_handler1, (MPI_Offset)0, &arrPresence, 6, MPI_CHAR, MPI_STATUS_IGNORE);
   MPI_File_close(&file_handler1);

   MPI_File file_handler2;
   std::string filename2 = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpArrays.bin";
   int rc2 = MPI_File_open(MPI_COMM_WORLD, filename2.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handler2);
   if (rc2 != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filename2);
   MPI_File_write_at(file_handler2, (MPI_Offset)0, &arrPresence, 6, MPI_CHAR, MPI_STATUS_IGNORE);
   MPI_File_sync(file_handler2);
   MPI_File_close(&file_handler2);

   std::string filenameRR = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step);
   std::string filenameWW = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step);

   if (arrPresence.isAverageDensityArrayPresent)
      convert___Array(step, procCount, filenameRR + "/cpAverageDensityArray.bin", filenameWW + "/cpAverageDensityArray.bin");

   if (arrPresence.isAverageVelocityArrayPresent)
      convert___Array(step, procCount, filenameRR + "/cpAverageVelocityArray.bin", filenameWW + "/cpAverageVelocityArray.bin");

   if (arrPresence.isAverageFluktuationsArrayPresent)
      convert___Array(step, procCount, filenameRR + "/cpAverageFluktuationsArray.bin", filenameWW + "/cpAverageFluktuationsArray.bin");

   if (arrPresence.isAverageTripleArrayPresent)
      convert___Array(step, procCount, filenameRR + "/cpAverageTripleArray.bin", filenameWW + "/cpAverageTripleArray.bin");

   if (arrPresence.isShearStressValArrayPresent)
      convert___Array(step, procCount, filenameRR + "/cpShearStressValArray.bin", filenameWW + "/cpShearStressValArray.bin");

   if (arrPresence.isRelaxationFactorPresent)
      convert___Array(step, procCount, filenameRR + "/cpRelaxationFactor.bin", filenameWW + "/cpRelaxationFactor.bin");

   finish = MPI_Wtime();
   UBLOG(logINFO, "UtilConvertor::convertDataSet time: " << finish - start << " s");

}

void CheckpointConverter::convert___Array(int step, int procCount, std::string filenameR, std::string filenameW)
{
   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   MPI_File file_handlerR;
   int rcR = MPI_File_open(MPI_COMM_WORLD, filenameR.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handlerR);
   if (rcR != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameR);

   MPI_File file_handlerW;
   int rcW = MPI_File_open(MPI_COMM_WORLD, filenameW.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handlerW);
   if (rcW != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameW);

   int blocksCount = 0;
   dataSetParam dataSetParamStr;
   memset(&dataSetParamStr, 0, sizeof(dataSetParam));
   DataSetSmallRead* dataSetSmallReadArray;
   DataSetSmallWrite* dataSetSmallWriteArray;
   int doubleCountInBlock;
   std::vector<double> doubleValuesArray;

   // calculate the read offset
   MPI_Offset read_offset = (MPI_Offset)(procCount * sizeof(int));
   MPI_Offset write_offset;
   size_t sizeofOneDataSet;

   for (int pc = 0; pc < procCount; pc++)
   {
      MPI_File_read_at(file_handlerR, (MPI_Offset)(pc * sizeof(int)), &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, read_offset, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

      dataSetSmallReadArray = new DataSetSmallRead[blocksCount];
      dataSetSmallWriteArray = new DataSetSmallWrite[blocksCount];
      doubleCountInBlock = dataSetParamStr.nx[0] * dataSetParamStr.nx[1] * dataSetParamStr.nx[2] * dataSetParamStr.nx[3];
      doubleValuesArray.resize(blocksCount * doubleCountInBlock); 
      MPI_Type_contiguous(doubleCountInBlock, MPI_DOUBLE, &dataSetDoubleType);
      MPI_Type_commit(&dataSetDoubleType);

      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + sizeof(dataSetParam)), dataSetSmallReadArray, blocksCount, dataSetSmallTypeRead, MPI_STATUS_IGNORE);
      if (doubleCountInBlock > 0)
         MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + sizeof(dataSetParam) + blocksCount * sizeof(DataSetSmallRead)), 
            &doubleValuesArray[0], blocksCount, dataSetDoubleType, MPI_STATUS_IGNORE);

      read_offset = read_offset + sizeof(dataSetParam) + blocksCount * (sizeof(DataSetSmallRead) + doubleCountInBlock * sizeof(double));

      sizeofOneDataSet = sizeof(DataSetSmallWrite) + doubleCountInBlock * sizeof(double);
      
      MPI_File_write_at(file_handlerW, 0, &dataSetParamStr, 1, dataSetParamType, MPI_STATUS_IGNORE);

      for (int nb = 0; nb < blocksCount; nb++)
      {
         SPtr<Block3D> block = grid->getBlock(dataSetSmallReadArray[nb].x1, dataSetSmallReadArray[nb].x2, dataSetSmallReadArray[nb].x3, dataSetSmallReadArray[nb].level);
         dataSetSmallWriteArray[nb].globalID = block->getGlobalID();    

         write_offset = (MPI_Offset)(sizeof(dataSetParam) + dataSetSmallWriteArray[nb].globalID * sizeofOneDataSet);
         MPI_File_write_at(file_handlerW, write_offset, &dataSetSmallWriteArray[nb], 1, dataSetSmallTypeWrite, MPI_STATUS_IGNORE);
         MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(DataSetSmallWrite)), 
            &doubleValuesArray[nb * doubleCountInBlock], 1, dataSetDoubleType, MPI_STATUS_IGNORE);
      }

      MPI_Type_free(&dataSetDoubleType);
      delete[] dataSetSmallReadArray;
      delete[] dataSetSmallWriteArray;
   }
   MPI_File_close(&file_handlerR);

   MPI_File_sync(file_handlerW);
   MPI_File_close(&file_handlerW);
   
   finish = MPI_Wtime();
   UBLOG(logINFO, "UtilConvertor::convert___Array time: " << finish - start << " s");

}

void CheckpointConverter::convertBC(int step, int procCount)
{
   // file to read from
   MPI_File file_handlerR;
   std::string filenameR = path + "/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   int rcR = MPI_File_open(MPI_COMM_WORLD, filenameR.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handlerR);
   if (rcR != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameR);

   // file to write to
   MPI_File file_handlerW;
   std::string filenameW = path + "/mig/mpi_io_cp/mpi_io_cp_" + UbSystem::toString(step) + "/cpBC.bin";
   int rcW = MPI_File_open(MPI_COMM_WORLD, filenameW.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handlerW);
   if (rcW != MPI_SUCCESS) throw UbException(UB_EXARGS, "couldn't open file " + filenameW);

   double start, finish;
   if (comm->isRoot()) start = MPI_Wtime();

   int blocksCount = 0;
   int dataCount1000 = 0;
   int dataCount2 = 0;
   size_t dataCount;
   BCAddRead* bcAddReadArray;
   BCAddWrite* bcAddWriteArray;
   BoundaryCondition* bcArray;
   BoundaryCondition* nullBouCond = new BoundaryCondition();
   memset(nullBouCond, 0, sizeof(BoundaryCondition));
   int* intArray1;
   int* intArray2;
   int indexBC;
   int indexC;

   MPI_Offset read_offset;
   MPI_Offset read_offset1 = (MPI_Offset)(procCount * (3 * sizeof(int) + sizeof(boundCondParam)));
   MPI_Offset write_offset = (MPI_Offset)(sizeof(boundCondParam) + grid->getNumberOfBlocks() * sizeof(size_t));
   MPI_Offset write_offsetIndex;

   for (int pc = 0; pc < procCount; pc++)
   {
      read_offset = (MPI_Offset)(pc * (3 * sizeof(int) + sizeof(boundCondParam)));

      // read count of blocks
      MPI_File_read_at(file_handlerR, read_offset, &blocksCount, 1, MPI_INT, MPI_STATUS_IGNORE);
      // read count of big BoundaryCondition blocks
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + sizeof(int)), &dataCount1000, 1, MPI_INT, MPI_STATUS_IGNORE);
      // read count of indexContainer values in all blocks
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 2 * sizeof(int)), &dataCount2, 1, MPI_INT, MPI_STATUS_IGNORE);
      // read count of bcindexmatrix values in every block
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset + 3 * sizeof(int)), &boundCondParamStr, 1, boundCondParamType, MPI_STATUS_IGNORE);

      bcAddReadArray = new BCAddRead[blocksCount];
      bcAddWriteArray = new BCAddWrite[blocksCount];
      dataCount = dataCount1000 * BLOCK_SIZE;
      bcArray = new BoundaryCondition[dataCount];
      intArray1 = new int[blocksCount * boundCondParamStr.bcindexmatrixCount];
      intArray2 = new int[dataCount2];
      MPI_Type_contiguous(boundCondParamStr.bcindexmatrixCount, MPI_INT, &bcindexmatrixType);
      MPI_Type_commit(&bcindexmatrixType);

      // read data arrays
      MPI_File_read_at(file_handlerR, read_offset1, bcAddReadArray, blocksCount * 6, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset1 + blocksCount * sizeof(BCAddRead)),
         &bcArray[0], dataCount1000, boundCondType1000, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset1 + blocksCount * sizeof(BCAddRead) + dataCount * sizeof(BoundaryCondition)),
         &intArray1[0], blocksCount * boundCondParamStr.bcindexmatrixCount, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read_at(file_handlerR, (MPI_Offset)(read_offset1 + blocksCount * sizeof(BCAddRead) + dataCount * sizeof(BoundaryCondition) + blocksCount * boundCondParamStr.bcindexmatrixCount * sizeof(int)),
         &intArray2[0], dataCount2, MPI_INT, MPI_STATUS_IGNORE);

      // offset to read the data of the next process
      read_offset1 = read_offset1 + blocksCount * sizeof(BCAddRead) + dataCount * sizeof(BoundaryCondition) + (blocksCount * boundCondParamStr.bcindexmatrixCount + dataCount2) * sizeof(int);

      MPI_File_write_at(file_handlerW, 0, &boundCondParamStr, 4, MPI_INT, MPI_STATUS_IGNORE);

      indexBC = 0;
      indexC = 0;
      // write blocks and their BC data arrays
      for (int nb = 0; nb < blocksCount; nb++)
      {
         SPtr<Block3D> block = grid->getBlock(bcAddReadArray[nb].x1, bcAddReadArray[nb].x2, bcAddReadArray[nb].x3, bcAddReadArray[nb].level);
         bcAddWriteArray[nb].globalID = block->getGlobalID(); 
         bcAddWriteArray[nb].boundCond_count = bcAddReadArray[nb].boundCond_count;             // how many BoundaryConditions in this block
         bcAddWriteArray[nb].indexContainer_count = bcAddReadArray[nb].indexContainer_count;        // how many indexContainer-values in this block

         write_offsetIndex = (MPI_Offset)(sizeof(boundCondParam) + bcAddWriteArray[nb].globalID * sizeof(size_t));
         MPI_File_write_at(file_handlerW, write_offsetIndex, &write_offset, 1, MPI_LONG_LONG_INT, MPI_STATUS_IGNORE);

         //MPI_File_write_at(file_handlerW, write_offset, &bcAddWriteArray[nb], 1, boundCondTypeAddWrite, MPI_STATUS_IGNORE);
         MPI_File_write_at(file_handlerW, write_offset, &bcAddWriteArray[nb], 3, MPI_INT, MPI_STATUS_IGNORE);
         if (bcAddWriteArray[nb].boundCond_count > 0)
            MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(BCAddWrite)), &bcArray[indexBC], bcAddWriteArray[nb].boundCond_count, boundCondType, MPI_STATUS_IGNORE);
         indexBC += bcAddWriteArray[nb].boundCond_count;

         if (boundCondParamStr.bcindexmatrixCount > 0)
            MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(BCAddWrite) + bcAddWriteArray[nb].boundCond_count * sizeof(BoundaryCondition)),
               &intArray1[nb * boundCondParamStr.bcindexmatrixCount], boundCondParamStr.bcindexmatrixCount, MPI_INT, MPI_STATUS_IGNORE);

         if (bcAddWriteArray[nb].indexContainer_count > 0)
            MPI_File_write_at(file_handlerW, (MPI_Offset)(write_offset + sizeof(BCAddWrite) + bcAddWriteArray[nb].boundCond_count * sizeof(BoundaryCondition) + boundCondParamStr.bcindexmatrixCount * sizeof(int)),
               &intArray2[indexC], bcAddWriteArray[nb].indexContainer_count, MPI_INT, MPI_STATUS_IGNORE);
         indexC += bcAddWriteArray[nb].indexContainer_count;

         write_offset += sizeof(BCAddWrite) + bcAddWriteArray[nb].boundCond_count * sizeof(BoundaryCondition) + boundCondParamStr.bcindexmatrixCount * sizeof(int) + bcAddWriteArray[nb].indexContainer_count * sizeof(int);
      }

      MPI_Type_free(&bcindexmatrixType);
      delete[] bcAddReadArray;
      delete[] bcAddWriteArray;
      delete[] bcArray;
      delete[] intArray1;
      delete[] intArray2;
   }

   delete nullBouCond;

   MPI_File_close(&file_handlerR);

   MPI_File_sync(file_handlerW);
   MPI_File_close(&file_handlerW);

   finish = MPI_Wtime();
   UBLOG(logINFO, "UtilConvertor::convertBC time: " << finish - start << " s");

}
