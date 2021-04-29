#include "Init/ReadGeometry.h"

////////////////////////////////////////////////////////////////////////////////
void readGeometry(Parameter* para, VF::GPU::Communicator* comm, int lev, std::string geometryFile)
{
   int dataSizePerGPU = para->getParH(lev)->gridNX * para->getParH(lev)->gridNY * para->getParH(lev)->gridNZ;
   unsigned int dataSizeTotal = para->getParH(lev)->gridNX * para->getParH(lev)->gridNY * para->getParH(lev)->gridNZ * para->getNumprocs();
   unsigned int * dataRoot = NULL;
   
   unsigned int m, l;
   if(para->getMyID() == 0)	// Are we the root node?
   {
      dataRoot = new unsigned int[dataSizeTotal];
      VtkGeometryReader::readFile(geometryFile, dataRoot);
   }

   std::cout << "dataSizePerGPU size: " << dataSizePerGPU << "\n";
   
   unsigned int *dataGPU = new unsigned int[dataSizePerGPU];
   
   std::cout  << "distributeGeometry: start \n";
   comm->distributeGeometry(dataRoot, dataGPU, dataSizePerGPU);
   std::cout  << "distributeGeometry: end \n";

   l=0;
   for(unsigned int k=STARTOFFZ ; k < para->getParH(lev)->gridNZ + STARTOFFZ ; k++)
   {
      for(unsigned int j=STARTOFFY ; j < para->getParH(lev)->gridNY + STARTOFFY ; j++)
      {
         for(unsigned int i=STARTOFFX ; i < para->getParH(lev)->gridNX + STARTOFFX ; i++)
         {
            m = para->getParH(lev)->nx* (para->getParH(lev)->ny * k + j) + i;
            para->getParH(lev)->geo[m] = dataGPU[l];
            l++;
         }
      }
   }

   if(para->getMyID() == 0)
      delete [] dataRoot;
   
   delete [] dataGPU;
}
