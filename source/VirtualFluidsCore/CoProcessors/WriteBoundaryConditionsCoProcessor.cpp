#include "WriteBoundaryConditionsCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include <vector>
#include <string>

#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "WbWriter.h"
#include "UbScheduler.h"
#include "CbArray3D.h"
#include "BCArray3D.h"

using namespace std;

WriteBoundaryConditionsCoProcessor::WriteBoundaryConditionsCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteBoundaryConditionsCoProcessor::WriteBoundaryConditionsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
   const std::string& path, WbWriter* const writer, SPtr<Communicator> comm)
   : CoProcessor(grid, s),
   path(path),
   writer(writer),
   comm(comm)
{
   gridRank = comm->getProcessID();
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   blockVector.resize(maxInitLevel+1);

   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);
   }
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
      collectData(step);

   UBLOG(logDEBUG3, "WriteBoundaryConditionsCoProcessor::update:"<<step);
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsCoProcessor::collectData(double step)
{
   int istep = static_cast<int>(step);

   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      for(SPtr<Block3D> block : blockVector[level])
      {
         if (block)
         {
            addDataGeo(block);
         }
      }
   }

   string pfilePath, partPath, subfolder, cfilePath;

   subfolder = "bc"+UbSystem::toString(istep);
   pfilePath = path+"/bc/"+subfolder;
   cfilePath = path+"/bc/bc_collection";
   partPath = pfilePath+"/bc"+UbSystem::toString(gridRank)+"_"+UbSystem::toString(istep);


   string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
   size_t found = partName.find_last_of("/");
   string piece = partName.substr(found+1);
   piece = subfolder+"/"+piece;

   vector<string> cellDataNames;
   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID()==comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
      found = pname.find_last_of("/");
      piece = pname.substr(found+1);

      vector<string> filenames;
      filenames.push_back(piece);
      if (step==CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
      }
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
      }
      UBLOG(logINFO, "WriteBoundaryConditionsCoProcessor step: "<<istep);
   }

   clearData();
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsCoProcessor::clearData()
{
   nodes.clear();
   cells.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void WriteBoundaryConditionsCoProcessor::addDataGeo(SPtr<Block3D> block)
{
   UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
   double         dx = grid->getDeltaX(block);

   double level = (double)block->getLevel();

   //Diese Daten werden geschrieben:
   datanames.resize(0);
   datanames.push_back("Boundary Conditions");
   datanames.push_back("Geometry");
   datanames.push_back("Level");
   //datanames.push_back("Interface CF");

   data.resize(datanames.size());

   SPtr<ILBMKernel> kernel = block->getKernel();
   SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

   //knotennummerierung faengt immer bei 0 an!
   int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int maxX1 = (int)bcArray->getNX1();
   int maxX2 = (int)bcArray->getNX2();
   int maxX3 = (int)bcArray->getNX3();

   //nummern vergeben und node vector erstellen + daten sammeln
   CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
   //D3Q27BoundaryConditionPtr bcPtr;
   int nr = (int)nodes.size();

   maxX1 -= 1;
   maxX2 -= 1;
   maxX3 -= 1;

   int s=0;

   

   for (size_t ix3 = minX3; ix3<=maxX3; ix3++)
   {
      for (size_t ix2 = minX2; ix2<=maxX2; ix2++)
      {
         for (size_t ix1 = minX1; ix1<=maxX1; ix1++)
         {
            if (!bcArray->isUndefined(ix1, ix2, ix3))
            {
               //int index = 0;
               nodeNumbers(ix1, ix2, ix3) = nr++;
               nodes.push_back(makeUbTuple(float(val<1>(org)-val<1>(nodeOffset)+ix1*dx),
                  float(val<2>(org)-val<2>(nodeOffset)+ix2*dx),
                  float(val<3>(org)-val<3>(nodeOffset)+ix3*dx)));



               if (!bcArray->hasBC(ix1, ix2, ix3))
               {
                  data[0].push_back(0.0);
               }
               else if (bcArray->getBC(ix1, ix2, ix3)->hasNoSlipBoundary())
                  data[0].push_back(1.0);
               else if (bcArray->getBC(ix1, ix2, ix3)->hasVelocityBoundary())
                  data[0].push_back(2.0);
               else if (bcArray->getBC(ix1, ix2, ix3)->hasDensityBoundary())
                  data[0].push_back(3.0);
               else if (bcArray->getBC(ix1, ix2, ix3)->hasSlipBoundary())
                  data[0].push_back(4.0);
               else
                  data[0].push_back(5.0);


               if (bcArray->isSolid(ix1, ix2, ix3))
               {
                  data[1].push_back(1.0);
               }
               else
               {
                  data[1].push_back(0.0);
               }
                  

               data[2].push_back(level);

               //if (bcArray->isInterfaceCF(ix1, ix2, ix3))
               //{
               //   data[3].push_back(1.0);
               //} 
               //else
               //{
               //   data[3].push_back(0.0);
               //}

            }
         }
      }
   }

   maxX1 -= 1;
   maxX2 -= 1;
   maxX3 -= 1;

   //cell vector erstellen
   for (int ix3 = minX3; ix3<=maxX3; ix3++)
   {
      for (int ix2 = minX2; ix2<=maxX2; ix2++)
      {
         for (int ix1 = minX1; ix1<=maxX1; ix1++)
         {
            if ((SWB = nodeNumbers(ix1, ix2, ix3))>=0
               &&(SEB = nodeNumbers(ix1+1, ix2, ix3))>=0
               &&(NEB = nodeNumbers(ix1+1, ix2+1, ix3))>=0
               &&(NWB = nodeNumbers(ix1, ix2+1, ix3))>=0
               &&(SWT = nodeNumbers(ix1, ix2, ix3+1))>=0
               &&(SET = nodeNumbers(ix1+1, ix2, ix3+1))>=0
               &&(NET = nodeNumbers(ix1+1, ix2+1, ix3+1))>=0
               &&(NWT = nodeNumbers(ix1, ix2+1, ix3+1))>=0)
            {
               cells.push_back(makeUbTuple(SWB, SEB, NEB, NWB, SWT, SET, NET, NWT));
            }
         }
      }
   }
}
