#include "WriteMacroscopicQuantitiesCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include "basics/writer/WbWriterVtkXmlASCII.h"

using namespace std;

WriteMacroscopicQuantitiesCoProcessor::WriteMacroscopicQuantitiesCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteMacroscopicQuantitiesCoProcessor::WriteMacroscopicQuantitiesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                                                                                 const std::string& path, WbWriter* const writer, 
                                                                                 LBMUnitConverterPtr conv, 
                                                                                 CommunicatorPtr comm)
                                                                                 : CoProcessor(grid, s),
                                                                                 path(path),
                                                                                 writer(writer),
                                                                                 conv(conv),
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
void WriteMacroscopicQuantitiesCoProcessor::init()
{

}
//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesCoProcessor::process(double step)
{
   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "WriteMacroscopicQuantitiesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesCoProcessor::collectData(double step)
{
   int istep = static_cast<int>(step);

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            addDataMQ(block);
         }
      }
   }

   string pfilePath, partPath, subfolder, cfilePath;

      subfolder = "mq"+UbSystem::toString(istep);
      pfilePath =  path+"/mq/"+subfolder;
      cfilePath =  path+"/mq/mq_collection";
      partPath = pfilePath+"/mq"+UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep);


   string partName = writer->writeOctsWithNodeData(partPath,nodes,cells,datanames,data);
   size_t found=partName.find_last_of("/");
   string piece = partName.substr(found+1);
   piece = subfolder + "/" + piece;

   vector<string> cellDataNames;
   CommunicatorPtr comm = Communicator::getInstance();
   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID() == comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath,pieces,datanames,cellDataNames);
      found=pname.find_last_of("/");
      piece = pname.substr(found+1);

      vector<string> filenames;
      filenames.push_back(piece);
      if (step == CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath,filenames,istep,false);
      } 
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath,filenames,istep,false);
      }
      UBLOG(logINFO,"WriteMacroscopicQuantitiesCoProcessor step: " << istep);
   }

   clearData();
}
//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesCoProcessor::clearData()
{
   nodes.clear();
   cells.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void WriteMacroscopicQuantitiesCoProcessor::addDataMQ(Block3DPtr block)
{
   UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
   double         dx           = grid->getDeltaX(block);

   double level = (double)block->getLevel();
   double blockID = (double)block->getGlobalID();

   //Diese Daten werden geschrieben:
   datanames.resize(0);
   datanames.push_back("Rho");
   datanames.push_back("Vx");
   datanames.push_back("Vy");
   datanames.push_back("Vz");
   //datanames.push_back("Press");
   //datanames.push_back("Level");
   //datanames.push_back("BlockID");

     

   data.resize(datanames.size());

   LBMKernelPtr kernel = block->getKernel();
   BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();     
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal vx1,vx2,vx3,rho;

   //knotennummerierung faengt immer bei 0 an!
   int SWB,SEB,NEB,NWB,SWT,SET,NET,NWT;

   //Funktionszeiger
   typedef void (*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/,LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);

   CalcMacrosFct calcMacros = NULL;

   if(block->getKernel()->getCompressible())
   {
      calcMacros = &D3Q27System::calcCompMacroscopicValues;
   }
   else
   {
      calcMacros = &D3Q27System::calcIncompMacroscopicValues;
   }

   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int maxX1 = (int)(distributions->getNX1());
   int maxX2 = (int)(distributions->getNX2());
   int maxX3 = (int)(distributions->getNX3());

   //int minX1 = 1;
   //int minX2 = 1;
   //int minX3 = 1;

   //int maxX1 = (int)(distributions->getNX1());
   //int maxX2 = (int)(distributions->getNX2());
   //int maxX3 = (int)(distributions->getNX3());

   //nummern vergeben und node vector erstellen + daten sammeln
   CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3,-1);
   maxX1 -= 2;
   maxX2 -= 2;
   maxX3 -= 2;

   //D3Q27BoundaryConditionPtr bcPtr;
   int nr = (int)nodes.size();
 
   for(size_t ix3=minX3; ix3<=maxX3; ix3++)
   {
      for(size_t ix2=minX2; ix2<=maxX2; ix2++)
      {
         for(size_t ix1=minX1; ix1<=maxX1; ix1++)
         {
            if(!bcArray.isUndefined(ix1,ix2,ix3) && !bcArray.isSolid(ix1,ix2,ix3))
            {
               int index = 0;
               nodeNumbers(ix1,ix2,ix3) = nr++;
               nodes.push_back( makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
                                            float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
                                            float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)) );

               distributions->getDistribution(f, ix1, ix2, ix3);
               calcMacros(f,rho,vx1,vx2,vx3);
               //double press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);

               if (UbMath::isNaN(rho) || UbMath::isInfinity(rho)) 
                  UB_THROW( UbException(UB_EXARGS,"rho is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                   ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                     //rho=999.0;
               //if (UbMath::isNaN(press) || UbMath::isInfinity(press)) 
               //   UB_THROW( UbException(UB_EXARGS,"press is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
               //   ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                  //press=999.0;
               if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1)) 
                  UB_THROW( UbException(UB_EXARGS,"vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                     //vx1=999.0;
               if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2)) 
                  UB_THROW( UbException(UB_EXARGS,"vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));
                     //vx2=999.0;
               if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3)) 
                  UB_THROW( UbException(UB_EXARGS,"vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                  ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));

               data[index++].push_back(rho);
               data[index++].push_back(vx1);
               data[index++].push_back(vx2);
               data[index++].push_back(vx3);
               
               //data[index++].push_back(rho * conv->getFactorDensityLbToW2() );
               //data[index++].push_back(vx1 * conv->getFactorVelocityLbToW2());
               //data[index++].push_back(vx2 * conv->getFactorVelocityLbToW2());
               //data[index++].push_back(vx3 * conv->getFactorVelocityLbToW2());
               //data[index++].push_back(press * conv->getFactorPressureLbToW2());
               //data[index++].push_back(level);
               //data[index++].push_back(blockID);
            }
         }
      }
   }
   maxX1 -= 1;
   maxX2 -= 1;
   maxX3 -= 1;
   //cell vector erstellen
   for(int ix3=minX3; ix3<=maxX3; ix3++)
   {
      for(int ix2=minX2; ix2<=maxX2; ix2++)
      {
         for(int ix1=minX1; ix1<=maxX1; ix1++)
         {
            if(   (SWB=nodeNumbers( ix1  , ix2,   ix3   )) >= 0
               && (SEB=nodeNumbers( ix1+1, ix2,   ix3   )) >= 0
               && (NEB=nodeNumbers( ix1+1, ix2+1, ix3   )) >= 0
               && (NWB=nodeNumbers( ix1  , ix2+1, ix3   )) >= 0 
               && (SWT=nodeNumbers( ix1  , ix2,   ix3+1 )) >= 0
               && (SET=nodeNumbers( ix1+1, ix2,   ix3+1 )) >= 0
               && (NET=nodeNumbers( ix1+1, ix2+1, ix3+1 )) >= 0
               && (NWT=nodeNumbers( ix1  , ix2+1, ix3+1 )) >= 0                )
            {
               cells.push_back( makeUbTuple(SWB,SEB,NEB,NWB,SWT,SET,NET,NWT) );
            }
         }
      }
   }
}
