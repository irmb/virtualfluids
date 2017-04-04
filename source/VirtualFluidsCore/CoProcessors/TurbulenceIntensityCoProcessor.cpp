#include "TurbulenceIntensityCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/utilities/UbMath.h"

using namespace std;

TurbulenceIntensityCoProcessor::TurbulenceIntensityCoProcessor(Grid3DPtr grid, const std::string& path, 
                                                                                       WbWriter* const writer,
                                                                                       UbSchedulerPtr s, CommunicatorPtr comm)
                                                                                     : CoProcessor(grid, s),
                                                                                       path(path),
                                                                                       comm(comm),
                                                                                       writer(writer)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::init()
{
   gridRank  = grid->getRank();
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   blockVector.resize(maxInitLevel+1);

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);

      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         UbTupleInt3 nx = grid->getBlockNX();
         AverageValuesArray3DPtr averageValues = AverageValuesArray3DPtr(new AverageValuesArray3D(val<1>(nx)+1, val<2>(nx)+1, val<3>(nx)+1, 4, 0.0));
         block->getKernel()->getDataSet()->setAverageValues(averageValues);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::process(double step)
{
   calculateAverageValues(int(step));

   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "TurbulenceIntensityCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::collectData(double step)
{
   int istep = int(step);

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            addData(block);
         }
      }
   }

   string partName = writer->writeOctsWithNodeData(path+ UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep),nodes,cells,datanames,data);
   size_t found=partName.find_last_of("//");
   string piece = partName.substr(found+1);

   vector<string> cellDataNames;

   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID() == comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(path+"_"+UbSystem::toString(istep),pieces,datanames,cellDataNames);

      vector<string> filenames;
      filenames.push_back(pname);
      if (step == CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"_collection",filenames,istep,false);
      } 
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path+"_collection",filenames,istep,false);
      }
      UBLOG(logINFO,"TurbulenceIntensityCoProcessor step: " << istep);
   }

   clearData();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::clearData()
{
   nodes.clear();
   cells.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::addData(const Block3DPtr block)
{
   UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
   double         dx           = grid->getDeltaX(block);

   //Diese Daten werden geschrieben:
   datanames.resize(0);
   datanames.push_back("TI");

   data.resize(datanames.size());

   LBMKernelPtr kernel = block->getKernel();
   BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
   AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();
   //int ghostLayerWidth = kernel->getGhostLayerWidth();

   //knotennummerierung faengt immer bei 0 an!
   int SWB,SEB,NEB,NWB,SWT,SET,NET,NWT;

   int minX1 = 0;
   int minX2 = 0;
   int minX3 = 0;

   int maxX1 = int(distributions->getNX1());
   int maxX2 = int(distributions->getNX2());
   int maxX3 = int(distributions->getNX3());

   //nummern vergeben und node vector erstellen + daten sammeln
   CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3,-1);
   //D3Q27BoundaryConditionPtr bcPtr;
   int nr = (int)nodes.size();
 
   for(int ix3=minX3; ix3<maxX3-1; ix3++)
   {
      for(int ix2=minX2; ix2<maxX2-1; ix2++)
      {
         for(int ix1=minX1; ix1<maxX1-1; ix1++)
         {
            if(!bcArray.isUndefined(ix1,ix2,ix3) && !bcArray.isSolid(ix1,ix2,ix3))
            {
               int index = 0;
               nodeNumbers(ix1,ix2,ix3) = nr++;
               nodes.push_back( makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
                                            float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
                                            float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)) );

               //compute turbulence intensity
               double temp = (*av)(ix1,ix2,ix3,AvVxxyyzz)/
                            ((*av)(ix1,ix2,ix3,AvVx)*(*av)(ix1,ix2,ix3,AvVx)+
                             (*av)(ix1,ix2,ix3,AvVy)*(*av)(ix1,ix2,ix3,AvVy)+
                             (*av)(ix1,ix2,ix3,AvVz)*(*av)(ix1,ix2,ix3,AvVz));

               LBMReal ti = sqrt(temp);

               if (UbMath::isNaN(ti)) UB_THROW( UbException(UB_EXARGS,"TI is not a number (nan or -1.#IND), sqrt(temp), where temp = "+UbSystem::toString(temp)+", AvVx = "+UbSystem::toString((*av)(ix1,ix2,ix3,AvVx))+" AvVy = "+UbSystem::toString((*av)(ix1,ix2,ix3,AvVy))+" AvVz = "+UbSystem::toString((*av)(ix1,ix2,ix3,AvVz))));

               data[index++].push_back(ti);
            }
         }
      }
   }
   //cell vector erstellen
   for(int ix3=minX3; ix3<maxX3-1; ix3++)
   {
      for(int ix2=minX2; ix2<maxX2-1; ix2++)
      {
         for(int ix1=minX1; ix1<maxX1-1; ix1++)
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
//////////////////////////////////////////////////////////////////////////
void TurbulenceIntensityCoProcessor::calculateAverageValues(double timeStep)
{
   using namespace D3Q27System;

   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   LBMReal f[27];
   LBMReal vx, vy, vz;

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            LBMKernelPtr kernel = block->getKernel();
            BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
            DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
            AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();

            int minX1 = 0;
            int minX2 = 0;
            int minX3 = 0;

            int maxX1 = int(distributions->getNX1());
            int maxX2 = int(distributions->getNX2());
            int maxX3 = int(distributions->getNX3());

            for(int ix3=minX3; ix3<maxX3-1; ix3++)
            {
               for(int ix2=minX2; ix2<maxX2-1; ix2++)
               {
                  for(int ix1=minX1; ix1<maxX1-1; ix1++)
                  {
                     if(!bcArray.isUndefined(ix1,ix2,ix3) && !bcArray.isSolid(ix1,ix2,ix3))
                     {
                        //////////////////////////////////////////////////////////////////////////
                        //read distribution
                        ////////////////////////////////////////////////////////////////////////////
                        distributions->getDistribution(f, ix1, ix2, ix3);
                        //////////////////////////////////////////////////////////////////////////
                        //compute velocity
                        //////////////////////////////////////////////////////////////////////////
                        vx = f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
                        + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
                        + f[BSE] - f[BNW]; 

                        vy = f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
                        - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
                        + f[BNW]; 

                        vz = f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
                        + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
                        - f[BNW];
                        //////////////////////////////////////////////////////////////////////////
                        //compute average values
                        //////////////////////////////////////////////////////////////////////////
                        (*av)(ix1,ix2,ix3,AvVx) = ((*av)(ix1,ix2,ix3,AvVx)*timeStep + vx)/(timeStep+1.0);
                        (*av)(ix1,ix2,ix3,AvVy) = ((*av)(ix1,ix2,ix3,AvVy)*timeStep + vy)/(timeStep+1.0);
                        (*av)(ix1,ix2,ix3,AvVz) = ((*av)(ix1,ix2,ix3,AvVz)*timeStep + vz)/(timeStep+1.0);

                        (*av)(ix1,ix2,ix3,AvVxxyyzz) = ((vx-(*av)(ix1,ix2,ix3,AvVx))*(vx-(*av)(ix1,ix2,ix3,AvVx)) +
                                                         (vy-(*av)(ix1,ix2,ix3,AvVy))*(vy-(*av)(ix1,ix2,ix3,AvVy)) +
                                                         (vz-(*av)(ix1,ix2,ix3,AvVz))*(vz-(*av)(ix1,ix2,ix3,AvVz)) +
                                                         (*av)(ix1,ix2,ix3,AvVxxyyzz)*timeStep)/(timeStep+1.0);
                        //////////////////////////////////////////////////////////////////////////
                     }
                  }
               }
            }
         }
      }
   }
}

