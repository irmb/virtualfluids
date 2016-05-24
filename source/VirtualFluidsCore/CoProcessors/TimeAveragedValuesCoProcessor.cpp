#include "TimeAveragedValuesCoProcessor.h"
#include "LBMKernel3D.h"
#include "SimulationParameters.h"
#include "D3Q27ETBCProcessor.h"
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <boost/foreach.hpp>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/utilities/UbMath.h"

using namespace std;

TimeAveragedValuesCoProcessor::TimeAveragedValuesCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
TimeAveragedValuesCoProcessor::TimeAveragedValuesCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer,
   UbSchedulerPtr s, int options)
   : CoProcessor(grid, s),
   path(path),
   writer(writer),
   options(options)
{
   init(s);
   volumeAveraging = false;
}
//////////////////////////////////////////////////////////////////////////
TimeAveragedValuesCoProcessor::TimeAveragedValuesCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer, UbSchedulerPtr s, int options,
   std::vector<int> levels, std::vector<double>& levelCoords, std::vector<double>& bounds)
   : CoProcessor(grid, s),
   path(path),
   writer(writer),
   options(options),
   levels(levels),
   levelCoords(levelCoords),
   bounds(bounds)
{
   init(s);
   volumeAveraging = true;
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::init(UbSchedulerPtr s)
{
   gridRank = grid->getRank();
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   blockVector.resize(maxInitLevel + 1);

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);

      if (blockVector[level].size() > 0)
         compressible = blockVector[level][0]->getKernel()->getCompressible();

      double begin = s->getMinBegin();
      double gridTimeStep = grid->getTimeStep();

      if (gridTimeStep == begin || gridTimeStep == 0)
      {
         BOOST_FOREACH(Block3DPtr block, blockVector[level])
         {
            UbTupleInt3 nx = grid->getBlockNX();
            if ((options&Velocity) == Velocity)
            {
               AverageVelocityArray3DPtr av = AverageVelocityArray3DPtr(new AverageVelocityArray3D(3, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
               block->getKernel()->getDataSet()->setAverageVelocity(av);
            }

            if ((options&Fluctuations) == Fluctuations)
            {
               AverageFluctuationsArray3DPtr af = AverageFluctuationsArray3DPtr(new AverageFluctuationsArray3D(6, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
               block->getKernel()->getDataSet()->setAverageFluctuations(af);
            }

            if ((options&Triplecorrelations) == Triplecorrelations)
            {
               AverageTriplecorrelationsArray3DPtr at = AverageTriplecorrelationsArray3DPtr(new AverageTriplecorrelationsArray3D(10, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
               block->getKernel()->getDataSet()->setAverageTriplecorrelations(at);
            }

         }
      }
   }

   //breakStep = scheduler->getMaxEnd() - scheduler->getMinBegin()+1;
   //UBLOG(logINFO, "breakSteps = " << breakStep);
   //breakStep = breakStep * (double)(1 << maxInitLevel);
   //breakStep = scheduler->getMaxEnd()*(double)(1 << maxInitLevel);
   //UBLOG(logINFO, "breakSteps = " << breakStep);

   iMinX1 = 1;
   iMinX2 = 1;
   iMinX3 = 1;
   
   lcounter = 0;
   
   levelFactor = 1 << maxInitLevel;
   minFineStep = (int)scheduler->getMinBegin() * levelFactor;
   maxFineStep = (int)scheduler->getMaxEnd() * levelFactor + levelFactor - 1;
   numberOfFineSteps = (maxFineStep - minFineStep) + 1;
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::process(double step)
{
   fineStep = (int)step * levelFactor + lcounter;

   if (scheduler->isDue(step))
   {
      //DEBUG/////////////////////
      //UBLOG(logINFO, "step = " << step << ", lcounter = " << lcounter << ", fineStep = " << fineStep << ", maxFineStep = " << maxFineStep << ", numberOfFineSteps = " << numberOfFineSteps);
      ////////////////////////////

      calculateSubtotal();

      if (fineStep == maxFineStep)
      {
         calculateAverageValues((double)numberOfFineSteps);
         collectData(step);
         if (volumeAveraging)
         {
            volumeAverage(step);
         }
         clearData();
      }

      if (lcounter == levelFactor-1)
      {
         lcounter = 0;
      }
      else
      {
         lcounter++;
      }
   }

   UBLOG(logDEBUG3, "AverageValuesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::collectData(double step)
{
   int istep = int(step);

   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            addData(block);
         }
      }
   }

   string pfilePath, partPath, subfolder, cfilePath;
   subfolder = "tav" + UbSystem::toString(istep);
   pfilePath = path + "/tav/" + subfolder;
   partPath = pfilePath + "/tav" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

   string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
   size_t found = partName.find_last_of("/");
   string piece = partName.substr(found + 1);
   piece = subfolder + "/" + piece;

   vector<string> cellDataNames;
   CommunicatorPtr comm = Communicator::getInstance();
   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID() == comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
      UBLOG(logINFO, "TimeAveragedValuesCoProcessor step: " << istep);
   }

   clearData();
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::clearData()
{
   nodes.clear();
   cells.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::addData(const Block3DPtr block)
{
   UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
   double         dx = grid->getDeltaX(block);

   //Diese Daten werden geschrieben:
   datanames.resize(0);

   if ((options&Velocity) == Velocity)
   {
      datanames.push_back("taVx");
      datanames.push_back("taVy");
      datanames.push_back("taVz");
   }

   if ((options&Fluctuations) == Fluctuations)
   {
      datanames.push_back("taVxx");
      datanames.push_back("taVyy");
      datanames.push_back("taVzz");
      datanames.push_back("taVxy");
      datanames.push_back("taVxz");
      datanames.push_back("taVyz");
   }

   if ((options&Triplecorrelations) == Triplecorrelations)
   {
      datanames.push_back("taVxxx");
      datanames.push_back("taVxxy");
      datanames.push_back("taVxxz");
      datanames.push_back("taVyyy");
      datanames.push_back("taVyyx");
      datanames.push_back("taVyyz");
      datanames.push_back("taVzzz");
      datanames.push_back("taVzzx");
      datanames.push_back("taVzzy");
      datanames.push_back("taVxyz");
   }


   //datanames.push_back("AvP");
   //datanames.push_back("AvPrms");


   data.resize(datanames.size());

   LBMKernel3DPtr kernel = block->getKernel();
   BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
   AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageVelocity();
   AverageValuesArray3DPtr af = kernel->getDataSet()->getAverageFluctuations();
   AverageValuesArray3DPtr at = kernel->getDataSet()->getAverageTriplecorrelations();
   //int ghostLayerWidth = kernel->getGhostLayerWidth();

   //knotennummerierung faengt immer bei 0 an!
   int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

   int minX1 = iMinX1;
   int minX2 = iMinX2;
   int minX3 = iMinX3;

   int maxX1 = int(distributions->getNX1());
   int maxX2 = int(distributions->getNX2());
   int maxX3 = int(distributions->getNX3());

   //nummern vergeben und node vector erstellen + daten sammeln
   CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);

   maxX1 -= 2;
   maxX2 -= 2;
   maxX3 -= 2;

   //D3Q27BoundaryConditionPtr bcPtr;

   int nr = (int)nodes.size();

   for (int ix3 = minX3; ix3 <= maxX3; ix3++)
   {
      for (int ix2 = minX2; ix2 <= maxX2; ix2++)
      {
         for (int ix1 = minX1; ix1 <= maxX1; ix1++)
         {
            if (!bcArray.isUndefined(ix1, ix2, ix3) && !bcArray.isSolid(ix1, ix2, ix3))
            {
               int index = 0;
               nodeNumbers(ix1, ix2, ix3) = nr++;
               nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
                  float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
                  float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)));

               if ((options&Velocity) == Velocity)
               {
                  data[index++].push_back((*av)(Vx, ix1, ix2, ix3));
                  data[index++].push_back((*av)(Vy, ix1, ix2, ix3));
                  data[index++].push_back((*av)(Vz, ix1, ix2, ix3));
               }

               if ((options&Fluctuations) == Fluctuations)
               {
                  data[index++].push_back((*af)(Vxx, ix1, ix2, ix3));
                  data[index++].push_back((*af)(Vyy, ix1, ix2, ix3));
                  data[index++].push_back((*af)(Vzz, ix1, ix2, ix3));
                  data[index++].push_back((*af)(Vxy, ix1, ix2, ix3));
                  data[index++].push_back((*af)(Vxz, ix1, ix2, ix3));
                  data[index++].push_back((*af)(Vyz, ix1, ix2, ix3));
               }

               if ((options&Triplecorrelations) == Triplecorrelations)
               {
                  data[index++].push_back((*at)(Vxxx, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vxxy, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vxxz, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vyyy, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vyyx, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vyyz, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vzzz, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vzzx, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vzzy, ix1, ix2, ix3));
                  data[index++].push_back((*at)(Vxyz, ix1, ix2, ix3));
               }

               //LBMReal vp=(*av)(P, ix1, ix2, ix3);
               //LBMReal vprms=(*av)(Prms, ix1, ix2, ix3);

            }
         }
      }
   }

   maxX1 -= 1;
   maxX2 -= 1;
   maxX3 -= 1;

   //cell vector erstellen
   for (int ix3 = minX3; ix3 <= maxX3; ix3++)
   {
      for (int ix2 = minX2; ix2 <= maxX2; ix2++)
      {
         for (int ix1 = minX1; ix1 <= maxX1; ix1++)
         {
            if ((SWB = nodeNumbers(ix1, ix2, ix3)) >= 0
               && (SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0
               && (NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0
               && (NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0
               && (SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0
               && (SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0
               && (NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0
               && (NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0)
            {
               cells.push_back(makeUbTuple(SWB, SEB, NEB, NWB, SWT, SET, NET, NWT));
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::calculateAverageValues(double timeSteps)
{
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      int i;
      //#ifdef _OPENMP
      //   #pragma omp parallel for 
      //#endif
            //BOOST_FOREACH(Block3DPtr block, blockVector[level])
      for (i = 0; i < blockVector[level].size(); i++)
      {
         Block3DPtr block = blockVector[level][i];
         if (block)
         {
            LBMKernel3DPtr kernel = block->getKernel();
            BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
            DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
            AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageVelocity();
            AverageValuesArray3DPtr af = kernel->getDataSet()->getAverageFluctuations();
            AverageValuesArray3DPtr at = kernel->getDataSet()->getAverageTriplecorrelations();

            int minX1 = iMinX1;
            int minX2 = iMinX2;
            int minX3 = iMinX3;

            int maxX1 = int(distributions->getNX1());
            int maxX2 = int(distributions->getNX2());
            int maxX3 = int(distributions->getNX3());

            maxX1 -= 2;
            maxX2 -= 2;
            maxX3 -= 2;

            LBMReal ux, uy, uz, uxx, uzz, uyy, uxy, uxz, uyz;

            for (int ix3 = minX3; ix3 <= maxX3; ix3++)
            {
               for (int ix2 = minX2; ix2 <= maxX2; ix2++)
               {
                  for (int ix1 = minX1; ix1 <= maxX1; ix1++)
                  {
                     if (!bcArray.isUndefined(ix1, ix2, ix3) && !bcArray.isSolid(ix1, ix2, ix3))
                     {
                        //////////////////////////////////////////////////////////////////////////
                        //compute average values
                        //////////////////////////////////////////////////////////////////////////

                        //mean velocity
                        if ((options&Velocity) == Velocity)
                        {
                           ux = (*av)(Vx, ix1, ix2, ix3) / timeSteps;
                           uy = (*av)(Vy, ix1, ix2, ix3) / timeSteps;
                           uz = (*av)(Vz, ix1, ix2, ix3) / timeSteps;

                           (*av)(Vx, ix1, ix2, ix3) = ux;
                           (*av)(Vy, ix1, ix2, ix3) = uy;
                           (*av)(Vz, ix1, ix2, ix3) = uz;
                        }

                        //fluctuations
                        if ((options&Fluctuations) == Fluctuations)
                        {
                           uxx = (*af)(Vxx, ix1, ix2, ix3) / timeSteps;
                           uyy = (*af)(Vyy, ix1, ix2, ix3) / timeSteps;
                           uzz = (*af)(Vzz, ix1, ix2, ix3) / timeSteps;
                           uxy = (*af)(Vxy, ix1, ix2, ix3) / timeSteps;
                           uxz = (*af)(Vxz, ix1, ix2, ix3) / timeSteps;
                           uyz = (*af)(Vyz, ix1, ix2, ix3) / timeSteps;

                           (*af)(Vxx, ix1, ix2, ix3) = uxx - ux*ux;
                           (*af)(Vyy, ix1, ix2, ix3) = uyy - uy*uy;
                           (*af)(Vzz, ix1, ix2, ix3) = uzz - uz*uz;
                           (*af)(Vxy, ix1, ix2, ix3) = uxy - ux*uy;
                           (*af)(Vxz, ix1, ix2, ix3) = uxz - ux*uz;
                           (*af)(Vyz, ix1, ix2, ix3) = uyz - uy*uz;
                        }

                        if ((options&Triplecorrelations) == Triplecorrelations)
                        {
                           //triple-correlations
                           (*at)(Vxxx, ix1, ix2, ix3) = (*at)(Vxxx, ix1, ix2, ix3) / timeSteps - 3 * uxx*ux + 2 * ux*ux*ux;
                           (*at)(Vxxy, ix1, ix2, ix3) = (*at)(Vxxy, ix1, ix2, ix3) / timeSteps - 2 * uxy*ux - uxx*uy + 2 * ux*ux*uy;
                           (*at)(Vxxz, ix1, ix2, ix3) = (*at)(Vxxz, ix1, ix2, ix3) / timeSteps - 2 * uxz*ux - uxx*uz + 2 * ux*ux*uz;
                           (*at)(Vyyy, ix1, ix2, ix3) = (*at)(Vyyy, ix1, ix2, ix3) / timeSteps - 3 * uyy*uy + 2 * uy*uy*uy;
                           (*at)(Vyyx, ix1, ix2, ix3) = (*at)(Vyyx, ix1, ix2, ix3) / timeSteps - 2 * uxy*uy - uyy*ux + 2 * uy*uy*ux;
                           (*at)(Vyyz, ix1, ix2, ix3) = (*at)(Vyyz, ix1, ix2, ix3) / timeSteps - 2 * uyz*uy - uyy*uz + 2 * uy*uy*uz;
                           (*at)(Vzzz, ix1, ix2, ix3) = (*at)(Vzzz, ix1, ix2, ix3) / timeSteps - 3 * uzz*uz + 2 * uz*uz*uz;
                           (*at)(Vzzx, ix1, ix2, ix3) = (*at)(Vzzx, ix1, ix2, ix3) / timeSteps - 2 * uxz*uz - uzz*ux + 2 * uz*uz*ux;
                           (*at)(Vzzy, ix1, ix2, ix3) = (*at)(Vzzy, ix1, ix2, ix3) / timeSteps - 2 * uyz*uz - uzz*uy + 2 * uz*uz*uy;
                           (*at)(Vxyz, ix1, ix2, ix3) = (*at)(Vxyz, ix1, ix2, ix3) / timeSteps - uxy*uz - uxz*uy - uyz*ux + 2 * ux*uy*uz;
                        }
                        //////////////////////////////////////////////////////////////////////////
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::calculateSubtotal()
{

   using namespace D3Q27System;

   //Funktionszeiger
   calcMacros = NULL;
   if (compressible)
   {
      calcMacros = &calcCompMacroscopicValues;
   }
   else
   {
      calcMacros = &calcIncompMacroscopicValues;
   }

   LBMReal f[27];
#ifdef _OPENMP
#pragma omp parallel private (f)
#endif
   {
      for (int level = minInitLevel; level <= maxInitLevel; level++)
      {
         int i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
         //BOOST_FOREACH(Block3DPtr block, blockVector[level])
         for (i = 0; i < blockVector[level].size(); i++)
         {
            Block3DPtr block = blockVector[level][i];
            if (block)
            {
               LBMKernel3DPtr kernel = block->getKernel();
               BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
               DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
               AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageVelocity();
               AverageValuesArray3DPtr af = kernel->getDataSet()->getAverageFluctuations();
               AverageValuesArray3DPtr at = kernel->getDataSet()->getAverageTriplecorrelations();

               int minX1 = iMinX1;
               int minX2 = iMinX2;
               int minX3 = iMinX3;

               int maxX1 = int(distributions->getNX1());
               int maxX2 = int(distributions->getNX2());
               int maxX3 = int(distributions->getNX3());

               maxX1 -= 2;
               maxX2 -= 2;
               maxX3 -= 2;

               for (int ix3 = minX3; ix3 <= maxX3; ix3++)
               {
                  for (int ix2 = minX2; ix2 <= maxX2; ix2++)
                  {
                     for (int ix1 = minX1; ix1 <= maxX1; ix1++)
                     {
                        if (!bcArray.isUndefined(ix1, ix2, ix3) && !bcArray.isSolid(ix1, ix2, ix3))
                        {
                           //////////////////////////////////////////////////////////////////////////
                           //read distribution
                           ////////////////////////////////////////////////////////////////////////////
                           distributions->getDistribution(f, ix1, ix2, ix3);
                           //////////////////////////////////////////////////////////////////////////
                           //compute velocity
                           //////////////////////////////////////////////////////////////////////////
                           LBMReal vx, vy, vz, rho;
                           calcMacros(f, rho, vx, vy, vz);
                           //double press = D3Q27System::calcPress(f, rho, vx, vy, vz);

                           //////////////////////////////////////////////////////////////////////////
                           //compute subtotals
                           //////////////////////////////////////////////////////////////////////////

                           //mean velocity
                           if ((options&Velocity) == Velocity)
                           {
                              (*av)(Vx, ix1, ix2, ix3) = (*av)(Vx, ix1, ix2, ix3) + vx;
                              (*av)(Vy, ix1, ix2, ix3) = (*av)(Vy, ix1, ix2, ix3) + vy;
                              (*av)(Vz, ix1, ix2, ix3) = (*av)(Vz, ix1, ix2, ix3) + vz;
                           }

                           //fluctuations
                           if ((options&Fluctuations) == Fluctuations)
                           {
                              (*af)(Vxx, ix1, ix2, ix3) = (*af)(Vxx, ix1, ix2, ix3) + vx*vx;
                              (*af)(Vyy, ix1, ix2, ix3) = (*af)(Vyy, ix1, ix2, ix3) + vy*vy;
                              (*af)(Vzz, ix1, ix2, ix3) = (*af)(Vzz, ix1, ix2, ix3) + vz*vz;
                              (*af)(Vxy, ix1, ix2, ix3) = (*af)(Vxy, ix1, ix2, ix3) + vx*vy;
                              (*af)(Vxz, ix1, ix2, ix3) = (*af)(Vxz, ix1, ix2, ix3) + vx*vz;
                              (*af)(Vyz, ix1, ix2, ix3) = (*af)(Vyz, ix1, ix2, ix3) + vy*vz;
                           }

                           //triple-correlations
                           if ((options&Triplecorrelations) == Triplecorrelations)
                           {
                              (*at)(Vxxx, ix1, ix2, ix3) = (*at)(Vxxx, ix1, ix2, ix3) + vx*vx*vx;
                              (*at)(Vxxy, ix1, ix2, ix3) = (*at)(Vxxy, ix1, ix2, ix3) + vx*vx*vy;
                              (*at)(Vxxz, ix1, ix2, ix3) = (*at)(Vxxz, ix1, ix2, ix3) + vx*vx*vz;
                              (*at)(Vyyy, ix1, ix2, ix3) = (*at)(Vyyy, ix1, ix2, ix3) + vy*vy*vy;
                              (*at)(Vyyx, ix1, ix2, ix3) = (*at)(Vyyx, ix1, ix2, ix3) + vy*vy*vx;
                              (*at)(Vyyz, ix1, ix2, ix3) = (*at)(Vyyz, ix1, ix2, ix3) + vy*vy*vz;
                              (*at)(Vzzz, ix1, ix2, ix3) = (*at)(Vzzz, ix1, ix2, ix3) + vz*vz*vz;
                              (*at)(Vzzx, ix1, ix2, ix3) = (*at)(Vzzx, ix1, ix2, ix3) + vz*vz*vx;
                              (*at)(Vzzy, ix1, ix2, ix3) = (*at)(Vzzy, ix1, ix2, ix3) + vz*vz*vy;
                              (*at)(Vxyz, ix1, ix2, ix3) = (*at)(Vxyz, ix1, ix2, ix3) + vx*vy*vz;
                           }
                           //////////////////////////////////////////////////////////////////////////
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesCoProcessor::volumeAverage(double step)
{
   int istep = int(step);
   string fname = path + "/tav/" + "tav" + UbSystem::toString(istep) + ".csv";

   std::ofstream ostr;
   ostr.open(fname.c_str(), std::ios_base::out);
   if (!ostr)
   {
      ostr.clear();
      string path = UbSystem::getPathFromString(fname);
      if (path.size() > 0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out); }
      if (!ostr) throw UbException(UB_EXARGS, "couldn't open file " + fname);
   }
   ostr << "z;Vx;Vy;Vz;Vxx;Vyy;Vzz;Vxy;Vxz;Vyz;Vxxx;Vxxy;Vxxz;Vyyy;Vyyx;Vyyz;Vzzz;Vzzx;Vzzy;Vxyz\n";

   CommunicatorPtr comm = Communicator::getInstance();
   int size = (int)levels.size();

   for (int i = 0; i < size; i++)
   {
      int level = levels[i];
      double dx = grid->getDeltaX(level);
      double start = levelCoords[i];
      double stop = levelCoords[i + 1] - dx*0.5;

      for (double j = start; j <= stop; j += dx)
      {
         D3Q27IntegrateValuesHelper intValHelp(grid, comm,
            bounds[0], bounds[1], j,
            bounds[3], bounds[4], j + dx, level);

         intValHelp.calculateAV2();

         double numberOfFluidsNodes = intValHelp.getNumberOfFluidsNodes();
         if (numberOfFluidsNodes > 0)
         {
            double Vx = intValHelp.getAVx() / numberOfFluidsNodes;
            double Vy = intValHelp.getAVy() / numberOfFluidsNodes;
            double Vz = intValHelp.getAVz() / numberOfFluidsNodes;

            double Vxx = intValHelp.getAVxx() / numberOfFluidsNodes;
            double Vyy = intValHelp.getAVyy() / numberOfFluidsNodes;
            double Vzz = intValHelp.getAVzz() / numberOfFluidsNodes;
            double Vxy = intValHelp.getAVxy() / numberOfFluidsNodes;
            double Vxz = intValHelp.getAVxz() / numberOfFluidsNodes;
            double Vyz = intValHelp.getAVyz() / numberOfFluidsNodes;

            double Vxxx = intValHelp.getAVxxx() / numberOfFluidsNodes;
            double Vxxy = intValHelp.getAVxxy() / numberOfFluidsNodes;
            double Vxxz = intValHelp.getAVxxz() / numberOfFluidsNodes;
            double Vyyy = intValHelp.getAVyyy() / numberOfFluidsNodes;
            double Vyyx = intValHelp.getAVyyx() / numberOfFluidsNodes;
            double Vyyz = intValHelp.getAVyyz() / numberOfFluidsNodes;
            double Vzzz = intValHelp.getAVzzz() / numberOfFluidsNodes;
            double Vzzx = intValHelp.getAVzzx() / numberOfFluidsNodes;
            double Vzzy = intValHelp.getAVzzy() / numberOfFluidsNodes;
            double Vxyz = intValHelp.getAVxyz() / numberOfFluidsNodes;


            ostr << j + 0.5*dx << ";" << Vx << ";" << Vy << ";" << Vz << ";";
            ostr << Vxx << ";" << Vyy << ";" << Vzz << ";" << Vxy << ";" << Vxz << ";" << Vyz << ";";
            ostr << Vxxx << ";" << Vxxy << ";" << Vxxz << ";" << Vyyy << ";" << Vyyx << ";" << Vyyz << ";" << Vzzz << ";" << Vzzx << ";" << Vzzy << ";" << Vxyz << "\n";
         }
      }
   }
   ostr.close();
}



