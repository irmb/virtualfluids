#include "IntegrateValuesHelper.h"


#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/CoordinateTransformation3D.h>
#include <vector>

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "DataSet3D.h"
#include "BCArray3D.h"

//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::IntegrateValuesHelper(SPtr<Grid3D> grid, SPtr<Communicator> comm,
   double minX1, double minX2,
   double minX3, double maxX1,
   double maxX2, double maxX3) :

   grid(grid),
   comm(comm),
   sVx1(0.0), sVx2(0.0), sVx3(0.0), sRho(0.0), sCellVolume(0.0),
   numberOfFluidsNodes(0),
   numberOfSolidNodes(0)
{
   boundingBox = GbCuboid3DPtr(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
   init(-1);
}
//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::IntegrateValuesHelper(SPtr<Grid3D> grid, SPtr<Communicator> comm,
   double minX1, double minX2,
   double minX3, double maxX1,
   double maxX2, double maxX3,
   int level) :

   grid(grid),
   comm(comm),
   sVx1(0.0), sVx2(0.0), sVx3(0.0), sRho(0.0), sCellVolume(0.0),
   numberOfFluidsNodes(0),
   numberOfSolidNodes(0)
{
   boundingBox = GbCuboid3DPtr(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
   init(level);
}
//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::~IntegrateValuesHelper()
= default;
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::init(int level)
{
   root = comm->isRoot();

   double orgX1, orgX2, orgX3;
   int gridRank = grid->getRank();
   int minInitLevel, maxInitLevel;
   if (level < 0)
   {
      minInitLevel = this->grid->getCoarsestInitializedLevel();
      maxInitLevel = this->grid->getFinestInitializedLevel();
   }
   else
   {
      minInitLevel = level;
      maxInitLevel = level;
   }

   double numSolids = 0.0;
   double numFluids = 0.0;
   for (int level_it = minInitLevel; level_it <= maxInitLevel; level_it++)
   {
      std::vector<SPtr<Block3D>> blockVector;
      grid->getBlocks(level_it, gridRank, blockVector);
      for(SPtr<Block3D> block : blockVector)
      {
         CalcNodes cn;
         cn.block = block;
         //Koords bestimmen
         UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

         orgX1 = val<1>(org);
         orgX2 = val<2>(org);
         orgX3 = val<3>(org);

         SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
         SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
         int ghostLayerWitdh = kernel->getGhostLayerWidth();
         SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
         double internX1, internX2, internX3;

         double         dx = grid->getDeltaX(block);
         UbTupleDouble3 orgDelta = grid->getNodeOffset(block);

         for (int ix3 = ghostLayerWitdh; ix3 < (int)distributions->getNX3() - ghostLayerWitdh; ix3++)
         {
            for (int ix2 = ghostLayerWitdh; ix2 < (int)distributions->getNX2() - ghostLayerWitdh; ix2++)
            {
               for (int ix1 = ghostLayerWitdh; ix1 < (int)distributions->getNX1() - ghostLayerWitdh; ix1++)
               {
                  internX1 = orgX1 - val<1>(orgDelta) + ix1 * dx;
                  internX2 = orgX2 - val<2>(orgDelta) + ix2 * dx;
                  internX3 = orgX3 - val<3>(orgDelta) + ix3 * dx;
                  if (boundingBox->isPointInGbObject3D(internX1, internX2, internX3))
                  {
                     if (!bcArray->isSolid(ix1, ix2, ix3) && !bcArray->isUndefined(ix1, ix2, ix3))
                     {
                        cn.nodes.emplace_back(ix1, ix2, ix3);
                        numFluids++;
                     }
                     else if (bcArray->isSolid(ix1, ix2, ix3))
                     {
                        numSolids++;
                     }
                  }
               }
            }
         }
         if (cn.nodes.size() > 0)
            cnodes.push_back(cn);
      }
   }
   std::vector<double> rvalues;
   std::vector<double> values;
   values.push_back(numSolids);
   values.push_back(numFluids);
   rvalues = comm->gather(values);

   if (root)
   {
      numberOfSolidNodes = 0.0;
      numberOfFluidsNodes = 0.0;
      int rsize = (int)rvalues.size();
      int vsize = (int)values.size();
      for (int i = 0; i < rsize; i += vsize)
      {
         numberOfSolidNodes += rvalues[i];
         numberOfFluidsNodes += rvalues[i + 1];
      }
   }

}
//////////////////////////////////////////////////////////////////////////
// calculation conventional rho, velocity and averaged data
void IntegrateValuesHelper::calculateAV()
{
   clearData();

   for(CalcNodes cn : cnodes)
   {
      SPtr<ILBMKernel> kernel = cn.block->getKernel();
      SPtr<AverageValuesArray3D> averagedValues = kernel->getDataSet()->getAverageValues();

      for(UbTupleInt3 node : cn.nodes)
      {
         double Avx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVx);
         double Avy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVy);
         double Avz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVz);

         double Avxx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxx);
         double Avyy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVyy);
         double Avzz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVzz);

         double Avxz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxz);
         sAvVx1 += abs(Avx);
         sAvVx2 += abs(Avy);
         sAvVx3 += abs(Avz);

         sTSx1 += sqrt(Avxx);
         sTSx2 += sqrt(Avyy);
         sTSx3 += sqrt(Avzz);

         sTSx1x3 += Avxz;
         numberOfFluidsNodes++;
      }
   }
   std::vector<double> values;
   std::vector<double> rvalues;
   values.push_back(sAvVx1);
   values.push_back(sAvVx2);
   values.push_back(sAvVx3);
   values.push_back(sTSx1);
   values.push_back(sTSx2);
   values.push_back(sTSx3);
   values.push_back(sTSx1x3);
   values.push_back(numberOfFluidsNodes);

   rvalues = comm->gather(values);
   if (root)
   {
      clearData();
      for (int i = 0; i < (int)rvalues.size(); i += 8)
      {
         sAvVx1 += rvalues[i];
         sAvVx2 += rvalues[i + 1];
         sAvVx3 += rvalues[i + 2];
         sTSx1 += rvalues[i + 3];
         sTSx2 += rvalues[i + 4];
         sTSx3 += rvalues[i + 5];
         sTSx1x3 += rvalues[i + 6];
         numberOfFluidsNodes += rvalues[i + 7];
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::calculateMQ()
{
   LBMReal f[D3Q27System::ENDF + 1];
   LBMReal vx1, vx2, vx3, rho;
   clearData();

   //Funktionszeiger
   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);

   CalcMacrosFct calcMacros = NULL;

   for(CalcNodes cn : cnodes)
   {
      SPtr<ILBMKernel> kernel = cn.block->getKernel();
      LBMReal dx = 1.0 / (LBMReal)(1 << cn.block->getLevel());
      LBMReal cellVolume = dx*dx*dx;

      if (kernel->getCompressible())
      {
         calcMacros = &D3Q27System::calcCompMacroscopicValues;
      }
      else
      {
         calcMacros = &D3Q27System::calcIncompMacroscopicValues;
      }

      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
      for(UbTupleInt3 node : cn.nodes)
      {
         distributions->getDistribution(f, val<1>(node), val<2>(node), val<3>(node));
         calcMacros(f, rho, vx1, vx2, vx3);
         sRho += rho*cellVolume;
         sVx1 += vx1*cellVolume;
         sVx2 += vx2*cellVolume;
         sVx3 += vx3*cellVolume;
         sCellVolume += cellVolume;
      }
   }
   std::vector<double> values(5);
   std::vector<double> rvalues;
   values[0] = sRho;
   values[1] = sVx1;
   values[2] = sVx2;
   values[3] = sVx3;
   values[4] = sCellVolume;

   rvalues = comm->gather(values);
   if (root)
   {
      clearData();
      int rsize = (int)rvalues.size();
      int vsize = (int)values.size();
      for (int i = 0; i < rsize; i += vsize)
      {
         sRho += rvalues[i];
         sVx1 += rvalues[i + 1];
         sVx2 += rvalues[i + 2];
         sVx3 += rvalues[i + 3];
         sCellVolume += rvalues[i + 4];
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void IntegrateValuesHelper::clearData()
{
   sRho = 0.0;
   sVx1 = 0.0;
   sVx2 = 0.0;
   sVx3 = 0.0;
   sCellVolume = 0.0;
   //sVm = 0.0;
   //sPress = 0.0;
   //numberOfFluidsNodes = 0.0;
   sAvVx1 = 0.0;
   sAvVx2 = 0.0;
   sAvVx3 = 0.0;
   sTSx1 = 0.0;
   sTSx2 = 0.0;
   sTSx3 = 0.0;
   sTSx1x3 = 0.0;
}
//////////////////////////////////////////////////////////////////////////
LBMReal IntegrateValuesHelper::getNumberOfFluidsNodes()
{
   return this->numberOfFluidsNodes;
}
//////////////////////////////////////////////////////////////////////////
LBMReal IntegrateValuesHelper::getNumberOfSolidNodes()
{
   return this->numberOfSolidNodes;
}
//////////////////////////////////////////////////////////////////////////
GbCuboid3DPtr IntegrateValuesHelper::getBoundingBox()
{
   return this->boundingBox;
}
//////////////////////////////////////////////////////////////////////////
std::vector<IntegrateValuesHelper::CalcNodes> IntegrateValuesHelper::getCNodes()
{
   return cnodes;
}
