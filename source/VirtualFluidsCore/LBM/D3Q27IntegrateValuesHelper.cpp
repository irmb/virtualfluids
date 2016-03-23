#include "D3Q27IntegrateValuesHelper.h"

#include <boost/foreach.hpp>
#include <numerics/geometry3d/GbCuboid3D.h>
#include <vector>

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"

using namespace std;
//////////////////////////////////////////////////////////////////////////
D3Q27IntegrateValuesHelper::D3Q27IntegrateValuesHelper(Grid3DPtr grid, CommunicatorPtr comm,
   double minX1, double minX2,
   double minX3, double maxX1,
   double maxX2, double maxX3) :

   grid(grid),
   comm(comm),
   sVx1(0.0), sVx2(0.0), sVx3(0.0), sRho(0.0), sCellVolume(0.0),
   numberOfFluidsNodes(0),
   numberOfSolidNodes(0)
{
   boundingBox =  GbCuboid3DPtr(new GbCuboid3D(minX1, minX2, minX3, maxX1, maxX2, maxX3));
   init();
}
//////////////////////////////////////////////////////////////////////////
D3Q27IntegrateValuesHelper::~D3Q27IntegrateValuesHelper()
{
}
//////////////////////////////////////////////////////////////////////////
void D3Q27IntegrateValuesHelper::init()
{
   double orgX1, orgX2, orgX3;
   int gridRank = grid->getRank();
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   double numSolids = 0.0;
   double numFluids = 0.0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
      {
         CalcNodes cn;
         cn.block = block;
         //Koords bestimmen
         UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

         orgX1 = val<1>(org);
         orgX2 = val<2>(org);
         orgX3 = val<3>(org);

         LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
         BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
         int ghostLayerWitdh = kernel->getGhostLayerWidth();
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
         double internX1, internX2, internX3;

         double         dx       = grid->getDeltaX(block);
         UbTupleDouble3 orgDelta = grid->getNodeOffset(block);

         for (int ix3=ghostLayerWitdh; ix3<(int)distributions->getNX3()-ghostLayerWitdh; ix3++)
         {
            for (int ix2=ghostLayerWitdh; ix2<(int)distributions->getNX2()-ghostLayerWitdh; ix2++)
            {
               for (int ix1=ghostLayerWitdh; ix1<(int)distributions->getNX1()-ghostLayerWitdh; ix1++)
               {
                  internX1 = orgX1 - val<1>(orgDelta) +ix1 * dx;
                  internX2 = orgX2 - val<2>(orgDelta) +ix2 * dx;
                  internX3 = orgX3 - val<3>(orgDelta) +ix3 * dx;
                  if (boundingBox->isPointInGbObject3D(internX1, internX2, internX3))
                  {
                     if (!bcArray.isSolid(ix1, ix2, ix3) && !bcArray.isUndefined(ix1, ix2, ix3))
                     {
                        cn.nodes.push_back(UbTupleInt3(ix1, ix2, ix3));
                        numFluids++;
                     }
                     else if (bcArray.isSolid(ix1, ix2, ix3))
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
   vector<double> rvalues;
   vector<double> values;
   values.push_back(numSolids);
   values.push_back(numFluids);
   rvalues = comm->gather(values);

   if (comm->getProcessID() == comm->getRoot())
   {
      numberOfSolidNodes = 0.0;
      numberOfFluidsNodes = 0.0;
      int rsize = (int)rvalues.size();
      int vsize = (int)values.size();
      for (int i = 0; i < rsize; i += vsize)
      {
         numberOfSolidNodes += rvalues[i];
         numberOfFluidsNodes += rvalues[i+1];
      }
   }

}
// calculation conventional rho, velocity and averaged data
void D3Q27IntegrateValuesHelper::calculateAV()
{
   clearData();
   //Funktionszeiger
   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros = NULL;

   BOOST_FOREACH(CalcNodes cn, cnodes)
   {
      LBMKernel3DPtr kernel = cn.block->getKernel();
      if (kernel->getCompressible())
      {
         calcMacros = &D3Q27System::calcCompMacroscopicValues;
      }
      else
      {
         calcMacros = &D3Q27System::calcIncompMacroscopicValues;
      }
      BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
      int ghostLayerWitdh = kernel->getGhostLayerWidth();
      AverageValuesArray3DPtr averagedValues = kernel->getDataSet()->getAverageValues();

      BOOST_FOREACH(UbTupleInt3 node, cn.nodes)
      {
         double Avx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVx);
         double Avy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVy);
         double Avz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVz);

         double Avxx = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxx);
         double Avyy = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVyy);
         double Avzz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVzz);

         double Avxz = (*averagedValues)(val<1>(node), val<2>(node), val<3>(node), AvVxz);
         sAvVx1 +=abs(Avx);
         sAvVx2 +=abs(Avy);
         sAvVx3 +=abs(Avz);

         sTSx1 += sqrt(Avxx);
         sTSx2 += sqrt(Avyy);
         sTSx3 += sqrt(Avzz);

         sTSx1x3 += Avxz;
         numberOfFluidsNodes++;
      }
   }
   vector<double> values;
   vector<double> rvalues;
   values.push_back(sAvVx1);
   values.push_back(sAvVx2);
   values.push_back(sAvVx3);
   values.push_back(sTSx1);
   values.push_back(sTSx2);
   values.push_back(sTSx3);
   values.push_back(sTSx1x3);
   values.push_back(numberOfFluidsNodes);

   rvalues = comm->gather(values);
   if (comm->getProcessID() == comm->getRoot())
   {
      clearData();
      for (int i = 0; i < (int)rvalues.size(); i+=8)
      {
         sAvVx1 += rvalues[i];
         sAvVx2 += rvalues[i+1];
         sAvVx3 += rvalues[i+2];
         sTSx1 += rvalues[i+3];
         sTSx2 += rvalues[i+4];
         sTSx3 += rvalues[i+5];
         sTSx1x3 += rvalues[i+6];
         numberOfFluidsNodes += rvalues[i+7];
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27IntegrateValuesHelper::calculateMQ()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal vx1, vx2, vx3, rho;
   clearData();

   //Funktionszeiger
   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);

   CalcMacrosFct calcMacros = NULL;

   BOOST_FOREACH(CalcNodes cn, cnodes)
   {
      LBMKernel3DPtr kernel = cn.block->getKernel();
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

      BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
      int ghostLayerWitdh = kernel->getGhostLayerWidth();
      DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
      BOOST_FOREACH(UbTupleInt3 node, cn.nodes)
      {
         distributions->getDistribution(f, val<1>(node), val<2>(node), val<3>(node));
         calcMacros(f, rho, vx1, vx2, vx3);
         //press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);
         sRho += rho*cellVolume;
         sVx1 += vx1*cellVolume;
         sVx2 += vx2*cellVolume;
         sVx3 += vx3*cellVolume;
         sCellVolume+=cellVolume;
         //sPress += press*area;
         //sVm += (sqrt(vx1*vx1 + vx2*vx2 + vx3*vx3)*area);
      }
   }
   vector<double> values;
   vector<double> rvalues;
   values.push_back(sRho);
   values.push_back(sVx1);
   values.push_back(sVx2);
   values.push_back(sVx3);
   values.push_back(sCellVolume);
   //values.push_back(sPress);
   //values.push_back(sVm);

   comm->allGather(values, rvalues);
   clearData();
   int rsize = (int)rvalues.size();
   int vsize = (int)values.size();
   for (int i = 0; i < rsize; i+=vsize)
   {
      sRho += rvalues[i];
      sVx1 += rvalues[i+1];
      sVx2 += rvalues[i+2];
      sVx3 += rvalues[i+3];
      sCellVolume += rvalues[i+4];
      //sPress += rvalues[i+5];
      //sVm += rvalues[i+6];
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27IntegrateValuesHelper::clearData()
{
   sRho = 0.0;
   sVx1 = 0.0;
   sVx2 = 0.0;
   sVx3 = 0.0;
   sCellVolume= 0.0;
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
LBMReal D3Q27IntegrateValuesHelper::getNumberOfFluidsNodes()
{
   return this->numberOfFluidsNodes;
}
//////////////////////////////////////////////////////////////////////////
LBMReal D3Q27IntegrateValuesHelper::getNumberOfSolidNodes()
{
   return this->numberOfSolidNodes;
}
//////////////////////////////////////////////////////////////////////////
GbCuboid3DPtr D3Q27IntegrateValuesHelper::getBoundingBox()
{
   return this->boundingBox;
}
//////////////////////////////////////////////////////////////////////////
std::vector<CalcNodes> D3Q27IntegrateValuesHelper::getCNodes()
{
   return cnodes;
}
