#include "IntegrateValuesHelper.h"


#include <numerics/geometry3d/GbCuboid3D.h>
#include <numerics/geometry3d/CoordinateTransformation3D.h>
#include <vector>

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "DataSet3D.h"
#include "BCArray3D.h"

//////////////////////////////////////////////////////////////////////////
IntegrateValuesHelper::IntegrateValuesHelper(Grid3DPtr grid, CommunicatorPtr comm,
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
IntegrateValuesHelper::IntegrateValuesHelper(Grid3DPtr grid, CommunicatorPtr comm,
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
{
}
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
   for (int level = minInitLevel; level <= maxInitLevel; level++)
   {
      std::vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      for(Block3DPtr block : blockVector)
      {
         CalcNodes cn;
         cn.block = block;
         //Koords bestimmen
         UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

         orgX1 = val<1>(org);
         orgX2 = val<2>(org);
         orgX3 = val<3>(org);

         LBMKernelPtr kernel = std::dynamic_pointer_cast<LBMKernel>(block->getKernel());
         BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();
         int ghostLayerWitdh = kernel->getGhostLayerWidth();
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
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
                        cn.nodes.push_back(UbTupleInt3(ix1, ix2, ix3));
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
void IntegrateValuesHelper::prepare2DMatrix(int level)
{
   root = comm->isRoot();

   double orgX1, orgX2, orgX3;
   int gridRank = grid->getRank();
   int minInitLevel, maxInitLevel;
   if (level<0)
   {
      minInitLevel = this->grid->getCoarsestInitializedLevel();
      maxInitLevel = this->grid->getFinestInitializedLevel();
   }
   else
   {
      minInitLevel = level;
      maxInitLevel = level;
   }
   double dx = grid->getDeltaX(level);
   CoordinateTransformation3D trafo(boundingBox->getX1Minimum(),boundingBox->getX2Minimum(),boundingBox->getX3Minimum(),dx,dx,dx);
   cnodes2DMatrix.resize(UbMath::integerRounding<double>(boundingBox->getX1Maximum()/dx),UbMath::integerRounding<double>(boundingBox->getX2Maximum()/dx));

   double numSolids = 0.0;
   double numFluids = 0.0;
   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      std::vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      for(Block3DPtr block : blockVector)
      {
         Node cn;
         cn.block = block;
         //Koords bestimmen
         UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

         orgX1 = val<1>(org);
         orgX2 = val<2>(org);
         orgX3 = val<3>(org);

         ILBMKernelPtr kernel = block->getKernel();
         BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();
         int ghostLayerWitdh = kernel->getGhostLayerWidth();
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
         double internX1, internX2, internX3;

         double         dx = grid->getDeltaX(block);
         UbTupleDouble3 orgDelta = grid->getNodeOffset(block);

         for (int ix3 = ghostLayerWitdh; ix3<(int)distributions->getNX3()-ghostLayerWitdh; ix3++)
         {
            for (int ix2 = ghostLayerWitdh; ix2<(int)distributions->getNX2()-ghostLayerWitdh; ix2++)
            {
               for (int ix1 = ghostLayerWitdh; ix1<(int)distributions->getNX1()-ghostLayerWitdh; ix1++)
               {
                  internX1 = orgX1-val<1>(orgDelta)+ix1 * dx;
                  internX2 = orgX2-val<2>(orgDelta)+ix2 * dx;
                  internX3 = orgX3-val<3>(orgDelta)+ix3 * dx;
                  if (boundingBox->isPointInGbObject3D(internX1, internX2, internX3))
                  {
                     if (!bcArray->isSolid(ix1, ix2, ix3)&&!bcArray->isUndefined(ix1, ix2, ix3))
                     {
                        cn.node = UbTupleInt3(ix1, ix2, ix3);
                        int x1 = (int)trafo.transformForwardToX1Coordinate(internX1, internX2, internX3);
                        int x2 = (int)trafo.transformForwardToX2Coordinate(internX1, internX2, internX3);
                        cnodes2DMatrix(x1,x2)=cn;
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
      for (int i = 0; i<rsize; i += vsize)
      {
         numberOfSolidNodes += rvalues[i];
         numberOfFluidsNodes += rvalues[i+1];
      }
   }

}
// calculation conventional rho, velocity and averaged data
void IntegrateValuesHelper::calculateAV()
{
   clearData();

   for(CalcNodes cn : cnodes)
   {
      ILBMKernelPtr kernel = cn.block->getKernel();
      AverageValuesArray3DPtr averagedValues = kernel->getDataSet()->getAverageValues();

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
// calculation conventional rho, velocity and averaged data
void IntegrateValuesHelper::calculateAV2()
{
   saVx = 0;
   saVy = 0;
   saVz = 0;

   saVxx = 0;
   saVyy = 0;
   saVzz = 0;
   saVxy = 0;
   saVxz = 0;
   saVyz = 0;

   saVxxx = 0;
   saVxxy = 0;
   saVxxz = 0;
   saVyyy = 0;
   saVyyx = 0;
   saVyyz = 0;
   saVzzz = 0;
   saVzzx = 0;
   saVzzy = 0;
   saVxyz = 0;

   double lsaVx = 0;
   double lsaVy = 0;
   double lsaVz = 0;

   double lsaVxx = 0;
   double lsaVyy = 0;
   double lsaVzz = 0;
   double lsaVxy = 0;
   double lsaVxz = 0;
   double lsaVyz = 0;

   double lsaVxxx = 0;
   double lsaVxxy = 0;
   double lsaVxxz = 0;
   double lsaVyyy = 0;
   double lsaVyyx = 0;
   double lsaVyyz = 0;
   double lsaVzzz = 0;
   double lsaVzzx = 0;
   double lsaVzzy = 0;
   double lsaVxyz = 0;

   for(CalcNodes cn : cnodes)
   {
      ILBMKernelPtr kernel = cn.block->getKernel();
      AverageValuesArray3DPtr averagedVelocity = kernel->getDataSet()->getAverageVelocity();
      AverageValuesArray3DPtr averagedFluctuations = kernel->getDataSet()->getAverageFluctuations();
      AverageValuesArray3DPtr averagedTriplecorrelations = kernel->getDataSet()->getAverageTriplecorrelations();

      for(UbTupleInt3 node : cn.nodes)
      {
         double aVx = (*averagedVelocity)(Vx, val<1>(node), val<2>(node), val<3>(node));
         double aVy = (*averagedVelocity)(Vy, val<1>(node), val<2>(node), val<3>(node));
         double aVz = (*averagedVelocity)(Vz, val<1>(node), val<2>(node), val<3>(node));

         double aVxx = (*averagedFluctuations)(Vxx, val<1>(node), val<2>(node), val<3>(node));
         double aVyy = (*averagedFluctuations)(Vyy, val<1>(node), val<2>(node), val<3>(node));
         double aVzz = (*averagedFluctuations)(Vzz, val<1>(node), val<2>(node), val<3>(node));
         double aVxy = (*averagedFluctuations)(Vxy, val<1>(node), val<2>(node), val<3>(node));
         double aVxz = (*averagedFluctuations)(Vxz, val<1>(node), val<2>(node), val<3>(node));
         double aVyz = (*averagedFluctuations)(Vyz, val<1>(node), val<2>(node), val<3>(node));

         double aVxxx = (*averagedTriplecorrelations)(Vxxx, val<1>(node), val<2>(node), val<3>(node));
         double aVxxy = (*averagedTriplecorrelations)(Vxxy, val<1>(node), val<2>(node), val<3>(node));
         double aVxxz = (*averagedTriplecorrelations)(Vxxz, val<1>(node), val<2>(node), val<3>(node));
         double aVyyy = (*averagedTriplecorrelations)(Vyyy, val<1>(node), val<2>(node), val<3>(node));
         double aVyyx = (*averagedTriplecorrelations)(Vyyx, val<1>(node), val<2>(node), val<3>(node));
         double aVyyz = (*averagedTriplecorrelations)(Vyyz, val<1>(node), val<2>(node), val<3>(node));
         double aVzzz = (*averagedTriplecorrelations)(Vzzz, val<1>(node), val<2>(node), val<3>(node));
         double aVzzx = (*averagedTriplecorrelations)(Vzzx, val<1>(node), val<2>(node), val<3>(node));
         double aVzzy = (*averagedTriplecorrelations)(Vzzy, val<1>(node), val<2>(node), val<3>(node));
         double aVxyz = (*averagedTriplecorrelations)(Vxyz, val<1>(node), val<2>(node), val<3>(node));

         lsaVx += aVx;
         lsaVy += aVy;
         lsaVz += aVz;

         lsaVxx += aVxx;
         lsaVyy += aVyy;
         lsaVzz += aVzz;
         lsaVxy += aVxy;
         lsaVxz += aVxz;
         lsaVyz += aVyz;

         lsaVxxx += aVxxx;
         lsaVxxy += aVxxy;
         lsaVxxz += aVxxz;
         lsaVyyy += aVyyy;
         lsaVyyx += aVyyx;
         lsaVyyz += aVyyz;
         lsaVzzz += aVzzz;
         lsaVzzx += aVzzx;
         lsaVzzy += aVzzy;
         lsaVxyz += aVxyz;

         //numberOfFluidsNodes++;
      }
   }
   std::vector<double> values;
   std::vector<double> rvalues;
   values.push_back(lsaVx);
   values.push_back(lsaVy);
   values.push_back(lsaVz);

   values.push_back(lsaVxx);
   values.push_back(lsaVyy);
   values.push_back(lsaVzz);
   values.push_back(lsaVxy);
   values.push_back(lsaVxz);
   values.push_back(lsaVyz);

   values.push_back(lsaVxxx);
   values.push_back(lsaVxxy);
   values.push_back(lsaVxxz);
   values.push_back(lsaVyyy);
   values.push_back(lsaVyyx);
   values.push_back(lsaVyyz);
   values.push_back(lsaVzzz);
   values.push_back(lsaVzzx);
   values.push_back(lsaVzzy);
   values.push_back(lsaVxyz);

   rvalues = comm->gather(values);
   if (root)
   {
      for (int i = 0; i < (int)rvalues.size(); i += 19)
      {
         saVx += rvalues[i];
         saVy += rvalues[i + 1];
         saVz += rvalues[i + 2];

         saVxx += rvalues[i + 3];
         saVyy += rvalues[i + 4];
         saVzz += rvalues[i + 5];
         saVxy += rvalues[i + 6];
         saVxz += rvalues[i + 7];
         saVyz += rvalues[i + 8];

         saVxxx += rvalues[i + 9];
         saVxxy += rvalues[i + 10];
         saVxxz += rvalues[i + 11];
         saVyyy += rvalues[i + 12];
         saVyyx += rvalues[i + 13];
         saVyyz += rvalues[i + 14];
         saVzzz += rvalues[i + 15];
         saVzzx += rvalues[i + 16];
         saVzzy += rvalues[i + 17];
         saVxyz += rvalues[i + 18];

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
      ILBMKernelPtr kernel = cn.block->getKernel();
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

      BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();
      DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
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
