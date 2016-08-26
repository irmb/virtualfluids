//#include "InitDistributionsWithVelocityProfileBlockVisitor.h"
//#include <basics/utilities/UbFileInputASCII.h>
//#include "LBMKernel3D.h"
//#include "D3Q27ETBCProcessor.h"
//#include "Grid3DSystem.h"
//
//using namespace std;
//
//InitDistributionsWithVelocityProfileBlockVisitor::InitDistributionsWithVelocityProfileBlockVisitor(LBMReal nu, LBMReal rho, double dx, int size, std::string filename)
//   : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), nu(nu), rho(rho)
//{
//   UbFileInputASCII in(filename);
//   if (!in)
//   {
//      throw UbException(UB_EXARGS, "could not open file " + filename);
//   }
//   
//   in.readLine();
//
//   while(!in.eof())
//   {
//      z.push_back(in.readDouble());
//      vx1.push_back(in.readDouble());
//      vx2.push_back(in.readDouble());
//      vx3.push_back(in.readDouble());
//   }
//
//   double startCoord = z[0];
//   int old_size = (int)z.size();
//   for (int i = 0; i < old_size; i++)
//   {
//      z[i] = (z[i] - startCoord) / dx;
//   }
//
//
//}
////////////////////////////////////////////////////////////////////////////
//InitDistributionsWithVelocityProfileBlockVisitor::~InitDistributionsWithVelocityProfileBlockVisitor()
//{
//}
////////////////////////////////////////////////////////////////////////////
//void InitDistributionsWithVelocityProfileBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
//{
//   using namespace D3Q27System;
//
//   int gridRank = oldGrid->getRank();
//
//   int minInitLevel = oldGrid->getCoarsestInitializedLevel();
//   int maxInitLevel = oldGrid->getFinestInitializedLevel();
//
//   for (int level = minInitLevel; level <= maxInitLevel; level++)
//   {
//      vector<Block3DPtr> blockVector;
//      grid->getBlocks(level, gridRank, true, blockVector);
//      BOOST_FOREACH(Block3DPtr block, blockVector)
//      {
//         if (block)
//         {
//            UbTupleDouble3 coords = oldGrid->getBlockWorldCoordinates(block);
//
//         }
//      }
//   }
//
//   grid->getBlockWorldCoordinates(0,0,0,)
//
//   int size = z.size();
//   for (int i = 0; i < size; i++)
//   {
//      for (int j; j < )
//      {
//      }
//      grid->getBlockIndexes(z[i], grid->, 0);
//   }
//
//   
//   
//   
//   if (!block) UB_THROW(UbException(UB_EXARGS, "block is not exist"));
//
//   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
//   UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
//   double dx = grid->getDeltaX(block);
//   LBMReal o = LBMSystem::calcCollisionFactor(nu, block->getLevel());
//
//
//   //Funktionszeiger
//   typedef void(*CalcFeqsFct)(LBMReal* const& /*feq[27]*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
//   CalcFeqsFct   calcFeqsFct = NULL;
//
//   LBMReal vx1, vx2, vx3, rho;
//
//   int gridRank = grid->getRank();
//   int blockRank = block->getRank();
//
//   if (blockRank == gridRank && block->isActive())
//   {
//      LBMKernel3DPtr kernel = block->getKernel();
//      if (!kernel)
//         throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: " + block->toString());
//
//      if (kernel->getCompressible())
//         calcFeqsFct = &D3Q27System::calcCompFeq;
//      else
//         calcFeqsFct = &D3Q27System::calcIncompFeq;
//
//      UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
//
//      BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
//      EsoTwist3DPtr           distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
//
//      LBMReal f[D3Q27System::ENDF + 1];
//
//      size_t nx1 = distributions->getNX1();
//      size_t nx2 = distributions->getNX2();
//      size_t nx3 = distributions->getNX3();
//
//      int minX1 = 0;
//      int minX2 = 0;
//      int minX3 = 0;
//
//      int maxX1 = (int)bcArray.getNX1();
//      int maxX2 = (int)bcArray.getNX2();
//      int maxX3 = (int)bcArray.getNX3();
//
//      int level = block->getLevel();
//      int lMax = grid->getFinestInitializedLevel();
//
//      for (int ix3 = minX1; ix3 < maxX3; ix3++)
//         for (int ix2 = minX2; ix2 < maxX2; ix2++)
//            for (int ix1 = minX1; ix1 < maxX1; ix1++)
//            {
//               int x3f = ix3 + pow(2, lMax - level);
//
//
//
//               UbTupleInt3 coords = grid->getNodeIndexes(block, ix1, ix2, ix3);
//               int x1 = val<1>(coords);
//               int x2 = val<2>(coords);
//               int x3 = val<3>(coords);
//
//               vx1 = matrix(Vx1, x1, x2, x3);
//               vx2 = matrix(Vx2, x1, x2, x3);
//               vx3 = matrix(Vx3, x1, x2, x3);
//
//
//               //x-derivative
//               //double deltaX = dx*0.5;
//               int deltaX = 1;
//               x1 = val<1>(coords) + deltaX;
//               if (x1 > maxX1) x1 = val<1>(coords);
//               double vx1Plusx1 = matrix(Vx1, x1, x2, x3);
//               double vx2Plusx1 = matrix(Vx2, x1, x2, x3);
//               double vx3Plusx1 = matrix(Vx3, x1, x2, x3);
//
//               x1 = val<1>(coords) - deltaX;
//               if (x1 < minX1) x1 = val<1>(coords);
//               double vx1Minusx1 = matrix(Vx1, x1, x2, x3);
//               double vx2Minusx1 = matrix(Vx2, x1, x2, x3);
//               double vx3Minusx1 = matrix(Vx3, x1, x2, x3);
//
//               //y-derivative
//               x1 = val<1>(coords);
//               x2 = val<2>(coords) + deltaX;
//               if (x2 > maxX2) x2 = val<2>(coords);
//               double vx1Plusx2 = matrix(Vx1, x1, x2, x3);
//               double vx2Plusx2 = matrix(Vx2, x1, x2, x3);
//               double vx3Plusx2 = matrix(Vx3, x1, x2, x3);
//
//               x2 = val<2>(coords) - deltaX;
//               if (x2 < minX2) x2 = val<2>(coords);
//               double vx1Minusx2 = matrix(Vx1, x1, x2, x3);
//               double vx2Minusx2 = matrix(Vx2, x1, x2, x3);
//               double vx3Minusx2 = matrix(Vx3, x1, x2, x3);
//
//               //z-derivative
//               x2 = val<2>(coords);
//               x3 = val<3>(coords) + deltaX;
//               if (x3 > maxX3) x3 = val<3>(coords);
//               double vx1Plusx3 = matrix(Vx1, x1, x2, x3);
//               double vx2Plusx3 = matrix(Vx2, x1, x2, x3);
//               double vx3Plusx3 = matrix(Vx3, x1, x2, x3);
//
//               x3 = val<3>(coords) - deltaX;
//               if (x3 < minX3) x3 = val<3>(coords);
//               double vx1Minusx3 = matrix(Vx1, x1, x2, x3);
//               double vx2Minusx3 = matrix(Vx2, x1, x2, x3);
//               double vx3Minusx3 = matrix(Vx3, x1, x2, x3);
//
//               double ax = (vx1Plusx1 - vx1Minusx1) / (2.0*deltaX);
//               double bx = (vx2Plusx1 - vx2Minusx1) / (2.0*deltaX);
//               double cx = (vx3Plusx1 - vx3Minusx1) / (2.0*deltaX);
//
//               double ay = (vx1Plusx2 - vx1Minusx2) / (2.0*deltaX);
//               double by = (vx2Plusx2 - vx2Minusx2) / (2.0*deltaX);
//               double cy = (vx3Plusx2 - vx3Minusx2) / (2.0*deltaX);
//
//               double az = (vx1Plusx3 - vx1Minusx3) / (2.0*deltaX);
//               double bz = (vx2Plusx3 - vx2Minusx3) / (2.0*deltaX);
//               double cz = (vx3Plusx3 - vx3Minusx3) / (2.0*deltaX);
//               double eps_new = 1.0;
//               LBMReal op = 1.;
//
//               LBMReal feq[27];
//
//               calcFeqsFct(feq, rho, vx1, vx2, vx3);
//
//               double f_E = eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op) / (54.*o*op));
//               double f_N = f_E + eps_new *((2.*(ax - by)) / (9.*o));
//               double f_T = f_E + eps_new *((2.*(ax - cz)) / (9.*o));
//               double f_NE = eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op)) / (54.*o*op));
//               double f_SE = f_NE + eps_new *((ay + bx) / (9.*o));
//               double f_TE = eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op)) / (54.*o*op));
//               double f_BE = f_TE + eps_new *((az + cx) / (9.*o));
//               double f_TN = eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op) / (54.*o*op));
//               double f_BN = f_TN + eps_new *((bz + cy) / (9.*o));
//               double f_ZERO = eps_new *((5.*(ax + by + cz)) / (9.*op));
//               double f_TNE = eps_new *(-(ay + az + bx + bz + cx + cy) / (72.*o));
//               double f_TSW = -eps_new *((ay + bx) / (36.*o)) - f_TNE;
//               double f_TSE = -eps_new *((az + cx) / (36.*o)) - f_TNE;
//               double f_TNW = -eps_new *((bz + cy) / (36.*o)) - f_TNE;
//
//
//               f[E] = f_E + feq[E];
//               f[W] = f_E + feq[W];
//               f[N] = f_N + feq[N];
//               f[S] = f_N + feq[S];
//               f[T] = f_T + feq[T];
//               f[B] = f_T + feq[B];
//               f[NE] = f_NE + feq[NE];
//               f[SW] = f_NE + feq[SW];
//               f[SE] = f_SE + feq[SE];
//               f[NW] = f_SE + feq[NW];
//               f[TE] = f_TE + feq[TE];
//               f[BW] = f_TE + feq[BW];
//               f[BE] = f_BE + feq[BE];
//               f[TW] = f_BE + feq[TW];
//               f[TN] = f_TN + feq[TN];
//               f[BS] = f_TN + feq[BS];
//               f[BN] = f_BN + feq[BN];
//               f[TS] = f_BN + feq[TS];
//               f[TNE] = f_TNE + feq[TNE];
//               f[TNW] = f_TNW + feq[TNW];
//               f[TSE] = f_TSE + feq[TSE];
//               f[TSW] = f_TSW + feq[TSW];
//               f[BNE] = f_TSW + feq[BNE];
//               f[BNW] = f_TSE + feq[BNW];
//               f[BSE] = f_TNW + feq[BSE];
//               f[BSW] = f_TNE + feq[BSW];
//               f[ZERO] = f_ZERO + feq[ZERO];
//
//               distributions->setDistribution(f, ix1, ix2, ix3);
//               distributions->setDistributionInv(f, ix1, ix2, ix3);
//            }
//   }
//}
////////////////////////////////////////////////////////////////////////////
