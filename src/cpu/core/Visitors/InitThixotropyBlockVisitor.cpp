 //=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file InitThixotropyBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#include "InitThixotropyBlockVisitor.h"
#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "BCArray3D.h"

InitThixotropyBlockVisitor::InitThixotropyBlockVisitor()
   : Block3DVisitor(0, D3Q27System::MAXLEVEL)
{
    using namespace vf::basics::constant;

   //this->setVx1(0.0);
   //this->setVx2(0.0);
   //this->setVx3(0.0);
   //this->setRho(0.0);
   //this->setf1(0.0);
   //this->setf2(0.0);
   //this->setf3(0.0);
   //this->setConcentration(0.0);
   this->setLambda(c0o1);
}
//////////////////////////////////////////////////////////////////////////
//InitThixotropyBlockVisitor::InitThixotropyBlockVisitor(real lambda /*real nu, real D, real rho, real vx1, real vx2, real vx3, real c, real f1, real f2, real f3*/)
//    : Block3DVisitor(0, D3Q27System::MAXLEVEL)
//{
//    //this->setVx1(vx1);
//    //this->setVx2(vx2);
//    //this->setVx3(vx3);
//    //this->setRho(rho);
//    //this->setf1(vx1);
//    //this->setf2(vx2);
//    //this->setf3(vx3);
//    //this->setConcentration(rho);
//    //this->setNu(nu);
//    //this->setD(D);
//    this->setLambda(lambda);
//}
//////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx1(const mu::Parser& parser)
//{
//   this->checkFunction(parser);
//   this->muVx1 = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx2(const mu::Parser& parser)
//{
//   this->checkFunction(parser);
//   this->muVx2 = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx3(const mu::Parser& parser)
//{
//   this->checkFunction(parser);
//   this->muVx3 = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setRho(const mu::Parser& parser)
//{
//    this->checkFunction(parser);
//    this->muRho = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf1(const mu::Parser& parser)
//{
//    this->checkFunction(parser);
//    this->muf1 = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf2(const mu::Parser& parser)
//{
//    this->checkFunction(parser);
//    this->muf2 = parser;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf3(const mu::Parser& parser)
//{
//    this->checkFunction(parser);
//    this->muf3 = parser;
//}
////////////////////////////////////////////////////////////////////////////
void InitThixotropyBlockVisitor::setLambda(const mu::Parser& parser)
{
   this->checkFunction(parser);
   this->muLambda = parser;
}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx1(const std::string& muParserString)
//{
//    this->muVx1.SetExpr(muParserString);
//    this->checkFunction(muVx1);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx2(const std::string& muParserString)
//{
//    this->muVx2.SetExpr(muParserString);
//    this->checkFunction(muVx2);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx3(const std::string& muParserString)
//{
//    this->muVx3.SetExpr(muParserString);
//    this->checkFunction(muVx3);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setRho(const std::string& muParserString)
//{
//    this->muRho.SetExpr(muParserString);
//    this->checkFunction(muRho);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf1(const std::string& muParserString)
//{
//    this->muf1.SetExpr(muParserString);
//    this->checkFunction(muf1);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf2(const std::string& muParserString)
//{
//    this->muf2.SetExpr(muParserString);
//    this->checkFunction(muf2);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf3(const std::string& muParserString)
//{
//    this->muf3.SetExpr(muParserString);
//    this->checkFunction(muf3);
//}
////////////////////////////////////////////////////////////////////////////
void InitThixotropyBlockVisitor::setLambda(const std::string& muParserString)
{
   this->muLambda.SetExpr(muParserString);
   this->checkFunction(muLambda);
}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx1(real vx1)
//{
//    this->muVx1.SetExpr(UbSystem::toString(vx1, realLim::digits10));
//    this->checkFunction(muVx1);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx2(real vx2)
//{
//    this->muVx2.SetExpr(UbSystem::toString(vx2, realLim::digits10));
//    this->checkFunction(muVx2);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setVx3(real vx3)
//{
//    this->muVx3.SetExpr(UbSystem::toString(vx3, realLim::digits10));
//    this->checkFunction(muVx3);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setRho(real rho)
//{
//    this->muRho.SetExpr(UbSystem::toString(rho, realLim::digits10));
//    this->checkFunction(muRho);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf1(real f1)
//{
//    this->muf1.SetExpr(UbSystem::toString(f1, realLim::digits10));
//    this->checkFunction(muf1);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf2(real f2)
//{
//    this->muf2.SetExpr(UbSystem::toString(f2, realLim::digits10));
//    this->checkFunction(muf2);
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setf3(real f3)
//{
//    this->muf3.SetExpr(UbSystem::toString(f3, realLim::digits10));
//    this->checkFunction(muf3);
//}
//////////////////////////////////////////////////////////////////////////
void InitThixotropyBlockVisitor::setLambda(real lambda)
{
   this->muLambda.SetExpr(UbSystem::toString(lambda, realLim::digits10));
   this->checkFunction(muLambda);
}
//////////////////////////////////////////////////////////////////////////
void InitThixotropyBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   using namespace D3Q27System;
   using namespace vf::basics::constant;

   if(!block) UB_THROW( UbException(UB_EXARGS,"block is not exist") );

//   double dx = grid->getDeltaX(block);

   //define vars for functions
   mu::value_type x1,x2,x3;
   //this->muVx1.DefineVar("x1",&x1); this->muVx1.DefineVar("x2",&x2); this->muVx1.DefineVar("x3",&x3);
   //this->muVx2.DefineVar("x1",&x1); this->muVx2.DefineVar("x2",&x2); this->muVx2.DefineVar("x3",&x3);
   //this->muVx3.DefineVar("x1",&x1); this->muVx3.DefineVar("x2",&x2); this->muVx3.DefineVar("x3",&x3);
   //this->muRho.DefineVar("x1",&x1); this->muRho.DefineVar("x2",&x2); this->muRho.DefineVar("x3",&x3);

   this->muLambda.DefineVar("x1",&x1); this->muLambda.DefineVar("x2",&x2); this->muLambda.DefineVar("x3",&x3);

   //Funktionszeiger
   typedef void (*CalcFeqsFct)(real* const& /*feq[27]*/,const real& /*(d)rho*/,const real& /*vx1*/,const real& /*vx2*/,const real& /*vx3*/);
   CalcFeqsFct   calcFeqsFct   = NULL;

   int gridRank = grid->getRank();
   int blockRank = block->getRank();

   if (blockRank == gridRank && block->isActive())
   {
       SPtr<ILBMKernel> kernel = block->getKernel();
      if (!kernel)
         throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

      if(kernel->getCompressible()) 
         calcFeqsFct   = &D3Q27System::calcCompFeq; 
      else                                                        
         calcFeqsFct   = &D3Q27System::calcIncompFeq; 

      SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getHdistributions();  

      real h[D3Q27System::ENDF+1];

      for(std::size_t ix3=0; ix3<bcArray->getNX3(); ix3++)
         for(std::size_t ix2=0; ix2<bcArray->getNX2(); ix2++)
            for(std::size_t ix1=0; ix1<bcArray->getNX1(); ix1++)
            {
               //UbTupleDouble3 coords = grid->getNodeCoordinates(block, ix1, ix2, ix3);
               //x1 = val<1>(coords);
               //x2 = val<2>(coords);
               //x3 = val<3>(coords);

               //vx1 = muVx1.Eval();
               //vx2 = muVx2.Eval();
               //vx3 = muVx3.Eval();
               //rho = muRho.Eval();

               //f1 = muf1.Eval();
               //f2 = muf2.Eval();
               //f3 = muf3.Eval();
               //conc = muConcentration.Eval();

               //initialize(f, x1, x2, x3, vx1, vx2, vx3, rho, coords, dx, o, true);
               //initialize(h, x1, x2, x3, f1, f2, f3, conc, coords, dx, oDiffusion, false);


               //distributionsf->setPostCollisionDistribution(f, ix1, ix2, ix3);
               //distributionsf->setPreCollisionDistribution(f, ix1, ix2, ix3);

               real lambda = muLambda.Eval();
               
               calcFeqsFct(h,lambda,c0o1,c0o1,c0o1);
               
               distributions->setPostCollisionDistribution(h, ix1, ix2, ix3);
               distributions->setPreCollisionDistribution(h, ix1, ix2, ix3);


            }
   }

   //variablen der functions loeschen, da die verwiesenen Objecte nach dem verlassen des scopes ungueltig sind!
   //this->muVx1.ClearVar();
   //this->muVx2.ClearVar();
   //this->muVx3.ClearVar();
   //this->muRho.ClearVar();

   this->muLambda.ClearVar();
}
//////////////////////////////////////////////////////////////////////////
void InitThixotropyBlockVisitor::checkFunction(mu::Parser fct)
{
   real x1 = 1.0, x2 = 1.0, x3 = 1.0;
   fct.DefineVar("x1", &x1);
   fct.DefineVar("x2", &x2);
   fct.DefineVar("x3", &x3);

   try
   {
      fct.Eval();
      fct.ClearVar();
   }
   catch (mu::ParserError & e)
   {
      throw UbException(UB_EXARGS, "function: " + e.GetExpr() + (std::string)"error: " + e.GetMsg()
         + (std::string)", only x1,x2,x3 are allowed as variables");
   }
}
//////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setNu(real nu)
//{
//    this->nu = nu;
//}
////////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::setD(real D)
//{
//    this->D = D;
//}
//////////////////////////////////////////////////////////////////////////
//void InitThixotropyBlockVisitor::initialize(double* f, double x1, double x2, double x3, double vx1, double vx2, double vx3, double rho, UbTupleDouble3 coords, double dx, double o, bool NSE)
//{
//   using namespace D3Q27System;
//   //Funktionszeiger
//   typedef void(*CalcFeqsFct)(real* const& /*feq[27]*/, const real& /*(d)rho*/, const real& /*vx1*/, const real& /*vx2*/, const real& /*vx3*/);
//   CalcFeqsFct   calcFeqsFct = NULL;
//
//
//   //if (NSE)
//   //{
//   //    calcFeqsFct = &D3Q27System::calcIncompFeq;
//
//   //}
//   //else
//   //{
//      //calcFeqsFct = &D3Q27System::calcCompFeq;
//      //muVx1 = muf1;
//      //muVx2 = muf2;
//      //muVx3 = muf3;
//   //}
//
//   if (kernel->getCompressible())
//      calcFeqsFct = &D3Q27System::calcCompFeq;
//   else
//      calcFeqsFct = &D3Q27System::calcIncompFeq;
//
//   //x-derivative
//   double deltaX = dx * 0.5;
//   x1 = val<1>(coords) + deltaX;
//   double vx1Plusx1 = muVx1.Eval();
//   double vx2Plusx1 = muVx2.Eval();
//   double vx3Plusx1 = muVx3.Eval();
//
//   x1 = val<1>(coords) - deltaX;
//   double vx1Minusx1 = muVx1.Eval();
//   double vx2Minusx1 = muVx2.Eval();
//   double vx3Minusx1 = muVx3.Eval();
//
//   //y-derivative
//   x1 = val<1>(coords);
//   x2 = val<2>(coords) + deltaX;
//   double vx1Plusx2 = muVx1.Eval();
//   double vx2Plusx2 = muVx2.Eval();
//   double vx3Plusx2 = muVx3.Eval();
//
//   x2 = val<2>(coords) - deltaX;
//   double vx1Minusx2 = muVx1.Eval();
//   double vx2Minusx2 = muVx2.Eval();
//   double vx3Minusx2 = muVx3.Eval();
//
//   //z-derivative
//   x2 = val<2>(coords);
//   x3 = val<3>(coords) + deltaX;
//   double vx1Plusx3 = muVx1.Eval();
//   double vx2Plusx3 = muVx2.Eval();
//   double vx3Plusx3 = muVx3.Eval();
//
//   x3 = val<3>(coords) - deltaX;
//   double vx1Minusx3 = muVx1.Eval();
//   double vx2Minusx3 = muVx2.Eval();
//   double vx3Minusx3 = muVx3.Eval();
//
//   double ax = (vx1Plusx1 - vx1Minusx1) / (2.0 * deltaX) * dx;
//   double bx = (vx2Plusx1 - vx2Minusx1) / (2.0 * deltaX) * dx;
//   double cx = (vx3Plusx1 - vx3Minusx1) / (2.0 * deltaX) * dx;
//
//   double ay = (vx1Plusx2 - vx1Minusx2) / (2.0 * deltaX) * dx;
//   double by = (vx2Plusx2 - vx2Minusx2) / (2.0 * deltaX) * dx;
//   double cy = (vx3Plusx2 - vx3Minusx2) / (2.0 * deltaX) * dx;
//
//   double az = (vx1Plusx3 - vx1Minusx3) / (2.0 * deltaX) * dx;
//   double bz = (vx2Plusx3 - vx2Minusx3) / (2.0 * deltaX) * dx;
//   double cz = (vx3Plusx3 - vx3Minusx3) / (2.0 * deltaX) * dx;
//   double eps_new = 1.0;
//   real op = 1.;
//
//   real feq[27];
//
//   calcFeqsFct(feq, rho, vx1, vx2, vx3);
//
//   double f_E = eps_new * ((5. * ax * o + 5. * by * o + 5. * cz * o - 8. * ax * op + 4. * by * op + 4. * cz * op) / (54. * o * op));
//   double f_N = f_E + eps_new * ((2. * (ax - by)) / (9. * o));
//   double f_T = f_E + eps_new * ((2. * (ax - cz)) / (9. * o));
//   double f_NE = eps_new * (-(5. * cz * o + 3. * (ay + bx) * op - 2. * cz * op + ax * (5. * o + op) + by * (5. * o + op)) / (54. * o * op));
//   double f_SE = f_NE + eps_new * ((ay + bx) / (9. * o));
//   double f_TE = eps_new * (-(5. * cz * o + by * (5. * o - 2. * op) + 3. * (az + cx) * op + cz * op + ax * (5. * o + op)) / (54. * o * op));
//   double f_BE = f_TE + eps_new * ((az + cx) / (9. * o));
//   double f_TN = eps_new * (-(5. * ax * o + 5. * by * o + 5. * cz * o - 2. * ax * op + by * op + 3. * bz * op + 3. * cy * op + cz * op) / (54. * o * op));
//   double f_BN = f_TN + eps_new * ((bz + cy) / (9. * o));
//   double f_ZERO = eps_new * ((5. * (ax + by + cz)) / (9. * op));
//   double f_TNE = eps_new * (-(ay + az + bx + bz + cx + cy) / (72. * o));
//   double f_TSW = -eps_new * ((ay + bx) / (36. * o)) - f_TNE;
//   double f_TSE = -eps_new * ((az + cx) / (36. * o)) - f_TNE;
//   double f_TNW = -eps_new * ((bz + cy) / (36. * o)) - f_TNE;
//
//
//   f[dP00] = f_E + feq[dP00];
//   f[W] = f_E + feq[W];
//   f[N] = f_N + feq[N];
//   f[S] = f_N + feq[S];
//   f[T] = f_T + feq[T];
//   f[B] = f_T + feq[B];
//   f[NE] = f_NE + feq[NE];
//   f[SW] = f_NE + feq[SW];
//   f[SE] = f_SE + feq[SE];
//   f[NW] = f_SE + feq[NW];
//   f[TE] = f_TE + feq[TE];
//   f[BW] = f_TE + feq[BW];
//   f[BE] = f_BE + feq[BE];
//   f[TW] = f_BE + feq[TW];
//   f[TN] = f_TN + feq[TN];
//   f[BS] = f_TN + feq[BS];
//   f[BN] = f_BN + feq[BN];
//   f[TS] = f_BN + feq[TS];
//   f[TNE] = f_TNE + feq[TNE];
//   f[TNW] = f_TNW + feq[TNW];
//   f[TSE] = f_TSE + feq[TSE];
//   f[TSW] = f_TSW + feq[TSW];
//   f[BNE] = f_TSW + feq[BNE];
//   f[BNW] = f_TSE + feq[BNW];
//   f[BSE] = f_TNW + feq[BSE];
//   f[BSW] = f_TNE + feq[BSW];
//   f[REST] = f_ZERO + feq[REST];
//
//
//
//
//}
