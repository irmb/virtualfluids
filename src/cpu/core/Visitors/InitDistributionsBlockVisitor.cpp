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
//! \file InitDistributionsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger
//=======================================================================================

#include "InitDistributionsBlockVisitor.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "LBMKernel.h"

InitDistributionsBlockVisitor::InitDistributionsBlockVisitor() : Block3DVisitor(0, D3Q27System::MAXLEVEL)
{
    using namespace vf::basics::constant;
    
    this->setVx1(c0o1);
    this->setVx2(c0o1);
    this->setVx3(c0o1);
    this->setRho(c0o1);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muVx1 = parser;
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muVx2 = parser;
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muVx3 = parser;
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho(const mu::Parser &parser)
{
    this->checkFunction(parser);
    this->muRho = parser;
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1(const std::string &muParserString)
{
    this->muVx1.SetExpr(muParserString);
    this->checkFunction(muVx1);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2(const std::string &muParserString)
{
    this->muVx2.SetExpr(muParserString);
    this->checkFunction(muVx2);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3(const std::string &muParserString)
{
    this->muVx3.SetExpr(muParserString);
    this->checkFunction(muVx3);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho(const std::string &muParserString)
{
    this->muRho.SetExpr(muParserString);
    this->checkFunction(muRho);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1(real vx1)
{
    this->muVx1.SetExpr(UbSystem::toString(vx1, realLim::digits10));
    this->checkFunction(muVx1);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2(real vx2)
{
    this->muVx2.SetExpr(UbSystem::toString(vx2, realLim::digits10));
    this->checkFunction(muVx2);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3(real vx3)
{
    this->muVx3.SetExpr(UbSystem::toString(vx3, realLim::digits10));
    this->checkFunction(muVx3);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho(real rho)
{
    this->muRho.SetExpr(UbSystem::toString(rho, realLim::digits10));
    this->checkFunction(muRho);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   if(!block) UB_THROW( UbException(UB_EXARGS,"block is not exist") );

   real dx = grid->getDeltaX(block);

   //define vars for functions
   mu::value_type x1,x2,x3;
   this->muVx1.DefineVar("x1",&x1); this->muVx1.DefineVar("x2",&x2); this->muVx1.DefineVar("x3",&x3);
   this->muVx2.DefineVar("x1",&x1); this->muVx2.DefineVar("x2",&x2); this->muVx2.DefineVar("x3",&x3);
   this->muVx3.DefineVar("x1",&x1); this->muVx3.DefineVar("x2",&x2); this->muVx3.DefineVar("x3",&x3);
   this->muRho.DefineVar("x1",&x1); this->muRho.DefineVar("x2",&x2); this->muRho.DefineVar("x3",&x3);

    using CalcFeqsFct = void (*)(real *const & /*feq[27]*/, const real & /*(d)rho*/, const real & /*vx1*/,
                                 const real & /*vx2*/, const real & /*vx3*/);
    CalcFeqsFct calcFeqsFct = NULL;
   
   real vx1, vx2, vx3, rho;

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
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();  

      real o  = kernel->getCollisionFactor();

      real f[D3Q27System::ENDF+1];

      for(std::size_t ix3=0; ix3<bcArray->getNX3(); ix3++)
         for(std::size_t ix2=0; ix2<bcArray->getNX2(); ix2++)
            for(std::size_t ix1=0; ix1<bcArray->getNX1(); ix1++)
            {
               Vector3D coords = grid->getNodeCoordinates(block, (int)ix1, (int)ix2, (int)ix3);
               x1 = coords[0];
               x2 = coords[1];
               x3 = coords[2];

               vx1 = muVx1.Eval();
               vx2 = muVx2.Eval();
               vx3 = muVx3.Eval();
               rho = muRho.Eval();

               //x-derivative
               real deltaX=dx*c1o2;
               x1 = coords[0]+deltaX;
               real vx1Plusx1 = muVx1.Eval();
               real vx2Plusx1 = muVx2.Eval();
               real vx3Plusx1 = muVx3.Eval();

               x1 = coords[0]-deltaX;
               real vx1Minusx1 = muVx1.Eval();
               real vx2Minusx1 = muVx2.Eval();
               real vx3Minusx1 = muVx3.Eval();

               //y-derivative
               x1 = coords[0];
               x2 = coords[1]+deltaX;
               real vx1Plusx2 = muVx1.Eval();
               real vx2Plusx2 = muVx2.Eval();
               real vx3Plusx2 = muVx3.Eval();

               x2 = coords[1]-deltaX;
               real vx1Minusx2 = muVx1.Eval();
               real vx2Minusx2 = muVx2.Eval();
               real vx3Minusx2 = muVx3.Eval();

               //z-derivative
               x2 = coords[1];
               x3 = coords[2]+deltaX;
               real vx1Plusx3 = muVx1.Eval();
               real vx2Plusx3 = muVx2.Eval();
               real vx3Plusx3 = muVx3.Eval();

               x3 = coords[2]-deltaX;
               real vx1Minusx3 = muVx1.Eval();
               real vx2Minusx3 = muVx2.Eval();
               real vx3Minusx3 = muVx3.Eval();

               real ax=(vx1Plusx1-vx1Minusx1)/(c2o1*deltaX)*dx;
               real bx=(vx2Plusx1-vx2Minusx1)/(c2o1*deltaX)*dx;
               real cx=(vx3Plusx1-vx3Minusx1)/(c2o1*deltaX)*dx;

               real ay=(vx1Plusx2-vx1Minusx2)/(c2o1*deltaX)*dx;
               real by=(vx2Plusx2-vx2Minusx2)/(c2o1*deltaX)*dx;
               real cy=(vx3Plusx2-vx3Minusx2)/(c2o1*deltaX)*dx;

               real az=(vx1Plusx3-vx1Minusx3)/(c2o1*deltaX)*dx;
               real bz=(vx2Plusx3-vx2Minusx3)/(c2o1*deltaX)*dx;
               real cz=(vx3Plusx3-vx3Minusx3)/(c2o1*deltaX)*dx;
               real eps_new=c1o1;
               real op = c1o1;

               real feq[27];

               calcFeqsFct(feq,rho,vx1,vx2,vx3);

               real f_E    = eps_new *((c5o1*ax*o + c5o1*by*o + c5o1*cz*o - c8o1*ax*op + c4o1*by*op + c4o1*cz*op)/(c54o1*o*op));
               real f_N    = f_E + eps_new *((c2o1*(ax - by))/(c9o1*o));
               real f_T    = f_E + eps_new *((c2o1*(ax - cz))/(c9o1*o));
               real f_NE   = eps_new *(-(c5o1*cz*o + c3o1*(ay + bx)*op - c2o1*cz*op + ax*(c5o1*o + op) + by*(c5o1*o + op))/(c54o1*o*op));
               real f_SE   = f_NE + eps_new *((  ay + bx )/(c9o1*o));
               real f_TE   = eps_new *(-(c5o1*cz*o + by*(c5o1*o - c2o1*op) + c3o1*(az + cx)*op + cz*op + ax*(c5o1*o + op))/(c54o1*o*op));
               real f_BE   = f_TE + eps_new *((  az + cx )/(c9o1*o));
               real f_TN   = eps_new *(-(c5o1*ax*o + c5o1*by*o + c5o1*cz*o - c2o1*ax*op + by*op + c3o1*bz*op + c3o1*cy*op + cz*op)/(c54o1*o*op));
               real f_BN   = f_TN + eps_new *((  bz + cy )/(c9o1*o));
               real f_ZERO = eps_new *((c5o1*(ax + by + cz))/(c9o1*op));
               real f_TNE  = eps_new *(-(ay + az + bx + bz + cx + cy)/(c72o1*o));
               real f_TSW  = - eps_new *((ay + bx)/(c36o1*o)) - f_TNE;
               real f_TSE  = - eps_new *((az + cx)/(c36o1*o)) - f_TNE;
               real f_TNW  = - eps_new *((bz + cy)/(c36o1*o)) - f_TNE;


               f[dP00]    = f_E    + feq[dP00];
               f[dM00]    = f_E    + feq[dM00];
               f[d0P0]    = f_N    + feq[d0P0];
               f[d0M0]    = f_N    + feq[d0M0];
               f[d00P]    = f_T    + feq[d00P];
               f[d00M]    = f_T    + feq[d00M];
               f[dPP0]   = f_NE   + feq[dPP0];
               f[dMM0]   = f_NE   + feq[dMM0];
               f[dPM0]   = f_SE   + feq[dPM0];
               f[dMP0]   = f_SE   + feq[dMP0];
               f[dP0P]   = f_TE   + feq[dP0P];
               f[dM0M]   = f_TE   + feq[dM0M];
               f[dP0M]   = f_BE   + feq[dP0M];
               f[dM0P]   = f_BE   + feq[dM0P];
               f[d0PP]   = f_TN   + feq[d0PP];
               f[d0MM]   = f_TN   + feq[d0MM];
               f[d0PM]   = f_BN   + feq[d0PM];
               f[d0MP]   = f_BN   + feq[d0MP];
               f[dPPP]  = f_TNE  + feq[dPPP];
               f[dMPP]  = f_TNW  + feq[dMPP];
               f[dPMP]  = f_TSE  + feq[dPMP];
               f[dMMP]  = f_TSW  + feq[dMMP];
               f[dPPM]  = f_TSW  + feq[dPPM];
               f[dMPM]  = f_TSE  + feq[dMPM];
               f[dPMM]  = f_TNW  + feq[dPMM];
               f[dMMM]  = f_TNE  + feq[dMMM];
               f[d000] = f_ZERO + feq[d000];

               //calcFeqsFct(f,rho,vx1,vx2,vx3);
               distributions->setPostCollisionDistribution(f, ix1, ix2, ix3);
               distributions->setPreCollisionDistribution(f, ix1, ix2, ix3);

               //distributions->swap();
               //distributions->setPostCollisionDistribution(f, ix1, ix2, ix3);
               //distributions->setPreCollisionDistribution(f, ix1, ix2, ix3);
               //distributions->swap();

            }
   }

   //variablen der functions loeschen, da die verwiesenen Objecte nach dem verlassen des scopes ungueltig sind!
   this->muVx1.ClearVar();
   this->muVx2.ClearVar();
   this->muVx3.ClearVar();
   this->muRho.ClearVar();

}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::checkFunction(mu::Parser fct)
{
    using namespace vf::basics::constant;
    
    double x1 = c1o1, x2 = c1o1, x3 = c1o1;
    fct.DefineVar("x1", &x1);
    fct.DefineVar("x2", &x2);
    fct.DefineVar("x3", &x3);

    try {
        fct.Eval();
        fct.ClearVar();
    } catch (mu::ParserError &e) {
        throw UbException(UB_EXARGS, "function: " + e.GetExpr() + (std::string) "error: " + e.GetMsg() +
                                         (std::string) ", only x1,x2,x3 are allowed as variables");
    }
}
//////////////////////////////////////////////////////////////////////////
