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
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "Grid3DSystem.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "BCArray3D.h"

InitDistributionsBlockVisitor::InitDistributionsBlockVisitor() 
   : Block3DVisitor(0, Grid3DSystem::MAXLEVEL)
{
   this->setVx1(0.0);
   this->setVx2(0.0);
   this->setVx3(0.0);
   this->setRho(0.0);
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muVx1 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2( const mu::Parser& parser)
{ 
   this->checkFunction(parser); 
   this->muVx2 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muVx3 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muRho = parser;  
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1( const std::string& muParserString)  
{ 
   this->muVx1.SetExpr(muParserString); 
   this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2( const std::string& muParserString) 
{ 
   this->muVx2.SetExpr(muParserString); 
   this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3( const std::string& muParserString)  
{ 
   this->muVx3.SetExpr(muParserString); 
   this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho( const std::string& muParserString)  
{ 
   this->muRho.SetExpr(muParserString); 
   this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx1( LBMReal vx1 ) 
{ 
   this->muVx1.SetExpr( UbSystem::toString(vx1,D3Q27RealLim::digits10) );  
   this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx2( LBMReal vx2 ) 
{ 
   this->muVx2.SetExpr( UbSystem::toString(vx2,D3Q27RealLim::digits10) );  
   this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setVx3( LBMReal vx3 ) 
{ 
   this->muVx3.SetExpr( UbSystem::toString(vx3,D3Q27RealLim::digits10) );  
   this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::setRho( LBMReal rho ) 
{ 
   this->muRho.SetExpr( UbSystem::toString(rho,D3Q27RealLim::digits10) );  
   this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void InitDistributionsBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block) 
{
   using namespace D3Q27System;

   if(!block) UB_THROW( UbException(UB_EXARGS,"block is not exist") );

   //UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
   //UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
   double dx = grid->getDeltaX(block);

   //define vars for functions
   mu::value_type x1,x2,x3;
   this->muVx1.DefineVar("x1",&x1); this->muVx1.DefineVar("x2",&x2); this->muVx1.DefineVar("x3",&x3);
   this->muVx2.DefineVar("x1",&x1); this->muVx2.DefineVar("x2",&x2); this->muVx2.DefineVar("x3",&x3);
   this->muVx3.DefineVar("x1",&x1); this->muVx3.DefineVar("x2",&x2); this->muVx3.DefineVar("x3",&x3);
   this->muRho.DefineVar("x1",&x1); this->muRho.DefineVar("x2",&x2); this->muRho.DefineVar("x3",&x3);

   //Funktionszeiger
   typedef void (*CalcFeqsFct)(LBMReal* const& /*feq[27]*/,const LBMReal& /*(d)rho*/,const LBMReal& /*vx1*/,const LBMReal& /*vx2*/,const LBMReal& /*vx3*/);
   CalcFeqsFct   calcFeqsFct   = NULL;
   
   LBMReal vx1,vx2,vx3,rho;

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

      //UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);

      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();  

      LBMReal o  = kernel->getCollisionFactor();

      LBMReal f[D3Q27System::ENDF+1];

      size_t nx1 = distributions->getNX1();
      size_t nx2 = distributions->getNX2();
      size_t nx3 = distributions->getNX3();

      for(int ix3=0; ix3<bcArray->getNX3(); ix3++)
         for(int ix2=0; ix2<bcArray->getNX2(); ix2++)
            for(int ix1=0; ix1<bcArray->getNX1(); ix1++)
            {
               Vector3D coords = grid->getNodeCoordinates(block, ix1, ix2, ix3);
               x1 = coords[0];
               x2 = coords[1];
               x3 = coords[2];

               vx1 = muVx1.Eval();
               vx2 = muVx2.Eval();
               vx3 = muVx3.Eval();
               rho = muRho.Eval();

               //x-derivative
               double deltaX=dx*0.5;
               x1 = coords[0]+deltaX;
               double vx1Plusx1 = muVx1.Eval();
               double vx2Plusx1 = muVx2.Eval();
               double vx3Plusx1 = muVx3.Eval();

               x1 = coords[0]-deltaX;
               double vx1Minusx1 = muVx1.Eval();
               double vx2Minusx1 = muVx2.Eval();
               double vx3Minusx1 = muVx3.Eval();

               //y-derivative
               x1 = coords[0];
               x2 = coords[1]+deltaX;
               double vx1Plusx2 = muVx1.Eval();
               double vx2Plusx2 = muVx2.Eval();
               double vx3Plusx2 = muVx3.Eval();

               x2 = coords[1]-deltaX;
               double vx1Minusx2 = muVx1.Eval();
               double vx2Minusx2 = muVx2.Eval();
               double vx3Minusx2 = muVx3.Eval();

               //z-derivative
               x2 = coords[1];
               x3 = coords[2]+deltaX;
               double vx1Plusx3 = muVx1.Eval();
               double vx2Plusx3 = muVx2.Eval();
               double vx3Plusx3 = muVx3.Eval();

               x3 = coords[2]-deltaX;
               double vx1Minusx3 = muVx1.Eval();
               double vx2Minusx3 = muVx2.Eval();
               double vx3Minusx3 = muVx3.Eval();

               double ax=(vx1Plusx1-vx1Minusx1)/(2.0*deltaX)*dx;
               double bx=(vx2Plusx1-vx2Minusx1)/(2.0*deltaX)*dx;
               double cx=(vx3Plusx1-vx3Minusx1)/(2.0*deltaX)*dx;

               double ay=(vx1Plusx2-vx1Minusx2)/(2.0*deltaX)*dx;
               double by=(vx2Plusx2-vx2Minusx2)/(2.0*deltaX)*dx;
               double cy=(vx3Plusx2-vx3Minusx2)/(2.0*deltaX)*dx;

               double az=(vx1Plusx3-vx1Minusx3)/(2.0*deltaX)*dx;
               double bz=(vx2Plusx3-vx2Minusx3)/(2.0*deltaX)*dx;
               double cz=(vx3Plusx3-vx3Minusx3)/(2.0*deltaX)*dx;
               double eps_new=1.0;
               LBMReal op = 1.;

               LBMReal feq[27];

               calcFeqsFct(feq,rho,vx1,vx2,vx3);

               double f_E    = eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op)/(54.*o*op));
               double f_N    = f_E + eps_new *((2.*(ax - by))/(9.*o));
               double f_T    = f_E + eps_new *((2.*(ax - cz))/(9.*o));
               double f_NE   = eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op))/(54.*o*op));
               double f_SE   = f_NE + eps_new *((  ay + bx )/(9.*o));
               double f_TE   = eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op))/(54.*o*op));
               double f_BE   = f_TE + eps_new *((  az + cx )/(9.*o));
               double f_TN   = eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op)/(54.*o*op));
               double f_BN   = f_TN + eps_new *((  bz + cy )/(9.*o));
               double f_ZERO = eps_new *((5.*(ax + by + cz))/(9.*op));
               double f_TNE  = eps_new *(-(ay + az + bx + bz + cx + cy)/(72.*o));
               double f_TSW  = - eps_new *((ay + bx)/(36.*o)) - f_TNE;
               double f_TSE  = - eps_new *((az + cx)/(36.*o)) - f_TNE;
               double f_TNW  = - eps_new *((bz + cy)/(36.*o)) - f_TNE;


               f[E]    = f_E    + feq[E];
               f[W]    = f_E    + feq[W];
               f[N]    = f_N    + feq[N];
               f[S]    = f_N    + feq[S];
               f[T]    = f_T    + feq[T];
               f[B]    = f_T    + feq[B];
               f[NE]   = f_NE   + feq[NE];
               f[SW]   = f_NE   + feq[SW];
               f[SE]   = f_SE   + feq[SE];
               f[NW]   = f_SE   + feq[NW];
               f[TE]   = f_TE   + feq[TE];
               f[BW]   = f_TE   + feq[BW];
               f[BE]   = f_BE   + feq[BE];
               f[TW]   = f_BE   + feq[TW];
               f[TN]   = f_TN   + feq[TN];
               f[BS]   = f_TN   + feq[BS];
               f[BN]   = f_BN   + feq[BN];
               f[TS]   = f_BN   + feq[TS];
               f[TNE]  = f_TNE  + feq[TNE];
               f[TNW]  = f_TNW  + feq[TNW];
               f[TSE]  = f_TSE  + feq[TSE];
               f[TSW]  = f_TSW  + feq[TSW];
               f[BNE]  = f_TSW  + feq[BNE];
               f[BNW]  = f_TSE  + feq[BNW];
               f[BSE]  = f_TNW  + feq[BSE];
               f[BSW]  = f_TNE  + feq[BSW];
               f[ZERO] = f_ZERO + feq[ZERO];

               //calcFeqsFct(f,rho,vx1,vx2,vx3);
               //distributions->setDistribution(f, ix1, ix2, ix3);
               distributions->setDistributionInv(f, ix1, ix2, ix3);

               //distributions->swap();
               //distributions->setDistribution(f, ix1, ix2, ix3);
               //distributions->setDistributionInv(f, ix1, ix2, ix3);
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
   double x1=1.0,x2=1.0,x3=1.0;
   fct.DefineVar("x1",&x1); 
   fct.DefineVar("x2",&x2); 
   fct.DefineVar("x3",&x3);

   try
   {
      fct.Eval();
      fct.ClearVar();
   }
   catch(mu::ParserError& e)
   {
      throw UbException(UB_EXARGS,"function: "+e.GetExpr() + (std::string)"error: "+e.GetMsg()
         +(std::string)", only x1,x2,x3 are allowed as variables" );
   }
}
//////////////////////////////////////////////////////////////////////////



