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
//! \file MultiphaseInitDistributionsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseInitDistributionsBlockVisitor.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "LBMKernel.h"

MultiphaseInitDistributionsBlockVisitor::MultiphaseInitDistributionsBlockVisitor() 
	: Block3DVisitor(0, D3Q27System::MAXLEVEL)
{
	this->setVx1(0.0);
	this->setVx2(0.0);
	this->setVx3(0.0);
	this->setRho(0.0);
}
//////////////////////////////////////////////////////////////////////////
MultiphaseInitDistributionsBlockVisitor::MultiphaseInitDistributionsBlockVisitor( LBMReal densityRatio, LBMReal vx1, LBMReal vx2, LBMReal vx3, LBMReal rho)
	: Block3DVisitor(0, D3Q27System::MAXLEVEL), densityRatio(densityRatio) 
{
	this->setVx1(vx1);
	this->setVx2(vx2);
	this->setVx3(vx3);
	this->setRho(rho);}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx1( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muVx1 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx2( const mu::Parser& parser)
{ 
	this->checkFunction(parser); 
	this->muVx2 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx3( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muVx3 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setRho( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muRho = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setPhi( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muPhi = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx1( const std::string& muParserString)  
{ 
	this->muVx1.SetExpr(muParserString); 
	this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx2( const std::string& muParserString) 
{ 
	this->muVx2.SetExpr(muParserString); 
	this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx3( const std::string& muParserString)  
{ 
	this->muVx3.SetExpr(muParserString); 
	this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setRho( const std::string& muParserString)  
{ 
	this->muRho.SetExpr(muParserString); 
	this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setPhi( const std::string& muParserString)  
{ 
	this->muPhi.SetExpr(muParserString); 
	this->checkFunction(muPhi); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx1( LBMReal vx1 ) 
{ 
	this->muVx1.SetExpr( UbSystem::toString(vx1,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx2( LBMReal vx2 ) 
{ 
	this->muVx2.SetExpr( UbSystem::toString(vx2,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx3( LBMReal vx3 ) 
{ 
	this->muVx3.SetExpr( UbSystem::toString(vx3,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setRho( LBMReal rho ) 
{ 
	this->muRho.SetExpr( UbSystem::toString(rho,D3Q27RealLim::digits10) );  
	this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setPhi( LBMReal phi ) 
{ 
	this->muPhi.SetExpr( UbSystem::toString(phi,D3Q27RealLim::digits10) );  
	this->checkFunction(muPhi); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block) 
{
	using namespace D3Q27System;

	if(!block) UB_THROW( UbException(UB_EXARGS,"block is not exist") );

	//define vars for functions
	mu::value_type x1,x2,x3;
	this->muVx1.DefineVar("x1",&x1); this->muVx1.DefineVar("x2",&x2); this->muVx1.DefineVar("x3",&x3);
	this->muVx2.DefineVar("x1",&x1); this->muVx2.DefineVar("x2",&x2); this->muVx2.DefineVar("x3",&x3);
	this->muVx3.DefineVar("x1",&x1); this->muVx3.DefineVar("x2",&x2); this->muVx3.DefineVar("x3",&x3);
	this->muRho.DefineVar("x1",&x1); this->muRho.DefineVar("x2",&x2); this->muRho.DefineVar("x3",&x3);
	this->muPhi.DefineVar("x1",&x1); this->muPhi.DefineVar("x2",&x2); this->muPhi.DefineVar("x3",&x3);

	LBMReal vx1, vx2, vx3, rho, /*p1,*/ phi;

	int gridRank = grid->getRank();
	int blockRank = block->getRank();

	if (blockRank == gridRank && block->isActive())
	{
        SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
		if (!kernel)
			throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

		SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
        SPtr<EsoTwist3D> distributionsF = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getFdistributions()); 
		SPtr<EsoTwist3D> distributionsH = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getHdistributions());
        SPtr<EsoTwist3D> distributionsH2 = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getH2distributions());

		LBMReal phiL = kernel->getPhiL();
		LBMReal phiH = kernel->getPhiH();

		LBMReal f[D3Q27System::ENDF+1];

		for(int ix3=0; ix3<(int)bcArray->getNX3(); ix3++)
            for (int ix2 = 0; ix2 < (int)bcArray->getNX2(); ix2++)
                for (int ix1 = 0; ix1 < (int)bcArray->getNX1(); ix1++)
				{
					Vector3D coords = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    x1              = coords[0];
                    x2              = coords[1];
                    x3              = coords[2];

					
					//p1  = 0.0;
					//p1 = muRho.Eval();
					vx1 = muVx1.Eval();
					vx2 = muVx2.Eval();
					vx3 = muVx3.Eval();
					phi = muPhi.Eval();
					
					//rho = phi*1.0 + (1.0-phi)/densityRatio;
					LBMReal rhoH = 1.0;
					LBMReal rhoL = 1.0/densityRatio;
					rho = rhoH + (rhoH - rhoL)*(phi - phiH)/(phiH - phiL);

			
					LBMReal feq[27];
					LBMReal geq[27];

					//calcFeqsFct(feq,rho,vx1,vx2,vx3);
					LBMReal vx1Sq = vx1*vx1;
					LBMReal vx2Sq = vx2*vx2;
					LBMReal vx3Sq = vx3*vx3;
					for (int dir = STARTF; dir < (ENDF+1); dir++)
					{
						LBMReal velProd = DX1[dir]*vx1 + DX2[dir]*vx2 + DX3[dir]*vx3;
						LBMReal velSq1 = velProd*velProd;
						LBMReal gamma = WEIGTH[dir]*(3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));

						feq[dir] = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));
						//geq[dir] = p1*WEIGTH[dir] + gamma;
						//geq[dir] = p1*WEIGTH[dir]/(rho*UbMath::c1o3) + gamma*rho;
						//geq[dir] = (p1*WEIGTH[dir]/(rho*UbMath::c1o3) + gamma*rho)*UbMath::c1o3;
						geq[dir] = (gamma*rho)*UbMath::c1o3;
					}


					f[E]    =  geq[E]    ;
					f[W]    =  geq[W]    ;
					f[N]    =  geq[N]    ;
					f[S]    =  geq[S]    ;
					f[T]    =  geq[T]    ;
					f[B]    =  geq[B]    ;
					f[NE]   =  geq[NE]   ;
					f[SW]   =  geq[SW]   ;
					f[SE]   =  geq[SE]   ;
					f[NW]   =  geq[NW]   ;
					f[TE]   =  geq[TE]   ;
					f[BW]   =  geq[BW]   ;
					f[BE]   =  geq[BE]   ;
					f[TW]   =  geq[TW]   ;
					f[TN]   =  geq[TN]   ;
					f[BS]   =  geq[BS]   ;
					f[BN]   =  geq[BN]   ;
					f[TS]   =  geq[TS]   ;
					f[TNE]  =  geq[TNE]  ;
					f[TNW]  =  geq[TNW]  ;
					f[TSE]  =  geq[TSE]  ;
					f[TSW]  =  geq[TSW]  ;
					f[BNE]  =  geq[BNE]  ;
					f[BNW]  =  geq[BNW]  ;
					f[BSE]  =  geq[BSE]  ;
					f[BSW]  =  geq[BSW]  ;
					f[DIR_000] =  geq[DIR_000] ;

					distributionsF->setDistribution(f, ix1, ix2, ix3);
					distributionsF->setDistributionInv(f, ix1, ix2, ix3);

					f[E]    =  phi * feq[E]    / rho;
					f[W]    =  phi * feq[W]    / rho;
					f[N]    =  phi * feq[N]    / rho;
					f[S]    =  phi * feq[S]    / rho;
					f[T]    =  phi * feq[T]    / rho;
					f[B]    =  phi * feq[B]    / rho;
					f[NE]   =  phi * feq[NE]   / rho;
					f[SW]   =  phi * feq[SW]   / rho;
					f[SE]   =  phi * feq[SE]   / rho;
					f[NW]   =  phi * feq[NW]   / rho;
					f[TE]   =  phi * feq[TE]   / rho;
					f[BW]   =  phi * feq[BW]   / rho;
					f[BE]   =  phi * feq[BE]   / rho;
					f[TW]   =  phi * feq[TW]   / rho;
					f[TN]   =  phi * feq[TN]   / rho;
					f[BS]   =  phi * feq[BS]   / rho;
					f[BN]   =  phi * feq[BN]   / rho;
					f[TS]   =  phi * feq[TS]   / rho;
					f[TNE]  =  phi * feq[TNE]  / rho;
					f[TNW]  =  phi * feq[TNW]  / rho;
					f[TSE]  =  phi * feq[TSE]  / rho;
					f[TSW]  =  phi * feq[TSW]  / rho;
					f[BNE]  =  phi * feq[BNE]  / rho;
					f[BNW]  =  phi * feq[BNW]  / rho;
					f[BSE]  =  phi * feq[BSE]  / rho;
					f[BSW]  =  phi * feq[BSW]  / rho;
					f[DIR_000] =  phi * feq[DIR_000] / rho;

					distributionsH->setDistribution(f, ix1, ix2, ix3);
					distributionsH->setDistributionInv(f, ix1, ix2, ix3);

					if (distributionsH2) {

						f[E]    = (1.-phi) * feq[E] / rho;
						f[W]    = (1.-phi) * feq[W] / rho;
						f[N]    = (1.-phi) * feq[N] / rho;
						f[S]    = (1.-phi) * feq[S] / rho;
						f[T]    = (1.-phi) * feq[T] / rho;
						f[B]    = (1.-phi) * feq[B] / rho;
						f[NE]   = (1.-phi) * feq[NE] / rho;
						f[SW]   = (1.-phi) * feq[SW] / rho;
						f[SE]   = (1.-phi) * feq[SE] / rho;
						f[NW]   = (1.-phi) * feq[NW] / rho;
						f[TE]   = (1.-phi) * feq[TE] / rho;
						f[BW]   = (1.-phi) * feq[BW] / rho;
						f[BE]   = (1.-phi) * feq[BE] / rho;
						f[TW]   = (1.-phi) * feq[TW] / rho;
						f[TN]   = (1.-phi) * feq[TN] / rho;
						f[BS]   = (1.-phi) * feq[BS] / rho;
						f[BN]   = (1.-phi) * feq[BN] / rho;
						f[TS]   = (1.-phi) * feq[TS] / rho;
						f[TNE]  = (1.-phi) * feq[TNE] / rho;
						f[TNW]  = (1.-phi) * feq[TNW] / rho;
						f[TSE]  = (1.-phi) * feq[TSE] / rho;
						f[TSW]  = (1.-phi) * feq[TSW] / rho;
						f[BNE]  = (1.-phi) * feq[BNE] / rho;
						f[BNW]  = (1.-phi) * feq[BNW] / rho;
						f[BSE]  = (1.-phi) * feq[BSE] / rho;
						f[BSW]  = (1.-phi) * feq[BSW] / rho;
						f[DIR_000] = (1.-phi) * feq[DIR_000] / rho;

                        distributionsH2->setDistribution(f, ix1, ix2, ix3);
                        distributionsH2->setDistributionInv(f, ix1, ix2, ix3);                    
					}
				}
	}

	//variablen der functions loeschen, da die verwiesenen Objecte nach dem verlassen des scopes ungueltig sind!
	this->muVx1.ClearVar();
	this->muVx2.ClearVar();
	this->muVx3.ClearVar();
	this->muRho.ClearVar();

}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::checkFunction(mu::Parser fct)
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
void MultiphaseInitDistributionsBlockVisitor::setNu( LBMReal nu )
{
	this->nu = nu;
}

