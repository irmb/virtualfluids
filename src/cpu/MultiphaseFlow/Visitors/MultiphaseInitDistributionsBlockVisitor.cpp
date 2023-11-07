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
#include "BCSet.h"
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
MultiphaseInitDistributionsBlockVisitor::MultiphaseInitDistributionsBlockVisitor( real densityRatio, real vx1, real vx2, real vx3, real rho)
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
void MultiphaseInitDistributionsBlockVisitor::setVx1( real vx1 ) 
{ 
	this->muVx1.SetExpr( UbSystem::toString(vx1,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx2( real vx2 ) 
{ 
	this->muVx2.SetExpr( UbSystem::toString(vx2,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setVx3( real vx3 ) 
{ 
	this->muVx3.SetExpr( UbSystem::toString(vx3,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setRho( real rho ) 
{ 
	this->muRho.SetExpr( UbSystem::toString(rho,D3Q27RealLim::digits10) );  
	this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::setPhi( real phi ) 
{ 
	this->muPhi.SetExpr( UbSystem::toString(phi,D3Q27RealLim::digits10) );  
	this->checkFunction(muPhi); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseInitDistributionsBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block) 
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	if(!block) UB_THROW( UbException(UB_EXARGS,"block is not exist") );

	//define vars for functions
	mu::value_type x1,x2,x3;
	this->muVx1.DefineVar("x1",&x1); this->muVx1.DefineVar("x2",&x2); this->muVx1.DefineVar("x3",&x3);
	this->muVx2.DefineVar("x1",&x1); this->muVx2.DefineVar("x2",&x2); this->muVx2.DefineVar("x3",&x3);
	this->muVx3.DefineVar("x1",&x1); this->muVx3.DefineVar("x2",&x2); this->muVx3.DefineVar("x3",&x3);
	this->muRho.DefineVar("x1",&x1); this->muRho.DefineVar("x2",&x2); this->muRho.DefineVar("x3",&x3);
	this->muPhi.DefineVar("x1",&x1); this->muPhi.DefineVar("x2",&x2); this->muPhi.DefineVar("x3",&x3);

	real vx1, vx2, vx3, rho, /*p1,*/ phi;

	int gridRank = grid->getRank();
	int blockRank = block->getRank();

	if (blockRank == gridRank && block->isActive())
	{
        SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
		if (!kernel)
			throw UbException(UB_EXARGS, "The LBM kernel isn't exist in block: "+block->toString());

		SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();
        SPtr<EsoTwist3D> distributionsF = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getFdistributions()); 
		SPtr<EsoTwist3D> distributionsH = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getHdistributions());
        SPtr<EsoTwist3D> distributionsH2 = dynamicPointerCast<EsoTwist3D>(kernel->getDataSet()->getH2distributions());

		real phiL = kernel->getPhiL();
		real phiH = kernel->getPhiH();

		real f[D3Q27System::ENDF+1];

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
					real rhoH = 1.0;
					real rhoL = 1.0/densityRatio;
					rho = rhoH + (rhoH - rhoL)*(phi - phiH)/(phiH - phiL);

			
					real feq[27];
					real geq[27];

					//calcFeqsFct(feq,rho,vx1,vx2,vx3);
					real vx1Sq = vx1*vx1;
					real vx2Sq = vx2*vx2;
					real vx3Sq = vx3*vx3;
					for (int dir = STARTF; dir < (ENDF+1); dir++)
					{
						real velProd = DX1[dir]*vx1 + DX2[dir]*vx2 + DX3[dir]*vx3;
						real velSq1 = velProd*velProd;
						real gamma = WEIGTH[dir]*(3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));

						feq[dir] = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));
						//geq[dir] = p1*WEIGTH[dir] + gamma;
						//geq[dir] = p1*WEIGTH[dir]/(rho*UbMath::c1o3) + gamma*rho;
						//geq[dir] = (p1*WEIGTH[dir]/(rho*UbMath::c1o3) + gamma*rho)*UbMath::c1o3;
						geq[dir] = (gamma*rho)* vf::basics::constant::c1o3;
					}


					f[dP00]    =  geq[dP00]    ;
					f[dM00]    =  geq[dM00]    ;
					f[DIR_0P0]    =  geq[DIR_0P0]    ;
					f[DIR_0M0]    =  geq[DIR_0M0]    ;
					f[DIR_00P]    =  geq[DIR_00P]    ;
					f[DIR_00M]    =  geq[DIR_00M]    ;
					f[DIR_PP0]   =  geq[DIR_PP0]   ;
					f[DIR_MM0]   =  geq[DIR_MM0]   ;
					f[DIR_PM0]   =  geq[DIR_PM0]   ;
					f[DIR_MP0]   =  geq[DIR_MP0]   ;
					f[DIR_P0P]   =  geq[DIR_P0P]   ;
					f[DIR_M0M]   =  geq[DIR_M0M]   ;
					f[DIR_P0M]   =  geq[DIR_P0M]   ;
					f[DIR_M0P]   =  geq[DIR_M0P]   ;
					f[DIR_0PP]   =  geq[DIR_0PP]   ;
					f[DIR_0MM]   =  geq[DIR_0MM]   ;
					f[DIR_0PM]   =  geq[DIR_0PM]   ;
					f[DIR_0MP]   =  geq[DIR_0MP]   ;
					f[DIR_PPP]  =  geq[DIR_PPP]  ;
					f[DIR_MPP]  =  geq[DIR_MPP]  ;
					f[DIR_PMP]  =  geq[DIR_PMP]  ;
					f[DIR_MMP]  =  geq[DIR_MMP]  ;
					f[DIR_PPM]  =  geq[DIR_PPM]  ;
					f[DIR_MPM]  =  geq[DIR_MPM]  ;
					f[DIR_PMM]  =  geq[DIR_PMM]  ;
					f[DIR_MMM]  =  geq[DIR_MMM]  ;
					f[d000] =  geq[d000] ;

					distributionsF->setPostCollisionDistribution(f, ix1, ix2, ix3);
					distributionsF->setPreCollisionDistribution(f, ix1, ix2, ix3);

					f[dP00]    =  phi * feq[dP00]    / rho;
					f[dM00]    =  phi * feq[dM00]    / rho;
					f[DIR_0P0]    =  phi * feq[DIR_0P0]    / rho;
					f[DIR_0M0]    =  phi * feq[DIR_0M0]    / rho;
					f[DIR_00P]    =  phi * feq[DIR_00P]    / rho;
					f[DIR_00M]    =  phi * feq[DIR_00M]    / rho;
					f[DIR_PP0]   =  phi * feq[DIR_PP0]   / rho;
					f[DIR_MM0]   =  phi * feq[DIR_MM0]   / rho;
					f[DIR_PM0]   =  phi * feq[DIR_PM0]   / rho;
					f[DIR_MP0]   =  phi * feq[DIR_MP0]   / rho;
					f[DIR_P0P]   =  phi * feq[DIR_P0P]   / rho;
					f[DIR_M0M]   =  phi * feq[DIR_M0M]   / rho;
					f[DIR_P0M]   =  phi * feq[DIR_P0M]   / rho;
					f[DIR_M0P]   =  phi * feq[DIR_M0P]   / rho;
					f[DIR_0PP]   =  phi * feq[DIR_0PP]   / rho;
					f[DIR_0MM]   =  phi * feq[DIR_0MM]   / rho;
					f[DIR_0PM]   =  phi * feq[DIR_0PM]   / rho;
					f[DIR_0MP]   =  phi * feq[DIR_0MP]   / rho;
					f[DIR_PPP]  =  phi * feq[DIR_PPP]  / rho;
					f[DIR_MPP]  =  phi * feq[DIR_MPP]  / rho;
					f[DIR_PMP]  =  phi * feq[DIR_PMP]  / rho;
					f[DIR_MMP]  =  phi * feq[DIR_MMP]  / rho;
					f[DIR_PPM]  =  phi * feq[DIR_PPM]  / rho;
					f[DIR_MPM]  =  phi * feq[DIR_MPM]  / rho;
					f[DIR_PMM]  =  phi * feq[DIR_PMM]  / rho;
					f[DIR_MMM]  =  phi * feq[DIR_MMM]  / rho;
					f[d000] =  phi * feq[d000] / rho;

					distributionsH->setPostCollisionDistribution(f, ix1, ix2, ix3);
					distributionsH->setPreCollisionDistribution(f, ix1, ix2, ix3);

					if (distributionsH2) {

						f[dP00]    = (1.-phi) * feq[dP00] / rho;
						f[dM00]    = (1.-phi) * feq[dM00] / rho;
						f[DIR_0P0]    = (1.-phi) * feq[DIR_0P0] / rho;
						f[DIR_0M0]    = (1.-phi) * feq[DIR_0M0] / rho;
						f[DIR_00P]    = (1.-phi) * feq[DIR_00P] / rho;
						f[DIR_00M]    = (1.-phi) * feq[DIR_00M] / rho;
						f[DIR_PP0]   = (1.-phi) * feq[DIR_PP0] / rho;
						f[DIR_MM0]   = (1.-phi) * feq[DIR_MM0] / rho;
						f[DIR_PM0]   = (1.-phi) * feq[DIR_PM0] / rho;
						f[DIR_MP0]   = (1.-phi) * feq[DIR_MP0] / rho;
						f[DIR_P0P]   = (1.-phi) * feq[DIR_P0P] / rho;
						f[DIR_M0M]   = (1.-phi) * feq[DIR_M0M] / rho;
						f[DIR_P0M]   = (1.-phi) * feq[DIR_P0M] / rho;
						f[DIR_M0P]   = (1.-phi) * feq[DIR_M0P] / rho;
						f[DIR_0PP]   = (1.-phi) * feq[DIR_0PP] / rho;
						f[DIR_0MM]   = (1.-phi) * feq[DIR_0MM] / rho;
						f[DIR_0PM]   = (1.-phi) * feq[DIR_0PM] / rho;
						f[DIR_0MP]   = (1.-phi) * feq[DIR_0MP] / rho;
						f[DIR_PPP]  = (1.-phi) * feq[DIR_PPP] / rho;
						f[DIR_MPP]  = (1.-phi) * feq[DIR_MPP] / rho;
						f[DIR_PMP]  = (1.-phi) * feq[DIR_PMP] / rho;
						f[DIR_MMP]  = (1.-phi) * feq[DIR_MMP] / rho;
						f[DIR_PPM]  = (1.-phi) * feq[DIR_PPM] / rho;
						f[DIR_MPM]  = (1.-phi) * feq[DIR_MPM] / rho;
						f[DIR_PMM]  = (1.-phi) * feq[DIR_PMM] / rho;
						f[DIR_MMM]  = (1.-phi) * feq[DIR_MMM] / rho;
						f[d000] = (1.-phi) * feq[d000] / rho;

                        distributionsH2->setPostCollisionDistribution(f, ix1, ix2, ix3);
                        distributionsH2->setPreCollisionDistribution(f, ix1, ix2, ix3);                    
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
	real x1=1.0,x2=1.0,x3=1.0;
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
void MultiphaseInitDistributionsBlockVisitor::setNu( real nu )
{
	this->nu = nu;
}

