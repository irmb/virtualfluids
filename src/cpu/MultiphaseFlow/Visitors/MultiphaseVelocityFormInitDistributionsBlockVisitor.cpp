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
//! \file MultiphaseVelocityFormInitDistributionsBlockVisitor.cpp
//! \ingroup Visitors
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseVelocityFormInitDistributionsBlockVisitor.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "LBMKernel.h"

MultiphaseVelocityFormInitDistributionsBlockVisitor::MultiphaseVelocityFormInitDistributionsBlockVisitor() 
	: Block3DVisitor(0, D3Q27System::MAXLEVEL)
{
	this->setVx1(0.0);
	this->setVx2(0.0);
	this->setVx3(0.0);
	this->setRho(0.0);
    this->setPhi(0.0);
	this->setPressure(0.0);
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx1( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muVx1 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx2( const mu::Parser& parser)
{ 
	this->checkFunction(parser); 
	this->muVx2 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx3( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muVx3 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setRho( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muRho = parser;  
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPhi( const mu::Parser& parser)  
{ 
	this->checkFunction(parser); 
	this->muPhi = parser;  
}
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPressure(const mu::Parser& parser)
{
	this->checkFunction(parser);
	this->muPressure = parser;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx1( const std::string& muParserString)  
{ 
	this->muVx1.SetExpr(muParserString); 
	this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx2( const std::string& muParserString) 
{ 
	this->muVx2.SetExpr(muParserString); 
	this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx3( const std::string& muParserString)  
{ 
	this->muVx3.SetExpr(muParserString); 
	this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setRho( const std::string& muParserString)  
{ 
	this->muRho.SetExpr(muParserString); 
	this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPhi( const std::string& muParserString)  
{ 
	this->muPhi.SetExpr(muParserString); 
	this->checkFunction(muPhi); 
}
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPressure(const std::string& muParserString)
{
	this->muPressure.SetExpr(muParserString);
	this->checkFunction(muPressure);
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx1( real vx1 ) 
{ 
	this->muVx1.SetExpr( UbSystem::toString(vx1,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx1); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx2( real vx2 ) 
{ 
	this->muVx2.SetExpr( UbSystem::toString(vx2,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx2); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setVx3( real vx3 ) 
{ 
	this->muVx3.SetExpr( UbSystem::toString(vx3,D3Q27RealLim::digits10) );  
	this->checkFunction(muVx3); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setRho( real rho ) 
{ 
	this->muRho.SetExpr( UbSystem::toString(rho,D3Q27RealLim::digits10) );  
	this->checkFunction(muRho); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPhi( real phi ) 
{ 
	this->muPhi.SetExpr( UbSystem::toString(phi,D3Q27RealLim::digits10) );  
	this->checkFunction(muPhi); 
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityFormInitDistributionsBlockVisitor::visit(const SPtr<Grid3D> grid, SPtr<Block3D> block) 
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
	this->muPressure.DefineVar("x1", &x1); this->muPressure.DefineVar("x2", &x2); this->muPressure.DefineVar("x3", &x3);

	

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
		SPtr<PhaseFieldArray3D> pressure = dynamicPointerCast<PhaseFieldArray3D>(kernel->getDataSet()->getPressureField());


		//LBMReal phiL = kernel->getPhiL();
		//LBMReal phiH = kernel->getPhiH();

		real f[D3Q27System::ENDF+1];

		for(int ix3=0; ix3<(int)bcArray->getNX3(); ix3++)
            for (int ix2 = 0; ix2 < (int)bcArray->getNX2(); ix2++)
                for (int ix1 = 0; ix1 < (int)bcArray->getNX1(); ix1++)
				{
					Vector3D coords = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    x1              = coords[0];
                    x2              = coords[1];
                    x3              = coords[2];

					real vx1 = 0, vx2 = 0, vx3 = 0, p1 = 0, phi = 0,pres=0;
					//p1  = 0.0;
					p1 = muRho.Eval();
					vx1 = muVx1.Eval();
					vx2 = muVx2.Eval();
					vx3 = muVx3.Eval();
					phi = muPhi.Eval();
					pres = muPressure.Eval();
                    if (pressure) 
						(*pressure)(ix1, ix2, ix3) = pres;

					
					//rho = phi*1.0 + (1.0-phi)/densityRatio;
					//LBMReal rhoH = 1.0;
					//LBMReal rhoL = 1.0/densityRatio;
					//LBMReal rho = rhoH + (rhoH - rhoL)*(phi - phiH)/(phiH - phiL);

			
					real feq[27];
					real geq[27];

					//calcFeqsFct(feq,rho,vx1,vx2,vx3);
					real vx1Sq = vx1*vx1;
					real vx2Sq = vx2*vx2;
					real vx3Sq = vx3*vx3;
					for (int dir = STARTF; dir <= ENDF; dir++)
					{
						real velProd = DX1[dir]*vx1 + DX2[dir]*vx2 + DX3[dir]*vx3;
						real velSq1 = velProd*velProd;
						real gamma = WEIGTH[dir]*(3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));

						//feq[dir] = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx1Sq+vx2Sq+vx3Sq));
						feq[dir] =  WEIGTH[dir] * (1 + 3 * velProd + 4.5 * velSq1 - 1.5 * (vx1Sq + vx2Sq + vx3Sq));
						//geq[dir] = p1*WEIGTH1[dir] + gamma;
						//geq[dir] = p1*WEIGTH[dir]/(rho*UbMath::c1o3) + gamma*rho;
						geq[dir] = p1 * WEIGTH[dir] / (vf::basics::constant::c1o3) + gamma ;
					}

					f[d000] = geq[d000];
					f[dP00] = geq[dP00];
					f[dM00] = geq[dM00];
					f[d0P0] = geq[d0P0];
					f[d0M0] = geq[d0M0];
					f[d00P] = geq[d00P];
					f[d00M] = geq[d00M];
					f[dPP0] = geq[dPP0];
					f[dMM0] = geq[dMM0];
					f[dPM0] = geq[dPM0];
					f[dMP0] = geq[dMP0];
					f[dP0P] = geq[dP0P];
					f[dM0M] = geq[dM0M];
					f[dP0M] = geq[dP0M];
					f[dM0P] = geq[dM0P];
					f[d0PP] = geq[d0PP];
					f[d0MM] = geq[d0MM];
					f[d0PM] = geq[d0PM];
					f[d0MP] = geq[d0MP];
					f[dPPP] = geq[dPPP];
					f[dMPP] = geq[dMPP];
					f[dPMP] = geq[dPMP];
					f[dMMP] = geq[dMMP];
					f[dPPM] = geq[dPPM];
					f[dMPM] = geq[dMPM];
					f[dPMM] = geq[dPMM];
					f[dMMM] = geq[dMMM];
					

					distributionsF->setPostCollisionDistribution(f, ix1, ix2, ix3);
					distributionsF->setPreCollisionDistribution(f, ix1, ix2, ix3);

					f[d000] = phi * feq[d000];// / rho;
					f[dP00] = phi * feq[dP00];// / rho;
					f[dM00] = phi * feq[dM00];// / rho;
					f[d0P0] = phi * feq[d0P0];// / rho;
					f[d0M0] = phi * feq[d0M0];// / rho;
					f[d00P] = phi * feq[d00P];// / rho;
					f[d00M] = phi * feq[d00M];// / rho;
					f[dPP0] = phi * feq[dPP0];// / rho;
					f[dMM0] = phi * feq[dMM0];// / rho;
					f[dPM0] = phi * feq[dPM0];// / rho;
					f[dMP0] = phi * feq[dMP0];// / rho;
					f[dP0P] = phi * feq[dP0P];// / rho;
					f[dM0M] = phi * feq[dM0M];// / rho;
					f[dP0M] = phi * feq[dP0M];// / rho;
					f[dM0P] = phi * feq[dM0P];// / rho;
					f[d0PP] = phi * feq[d0PP];// / rho;
					f[d0MM] = phi * feq[d0MM];// / rho;
					f[d0PM] = phi * feq[d0PM];// / rho;
					f[d0MP] = phi * feq[d0MP];// / rho;
					f[dPPP] = phi * feq[dPPP];// / rho;
					f[dMPP] = phi * feq[dMPP];// / rho;
					f[dPMP] = phi * feq[dPMP];// / rho;
					f[dMMP] = phi * feq[dMMP];// / rho;
					f[dPPM] = phi * feq[dPPM];// / rho;
					f[dMPM] = phi * feq[dMPM];// / rho;
					f[dPMM] = phi * feq[dPMM];// / rho;
					f[dMMM] = phi * feq[dMMM];// / rho;
					

					distributionsH->setPostCollisionDistribution(f, ix1, ix2, ix3);
					distributionsH->setPreCollisionDistribution(f, ix1, ix2, ix3);

					if (distributionsH2) {

						f[d000] = 0;//(1. - phi) * feq[d000]; //  / rho;
						f[dP00] = 0;//(1.-phi) * feq[dP00]   ;// / rho;
						f[dM00] = 0;//(1.-phi) * feq[dM00]   ;// / rho;
						f[d0P0] = 0;//(1.-phi) * feq[d0P0]   ;// / rho;
						f[d0M0] = 0;//(1.-phi) * feq[d0M0]   ;// / rho;
						f[d00P] = 0;//(1.-phi) * feq[d00P]   ;// / rho;
						f[d00M] = 0;//(1.-phi) * feq[d00M]   ;// / rho;
						f[dPP0] = 0;//(1.-phi) * feq[dPP0]  ;// / rho;
						f[dMM0] = 0;//(1.-phi) * feq[dMM0]  ;// / rho;
						f[dPM0] = 0;//(1.-phi) * feq[dPM0]  ;// / rho;
						f[dMP0] = 0;//(1.-phi) * feq[dMP0]  ;// / rho;
						f[dP0P] = 0;//(1.-phi) * feq[dP0P]  ;// / rho;
						f[dM0M] = 0;//(1.-phi) * feq[dM0M]  ;// / rho;
						f[dP0M] = 0;//(1.-phi) * feq[dP0M]  ;// / rho;
						f[dM0P] = 0;//(1.-phi) * feq[dM0P]  ;// / rho;
						f[d0PP] = 0;//(1.-phi) * feq[d0PP]  ;// / rho;
						f[d0MM] = 0;//(1.-phi) * feq[d0MM]  ;// / rho;
						f[d0PM] = 0;//(1.-phi) * feq[d0PM]  ;// / rho;
						f[d0MP] = 0;//(1.-phi) * feq[d0MP]  ;// / rho;
						f[dPPP] = 0;//(1.-phi) * feq[dPPP] ;// / rho;
						f[dMPP] = 0;//(1.-phi) * feq[dMPP] ;// / rho;
						f[dPMP] = 0;//(1.-phi) * feq[dPMP] ;// / rho;
						f[dMMP] = 0;//(1.-phi) * feq[dMMP] ;// / rho;
						f[dPPM] = 0;//(1.-phi) * feq[dPPM] ;// / rho;
						f[dMPM] = 0;//(1.-phi) * feq[dMPM] ;// / rho;
						f[dPMM] = 0;//(1.-phi) * feq[dPMM] ;// / rho;
						f[dMMM] = 0;//(1.-phi) * feq[dMMM] ;// / rho;
						

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
void MultiphaseVelocityFormInitDistributionsBlockVisitor::checkFunction(mu::Parser fct)
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
void MultiphaseVelocityFormInitDistributionsBlockVisitor::setNu( real nu )
{
	this->nu = nu;
}

void MultiphaseVelocityFormInitDistributionsBlockVisitor::setPressure(real pres)
{
	this->muPressure.SetExpr(UbSystem::toString(pres, D3Q27RealLim::digits10));
	this->checkFunction(muPressure);

}

