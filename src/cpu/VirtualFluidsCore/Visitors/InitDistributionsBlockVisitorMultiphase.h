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
//! \file InitDistributionsBlockVisitorMultiphase.h
//! \ingroup Visitors
//! \author Hesameddin Safari
//=======================================================================================

#ifndef InitDistributionsBlockVisitorMultiphase_H
#define InitDistributionsBlockVisitorMultiphase_H

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Block3D.h"

#include <muParser.h>



class InitDistributionsBlockVisitorMultiphase : public Block3DVisitor
{
public:
	typedef std::numeric_limits<LBMReal> D3Q27RealLim;

public:
	InitDistributionsBlockVisitorMultiphase();
	//D3Q27ETInitDistributionsBlockVisitor(LBMReal rho, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
	//! Constructor
	//! \param nu - viscosity
	//! \param rho - density
	//! \param vx1 - velocity in x
	//! \param vx2 - velocity in y
	//! \param vx3 - velocity in z
	InitDistributionsBlockVisitorMultiphase( LBMReal densityRatio, LBMReal intThickness, LBMReal radius, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
	//////////////////////////////////////////////////////////////////////////
	//automatic vars are: x1,x2, x3
	//ussage example: setVx1("x1*0.01+x2*0.003")
	//////////////////////////////////////////////////////////////////////////
	void setVx1( const mu::Parser& parser);
	void setVx2( const mu::Parser& parser);
	void setVx3( const mu::Parser& parser);
	void setRho( const mu::Parser& parser);
	void setPhi( const mu::Parser& parser);

	void setVx1( const std::string& muParserString);
	void setVx2( const std::string& muParserString);
	void setVx3( const std::string& muParserString);
	void setRho( const std::string& muParserString);
	void setPhi( const std::string& muParserString);

	//////////////////////////////////////////////////////////////////////////
	void setVx1( LBMReal vx1 );
	void setVx2( LBMReal vx2 );
	void setVx3( LBMReal vx3 );
	void setRho( LBMReal rho );
	void setPhi( LBMReal rho );
	void setNu( LBMReal nu );

	void visit(SPtr<Grid3D> grid, SPtr<Block3D> block);

protected:
	void checkFunction(mu::Parser fct);

private:
	mu::Parser muVx1;
	mu::Parser muVx2;
	mu::Parser muVx3;
	mu::Parser muRho;
	mu::Parser muPhi;

	LBMReal nu;
	LBMReal densityRatio;
	LBMReal intThickness;
	LBMReal radius;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H
