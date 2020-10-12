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
//! \file InitDistributionsBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher, Soeren Freudiger
//=======================================================================================

#ifndef InitDistributionsBlockVisitor_H
#define InitDistributionsBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"

#include <muParser.h>

class Grid3D;
class Block3D;

//! \brief A class implements an initialization of the flow area.
//! \details 
//! It is more flexible way to initialize flow area.
//! You can define functions to calculate macroscopic values for feq. 
//! x1,x2,x3 are automatically defined via this adapter and are the real world
//! vertex coordinates.
//!
//!if function is invalid an UbException with detailed information is thrown
//!
//! Example:
//! \code
//! InitDistributionsBlockVisitor init;
//! init.setVx1("0.01*x2");
//! init.setVx2("0.01*x2^2");
//! \endcode

class InitDistributionsBlockVisitor : public Block3DVisitor
{
public:
   using D3Q27RealLim = std::numeric_limits<LBMReal>;

public:
   InitDistributionsBlockVisitor();
   //////////////////////////////////////////////////////////////////////////
   //automatic vars are: x1,x2, x3
   //ussage example: setVx1("x1*0.01+x2*0.003")
   //////////////////////////////////////////////////////////////////////////
   void setVx1( const mu::Parser& parser);
   void setVx2( const mu::Parser& parser);
   void setVx3( const mu::Parser& parser);
   void setRho( const mu::Parser& parser);

   void setVx1( const std::string& muParserString);
   void setVx2( const std::string& muParserString);
   void setVx3( const std::string& muParserString);
   void setRho( const std::string& muParserString);
   //////////////////////////////////////////////////////////////////////////
   void setVx1( LBMReal vx1 );
   void setVx2( LBMReal vx2 );
   void setVx3( LBMReal vx3 );
   void setRho( LBMReal rho );

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

protected:
   void checkFunction(mu::Parser fct);

private:
   mu::Parser muVx1;
   mu::Parser muVx2;
   mu::Parser muVx3;
   mu::Parser muRho;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H
