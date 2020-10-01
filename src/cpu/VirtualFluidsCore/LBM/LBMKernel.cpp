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
//! \file LBMKernel.cpp
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================

#include "LBMKernel.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "BCArray3D.h"

LBMKernel::LBMKernel() : ghostLayerWidth(1),
                             deltaT(1.0),
                             withForcing(false),
                             withSpongeLayer(false),
                             compressible(false)
{
   this->setForcingX1(0.0);
   this->setForcingX2(0.0);
   this->setForcingX3(0.0);
   dataSet = SPtr<DataSet3D>(new DataSet3D());
   this->nx[0] = 0;
   this->nx[1] = 0;
   this->nx[2] = 0;
}
//////////////////////////////////////////////////////////////////////////
LBMKernel::~LBMKernel()
{

}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setBCProcessor(SPtr<BCProcessor> bcp)
{
   bcProcessor = bcp;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCProcessor> LBMKernel::getBCProcessor() const
{
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setCollisionFactor(double collFactor) 
{
   this->collFactor = collFactor;
}
//////////////////////////////////////////////////////////////////////////
double LBMKernel::getCollisionFactor() const
{
   return collFactor;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1(LBMReal forcingX1)
{
    this->muForcingX1.SetExpr( UbSystem::toString(forcingX1,LBMRealLim::digits10) );  
    this->checkFunction(muForcingX1); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2(LBMReal forcingX2)
{
   this->muForcingX2.SetExpr( UbSystem::toString(forcingX2,LBMRealLim::digits10) );  
   this->checkFunction(muForcingX2);
}
void LBMKernel::setForcingX3(LBMReal forcingX3)
{
   this->muForcingX3.SetExpr( UbSystem::toString(forcingX3,LBMRealLim::digits10) );  
   this->checkFunction(muForcingX3); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1( const mu::Parser& parser)
{ 
   this->checkFunction(parser); 
   this->muForcingX1 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muForcingX2 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX3( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muForcingX3 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX1( const std::string& muParserString)  
{ 
   this->muForcingX1.SetExpr(muParserString); 
   this->checkFunction(muForcingX1); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setForcingX2( const std::string& muParserString)  
{ 
   this->muForcingX2.SetExpr(muParserString); 
   this->checkFunction(muForcingX2); 
}
void LBMKernel::setForcingX3( const std::string& muParserString)  
{ 
   this->muForcingX3.SetExpr(muParserString); 
   this->checkFunction(muForcingX3); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::checkFunction(mu::Parser fct)
{
   double x1=1.0,x2=1.0,x3=1.0, dt=1.0, nue=1.0;
   fct.DefineVar("x1",&x1); 
   fct.DefineVar("x2",&x2); 
   fct.DefineVar("x3",&x3);
   fct.DefineVar("dt",&dt);
   fct.DefineVar("nue",&nue);

   try
   {
      fct.Eval();
      fct.ClearVar();
   }
   catch(mu::ParserError& e)
   {
      throw UbException(UB_EXARGS,"function: "+e.GetExpr() + (std::string)"error: "+e.GetMsg()
         +(std::string)", only x1,x2,x3,dx are allowed as variables" );
   }
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setGhostLayerWidth(int witdh)
{
   ghostLayerWidth = witdh;
}
//////////////////////////////////////////////////////////////////////////
int  LBMKernel::getGhostLayerWidth() const
{
   return ghostLayerWidth;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setIndex( int x1, int x2, int x3 )
{
   this->ix1 = x1;
   this->ix2 = x2;
   this->ix3 = x3;
}
//////////////////////////////////////////////////////////////////////////
SPtr<DataSet3D> LBMKernel::getDataSet() const
{
   return this->dataSet;
}
//////////////////////////////////////////////////////////////////////////
LBMReal LBMKernel::getDeltaT() const
{
   return this->deltaT;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setDeltaT( LBMReal dt )
{
   deltaT = dt;
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getCompressible() const 
{ 
   return compressible; 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setCompressible(bool val) 
{ 
   compressible = val; 
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getWithForcing() const
{
   return withForcing;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setWithForcing( bool val )
{
   withForcing = val;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setBlock( SPtr<Block3D> block )
{
   this->block = block;
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> LBMKernel::getBlock() const
{
   return block.lock();
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::getWithSpongeLayer() const
{
   return withSpongeLayer;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setWithSpongeLayer( bool val )
{
   withSpongeLayer = val;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setSpongeLayer( const mu::Parser& parser )
{
   this->checkFunction(parser); 
   this->muSpongeLayer = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setSpongeLayer( const std::string& muParserString )
{
   this->muSpongeLayer.SetExpr(muParserString); 
   this->checkFunction(muSpongeLayer); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setDataSet(SPtr<DataSet3D> dataSet)
{
   this->dataSet = dataSet;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::swapDistributions()
{
   dataSet->getFdistributions()->swap();
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel::setNX(std::array<int, 3> nx)
{
   this->nx = nx;
}
//////////////////////////////////////////////////////////////////////////
std::array<int, 3> LBMKernel::getNX()
{
   return nx;
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel::isInsideOfDomain(const int& x1, const int& x2, const int& x3) const
{
    const SPtr<BCArray3D> bcArray = this->bcProcessor->getBCArray();
    return bcArray->isInsideOfDomain(x1, x2, x3, ghostLayerWidth);
}
