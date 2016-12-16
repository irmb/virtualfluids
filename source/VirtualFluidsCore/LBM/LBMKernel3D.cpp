#include "LBMKernel3D.h"


LBMKernel3D::LBMKernel3D() : ghostLayerWidth(1),
                             deltaT(1.0),
                             withForcing(false),
                             withSpongeLayer(false),
                             compressible(false)
{
   this->setForcingX1(0.0);
   this->setForcingX2(0.0);
   this->setForcingX3(0.0);
   dataSet = DataSet3DPtr(new DataSet3D());
}
//////////////////////////////////////////////////////////////////////////
LBMKernel3D::~LBMKernel3D()
{

}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setBCProcessor(BCProcessorPtr bcp)
{
   bcProcessor = bcp;
}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr LBMKernel3D::getBCProcessor() 
{
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setCollisionFactor(double collFactor) 
{
   this->collFactor = collFactor;
}
//////////////////////////////////////////////////////////////////////////
double LBMKernel3D::getCollisionFactor() const
{
   return collFactor;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX1(LBMReal forcingX1)
{
    this->muForcingX1.SetExpr( UbSystem::toString(forcingX1,LBMRealLim::digits10) );  
    this->checkFunction(muForcingX1); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX2(LBMReal forcingX2)
{
   this->muForcingX2.SetExpr( UbSystem::toString(forcingX2,LBMRealLim::digits10) );  
   this->checkFunction(muForcingX2);
}
void LBMKernel3D::setForcingX3(LBMReal forcingX3)
{
   this->muForcingX3.SetExpr( UbSystem::toString(forcingX3,LBMRealLim::digits10) );  
   this->checkFunction(muForcingX3); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX1( const mu::Parser& parser)
{ 
   this->checkFunction(parser); 
   this->muForcingX1 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX2( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muForcingX2 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX3( const mu::Parser& parser)  
{ 
   this->checkFunction(parser); 
   this->muForcingX3 = parser;  
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX1( const std::string& muParserString)  
{ 
   this->muForcingX1.SetExpr(muParserString); 
   this->checkFunction(muForcingX1); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setForcingX2( const std::string& muParserString)  
{ 
   this->muForcingX2.SetExpr(muParserString); 
   this->checkFunction(muForcingX2); 
}
void LBMKernel3D::setForcingX3( const std::string& muParserString)  
{ 
   this->muForcingX3.SetExpr(muParserString); 
   this->checkFunction(muForcingX3); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::checkFunction(mu::Parser fct)
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
void LBMKernel3D::setGhostLayerWidth(int witdh)
{
   ghostLayerWidth = witdh;
}
//////////////////////////////////////////////////////////////////////////
int  LBMKernel3D::getGhostLayerWidth() const
{
   return ghostLayerWidth;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setIndex( int x1, int x2, int x3 )
{
   this->ix1 = x1;
   this->ix2 = x2;
   this->ix3 = x3;
}
//////////////////////////////////////////////////////////////////////////
DataSet3DPtr LBMKernel3D::getDataSet() const
{
   return this->dataSet;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setDeltaT( LBMReal dt )
{
   deltaT = dt;
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel3D::getCompressible() const 
{ 
   return compressible; 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setCompressible(bool val) 
{ 
   compressible = val; 
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel3D::getWithForcing() const
{
   return withForcing;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setWithForcing( bool val )
{
   withForcing = val;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setBlock( Block3DPtr block )
{
   this->block = block;
}
//////////////////////////////////////////////////////////////////////////
Block3DPtr LBMKernel3D::getBlock() const
{
   return block.lock();
}
//////////////////////////////////////////////////////////////////////////
bool LBMKernel3D::getWithSpongeLayer() const
{
   return withSpongeLayer;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setWithSpongeLayer( bool val )
{
   withSpongeLayer = val;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setSpongeLayer( const mu::Parser& parser )
{
   this->checkFunction(parser); 
   this->muSpongeLayer = parser;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setSpongeLayer( const std::string& muParserString )
{
   this->muSpongeLayer.SetExpr(muParserString); 
   this->checkFunction(muSpongeLayer); 
}
//////////////////////////////////////////////////////////////////////////
void LBMKernel3D::setDataSet(DataSet3DPtr dataSet)
{
   this->dataSet = dataSet;
}

