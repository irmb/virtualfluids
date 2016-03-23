#include "SimulationParameters.h"

SimulationParametersPtr SimulationParameters::instanz;

SimulationParameters::SimulationParameters()
{
   _collisionModel = UNDEFINED;
   _viscosity = 0.0;
   _rho = 0.0;
   _vx1 = 0.0;
   _vx3 = 0.0;
}
//////////////////////////////////////////////////////////////////////////
SimulationParameters::~SimulationParameters()
{}
//////////////////////////////////////////////////////////////////////////
SimulationParametersPtr SimulationParameters::getInstanz()
{
   if( instanz == 0 )
      instanz = SimulationParametersPtr(new SimulationParameters());
   return instanz;
}
//////////////////////////////////////////////////////////////////////////
SimulationParameters::CollisionModelType  SimulationParameters::getCollisionModelType()
{ 
   return _collisionModel;           
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setCollisionModelType(SimulationParameters::CollisionModelType collisionModel)  
{ 
   _collisionModel = collisionModel; 
}
//////////////////////////////////////////////////////////////////////////
bool SimulationParameters::isCompressibleModel()
{
   if(_collisionModel == COMPRESSIBLE) 
      return true;
   else
      return false;
}
//////////////////////////////////////////////////////////////////////////
bool SimulationParameters::isIncompressibleModel()
{
   if(_collisionModel == INCOMPRESSIBLE) 
      return true;
   else
      return false;
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setViscosity( double viscosity)
{
   _viscosity = viscosity;
}
//////////////////////////////////////////////////////////////////////////
double SimulationParameters::getViscosity()
{
   return _viscosity;
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setRho( double rho)
{
   _rho = rho;
}
//////////////////////////////////////////////////////////////////////////
double SimulationParameters::getRho()
{
   return _rho;
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setVelocityX( double velocity)
{
   _vx1 = velocity;
}
//////////////////////////////////////////////////////////////////////////
double SimulationParameters::getVelocityX()
{
   return _vx1;
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setVelocityY( double velocity)
{
   _vx2 = velocity;
}
//////////////////////////////////////////////////////////////////////////
double SimulationParameters::getVelocityY()
{
   return _vx2;
}
//////////////////////////////////////////////////////////////////////////
void SimulationParameters::setVelocityZ( double velocity)
{
   _vx3 = velocity;
}
//////////////////////////////////////////////////////////////////////////
double SimulationParameters::getVelocityZ()
{
   return _vx3;
}
