#ifndef SIMULATIONPARAMETERS_H 
#define SIMULATIONPARAMETERS_H

#include <boost/shared_ptr.hpp>
class SimulationParameters;
typedef boost::shared_ptr<SimulationParameters> SimulationParametersPtr;

//TODO move to LBMKernel

class SimulationParameters
{
public:
   enum CollisionModelType { UNDEFINED, INCOMPRESSIBLE, COMPRESSIBLE};
public:
   ~SimulationParameters();
   static SimulationParametersPtr getInstanz();
   CollisionModelType  getCollisionModelType();
   void setCollisionModelType(CollisionModelType collisionModel);
   bool isCompressibleModel();
   bool isIncompressibleModel();
   void setViscosity(double viscosity);
   double getViscosity();
   void setRho(double rho);
   double getRho();
   void setVelocityX(double velocity);
   double getVelocityX();
   void setVelocityY(double velocity);
   double getVelocityY();
   void setVelocityZ(double velocity);
   double getVelocityZ();
protected:
private:
   static SimulationParametersPtr instanz;
   SimulationParameters();
   //SimulationParameters(const SimulationParameters&);
   CollisionModelType _collisionModel;
   double _viscosity;
   double _rho;
   double _vx1;
   double _vx2;
   double _vx3;
};

#endif
