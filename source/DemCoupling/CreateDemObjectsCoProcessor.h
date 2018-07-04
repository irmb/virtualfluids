#ifndef CreateSphereCoProcessor_h__
#define CreateSphereCoProcessor_h__

#include "CoProcessor.h"
#include "Vector3D.h"
#include <vector>
#include <array>

class Grid3D;
class UbScheduler;
class DemCoProcessor;
class GbObject3D;
class BCAdapter;
class Reconstructor;
class PhysicsEngineMaterialAdapter;


class CreateDemObjectsCoProcessor : public CoProcessor
{
public:
   CreateDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<DemCoProcessor> demCoProcessor,  SPtr<PhysicsEngineMaterialAdapter> geoObjectMaterial);
   void process(double step) override;
   void addGeoObject(SPtr<GbObject3D> geoObjectPrototype, Vector3D  initalVelocity);
   void clearGeoObjects();
   void createGeoObjects();
protected:
private:
   SPtr<DemCoProcessor> demCoProcessor;
   std::vector< SPtr<GbObject3D> > geoObjectPrototypeVector;
   SPtr<PhysicsEngineMaterialAdapter> demObjectMaterial; 
   std::vector<Vector3D>  initalVelocity;
   SPtr<BCAdapter> velocityBcParticleAdapter;
   SPtr<Reconstructor> extrapolationReconstructor;
};
#endif // CreateSphereCoProcessor_h__
