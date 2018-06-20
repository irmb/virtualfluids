#ifndef CreateSphereCoProcessor_h__
#define CreateSphereCoProcessor_h__

#include "CoProcessor.h"
#include "Vector3D.h"
#include <vector>

class Grid3D;
class UbScheduler;
class DemCoProcessor;
class GbObject3D;
class PhysicsEngineMaterialAdapter;


class CreateGeoObjectsCoProcessor : public CoProcessor
{
public:
   CreateGeoObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<DemCoProcessor> demCoProcessor,  SPtr<PhysicsEngineMaterialAdapter> sphereMaterial, Vector3D initalVelocity);
   void process(double step) override;
   void addGeoObject( SPtr<GbObject3D> geoObjectPrototype);
protected:
private:
   SPtr<DemCoProcessor> demCoProcessor;
   std::vector< SPtr<GbObject3D> > geoObjectPrototypeVector;
   SPtr<PhysicsEngineMaterialAdapter> geoObjectMaterial; 
   Vector3D initalVelocity;
};
#endif // CreateSphereCoProcessor_h__
