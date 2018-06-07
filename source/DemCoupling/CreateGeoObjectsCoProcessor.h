#ifndef CreateSphereCoProcessor_h__
#define CreateSphereCoProcessor_h__

#include "CoProcessor.h"

class Grid3D;
class UbScheduler;
class DemCoProcessor;
class GbObject3D;
class PhysicsEngineMaterialAdapter;


class CreateGeoObjectsCoProcessor : public CoProcessor
{
public:
   CreateGeoObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, 
      SPtr<DemCoProcessor> demCoProcessor,  SPtr<GbObject3D> spherePrototype, 
      SPtr<PhysicsEngineMaterialAdapter> sphereMaterial);
   void process(double step) override;
protected:
private:
   SPtr<DemCoProcessor> demCoProcessor;
   SPtr<GbObject3D> geoObjectPrototype;
   SPtr<PhysicsEngineMaterialAdapter> geoObjectMaterial; 
};
#endif // CreateSphereCoProcessor_h__
