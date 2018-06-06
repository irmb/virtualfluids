#include "CreateGeoObjectsCoProcessor.h"
#include "UbScheduler.h"
#include "DemCoProcessor.h"
#include "GbSphere3D.h"
#include "VelocityBCAdapter.h"
#include "VelocityWithDensityBCAlgorithm.h"
#include "BoundaryConditionsBlockVisitor.h"
#include "MovableObjectInteractor.h"
#include "VelocityBcReconstructor.h"
#include "ExtrapolationReconstructor.h"
#include "PePhysicsEngineMaterialAdapter.h"
#include "muParser.h"
#include "PhysicsEngineMaterialAdapter.h"

CreateGeoObjectsCoProcessor::CreateGeoObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<DemCoProcessor> demCoProcessor,  
   SPtr<GbObject3D> geoObjectPrototype, SPtr<PhysicsEngineMaterialAdapter> geoObjectMaterial, SPtr<BoundaryConditionsBlockVisitor> bcVisitor) : 
   CoProcessor(grid, s),
   demCoProcessor(demCoProcessor), 
   geoObjectPrototype(geoObjectPrototype),
   geoObjectMaterial(geoObjectMaterial),
   bcVisitor(bcVisitor)
{

}
//////////////////////////////////////////////////////////////////////////
void CreateGeoObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      mu::Parser fct2;
      fct2.SetExpr("U");
      fct2.DefineConst("U", 0.0);
      SPtr<BCAdapter> velocityBcParticleAdapter(new VelocityBCAdapter(true, false, false, fct2, 0, BCFunction::INFCONST));
      velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

      SPtr<MovableObjectInteractor> sphereInt;
      const std::shared_ptr<Reconstructor> velocityReconstructor(new VelocityBcReconstructor());
      const std::shared_ptr<Reconstructor> extrapolationReconstructor(new ExtrapolationReconstructor(velocityReconstructor));
      //SPtr<GbObject3D> geoObject(dynamicPointerCast<GbSphere3D>(geoObjectPrototype)->clone());
      SPtr<GbObject3D> geoObject((GbObject3D*)(geoObjectPrototype->clone()));
      sphereInt = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(geoObject, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));
      sphereInt->setBlockVisitor(bcVisitor);

      demCoProcessor->addInteractor(sphereInt, geoObjectMaterial, Vector3D(0.0, 0.0, 0.0));
      demCoProcessor->distributeIDs();
   }
}

