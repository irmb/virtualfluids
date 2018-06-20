#include "CreateGeoObjectsCoProcessor.h"
#include "UbScheduler.h"
#include "DemCoProcessor.h"
#include "GbSphere3D.h"
#include "NoSlipBCAlgorithm.h"
#include "VelocityBCAdapter.h"
#include "VelocityWithDensityBCAlgorithm.h"
#include "VelocityBCAlgorithm.h"
#include "MovableObjectInteractor.h"
#include "EquilibriumReconstructor.h"
#include "VelocityBcReconstructor.h"
#include "ExtrapolationReconstructor.h"
#include "PePhysicsEngineMaterialAdapter.h"
#include "muParser.h"
#include "PhysicsEngineMaterialAdapter.h"
#include "SetSolidOrBoundaryBlockVisitor.h"
#include "Grid3D.h"

CreateGeoObjectsCoProcessor::CreateGeoObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<DemCoProcessor> demCoProcessor, SPtr<PhysicsEngineMaterialAdapter> geoObjectMaterial, Vector3D initalVelocity) : 
   CoProcessor(grid, s),
   demCoProcessor(demCoProcessor), 
   geoObjectMaterial(geoObjectMaterial),
   initalVelocity(initalVelocity)
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
      //velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new NoSlipBCAlgorithm()));

      SPtr<MovableObjectInteractor> geoObjectInt;
      //const std::shared_ptr<Reconstructor> velocityReconstructor(new VelocityBcReconstructor());
      const std::shared_ptr<Reconstructor> equilibriumReconstructor(new EquilibriumReconstructor());
      const std::shared_ptr<Reconstructor> extrapolationReconstructor(new ExtrapolationReconstructor(equilibriumReconstructor));
      
      for(SPtr<GbObject3D> proto : geoObjectPrototypeVector)
      {
         SPtr<GbObject3D> geoObject((GbObject3D*)(proto->clone()));
         geoObjectInt = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(geoObject, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));
         //SetSolidOrBoundaryBlockVisitor setSolidVisitor(geoObjectInt, SetSolidOrBoundaryBlockVisitor::SOLID);
         //grid->accept(setSolidVisitor);
         SetSolidOrBoundaryBlockVisitor setBcVisitor(geoObjectInt, SetSolidOrBoundaryBlockVisitor::BC);
         grid->accept(setBcVisitor);
         geoObjectInt->initInteractor();
         demCoProcessor->addInteractor(geoObjectInt, geoObjectMaterial, initalVelocity);
         //demCoProcessor->distributeIDs();
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void CreateGeoObjectsCoProcessor::addGeoObject(SPtr<GbObject3D> geoObjectPrototype)
{
   geoObjectPrototypeVector.push_back(geoObjectPrototype);
}

