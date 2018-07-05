#include "CreateDemObjectsCoProcessor.h"
#include "UbScheduler.h"
#include "DemCoProcessor.h"
#include "GbSphere3D.h"
#include "NoSlipBCAlgorithm.h"
#include "VelocityBCAdapter.h"
#include "VelocityWithDensityBCAlgorithm.h"
#include "VelocityBCAlgorithm.h"
#include "MovableObjectInteractor.h"
#include "LBMReconstructor.h"
#include "EquilibriumReconstructor.h"
#include "VelocityBcReconstructor.h"
#include "ExtrapolationReconstructor.h"
#include "PePhysicsEngineMaterialAdapter.h"
#include "muParser.h"
#include "PhysicsEngineMaterialAdapter.h"
#include "SetBcBlocksBlockVisitor.h"
#include "Grid3D.h"

CreateDemObjectsCoProcessor::CreateDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<DemCoProcessor> demCoProcessor, SPtr<PhysicsEngineMaterialAdapter> demObjectMaterial) : 
   CoProcessor(grid, s),
   demCoProcessor(demCoProcessor), 
   demObjectMaterial(demObjectMaterial)
{
   mu::Parser fct;
   fct.SetExpr("U");
   fct.DefineConst("U", 0.0);
   velocityBcParticleAdapter = SPtr<BCAdapter>(new VelocityBCAdapter(true, false, false, fct, 0, BCFunction::INFCONST));
   velocityBcParticleAdapter->setBcAlgorithm(SPtr<BCAlgorithm>(new VelocityWithDensityBCAlgorithm()));

   //const std::shared_ptr<Reconstructor> velocityReconstructor(new VelocityBcReconstructor());
   std::shared_ptr<Reconstructor> equilibriumReconstructor(new EquilibriumReconstructor());
   //const std::shared_ptr<Reconstructor> lbmReconstructor(new LBMReconstructor(false));
   extrapolationReconstructor = SPtr<Reconstructor>(new ExtrapolationReconstructor(equilibriumReconstructor));
}
//////////////////////////////////////////////////////////////////////////
void CreateDemObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      createGeoObjects();
   }
}
//////////////////////////////////////////////////////////////////////////
void CreateDemObjectsCoProcessor::addGeoObject(SPtr<GbObject3D> geoObjectPrototype,  Vector3D  initalVelocity)
{
   geoObjectPrototypeVector.push_back(geoObjectPrototype);
   this->initalVelocity.push_back(initalVelocity);
}

void CreateDemObjectsCoProcessor::clearGeoObjects()
{
   geoObjectPrototypeVector.clear();
   initalVelocity.clear();
}

void CreateDemObjectsCoProcessor::createGeoObjects()
{
   int size =  geoObjectPrototypeVector.size();

   for (int i = 0; i < size; i++)
   {
      std::array<double, 6> AABB ={ geoObjectPrototypeVector[i]->getX1Minimum(),geoObjectPrototypeVector[i]->getX2Minimum(),geoObjectPrototypeVector[i]->getX3Minimum(),geoObjectPrototypeVector[i]->getX1Maximum(),geoObjectPrototypeVector[i]->getX2Maximum(),geoObjectPrototypeVector[i]->getX3Maximum() };
      //UBLOG(logINFO, "demCoProcessor->isGeoObjectInAABB(AABB) = " << demCoProcessor->isGeoObjectInAABB(AABB));
      if (demCoProcessor->isDemObjectInAABB(AABB))
      {
         continue;
      }
      SPtr<GbObject3D> geoObject((GbObject3D*)(geoObjectPrototypeVector[i]->clone()));
      SPtr<MovableObjectInteractor> geoObjectInt = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(geoObject, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));
      SetBcBlocksBlockVisitor setBcVisitor(geoObjectInt);
      grid->accept(setBcVisitor);
      geoObjectInt->initInteractor();
      demCoProcessor->addInteractor(geoObjectInt, demObjectMaterial, initalVelocity[i]);
   }
}

