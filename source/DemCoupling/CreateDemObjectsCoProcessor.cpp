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
#include "Communicator.h"



CreateDemObjectsCoProcessor::CreateDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,  std::shared_ptr<Communicator> comm, SPtr<DemCoProcessor> demCoProcessor, SPtr<PhysicsEngineMaterialAdapter> demObjectMaterial, double toleranz) : 
   CoProcessor(grid, s),
   comm(comm),
   demCoProcessor(demCoProcessor), 
   demObjectMaterial(demObjectMaterial),
   toleranz(toleranz)
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
   demCounter = 0;
}
//////////////////////////////////////////////////////////////////////////
void CreateDemObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      int istep = static_cast<int>(step);
      

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "CreateDemObjectsCoProcessor::process start step: " << istep);
      timer.resetAndStart();
#endif
      
      createGeoObjects();

#ifdef TIMING
//      if (comm->isRoot()) UBLOG(logINFO, "createGeoObjects() time = "<<timer.stop()<<" s");
//      if (comm->isRoot()) UBLOG(logINFO, "number of objects = "<<(int)(geoObjectPrototypeVector.size()));
//      if (comm->isRoot()) UBLOG(logINFO, "total number of objects = "<<demCounter);
      if (comm->isRoot()) UBLOG(logINFO, "CreateDemObjectsCoProcessor::process stop step: " << istep);
#endif
      
      //demCoProcessor->distributeIDs();

//#ifdef TIMING
//      if (comm->isRoot()) UBLOG(logINFO, "demCoProcessor->distributeIDs() time = "<<timer.stop()<<" s");
//#endif

      
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
   int size =  (int)(geoObjectPrototypeVector.size());

   std::vector< std::shared_ptr<Block3D> > blockVector;

   for (int i = 0; i < size; i++)
   {
      SPtr<GbSphere3D> sphere = std::dynamic_pointer_cast<GbSphere3D>(geoObjectPrototypeVector[i]);
      if (demCoProcessor->isSpheresIntersection(sphere->getX1Centroid(), sphere->getX2Centroid(), sphere->getX3Centroid(), sphere->getRadius()*2.0*(1.0-toleranz)))
      {
         continue;
      }

      SPtr<GbObject3D> geoObject((GbObject3D*)(geoObjectPrototypeVector[i]->clone()));
      SPtr<MovableObjectInteractor> geoObjectInt = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(geoObject, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));
      demCoProcessor->addInteractor(geoObjectInt, demObjectMaterial, initalVelocity[i]);
      demCounter++;
   }

#ifdef TIMING
   if (comm->isRoot()) UBLOG(logINFO, "createGeoObjects() time = "<<timer.stop()<<" s");
   if (comm->isRoot()) UBLOG(logINFO, "number of objects = "<<(int)(geoObjectPrototypeVector.size()));
   if (comm->isRoot()) UBLOG(logINFO, "total number of objects = "<<demCounter);
   //if (comm->isRoot()) UBLOG(logINFO, "CreateDemObjectsCoProcessor::process stop step: " << istep);
#endif

   demCoProcessor->distributeIDs();

#ifdef TIMING
   if (comm->isRoot()) UBLOG(logINFO, "demCoProcessor->distributeIDs() time = "<<timer.stop()<<" s");
#endif
}

