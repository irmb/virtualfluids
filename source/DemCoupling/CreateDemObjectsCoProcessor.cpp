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



CreateDemObjectsCoProcessor::CreateDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,  std::shared_ptr<Communicator> comm, SPtr<DemCoProcessor> demCoProcessor, SPtr<PhysicsEngineMaterialAdapter> demObjectMaterial) : 
   CoProcessor(grid, s),
   comm(comm),
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
      int istep = static_cast<int>(step);
      

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "CreateDemObjectsCoProcessor::process start step: " << istep);
      timer.resetAndStart();
#endif

      createGeoObjects();

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "createGeoObjects() time = "<<timer.stop()<<" s");
      if (comm->isRoot()) UBLOG(logINFO, "number of objects = "<<(int)(geoObjectPrototypeVector.size()));
      if (comm->isRoot()) UBLOG(logINFO, "CreateDemObjectsCoProcessor::process stop step: " << istep);
#endif
      
      demCoProcessor->distributeIDs();
//
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
      //timer.resetAndStart();
      std::array<double, 6> AABB ={ geoObjectPrototypeVector[i]->getX1Minimum(),geoObjectPrototypeVector[i]->getX2Minimum(),geoObjectPrototypeVector[i]->getX3Minimum(),geoObjectPrototypeVector[i]->getX1Maximum(),geoObjectPrototypeVector[i]->getX2Maximum(),geoObjectPrototypeVector[i]->getX3Maximum() };
      //UBLOG(logINFO, "demCoProcessor->isGeoObjectInAABB(AABB) = " << demCoProcessor->isGeoObjectInAABB(AABB));
      //if (demCoProcessor->isDemObjectInAABB(AABB))
      //{
      //   continue;
      //}
      SPtr<GbSphere3D> sphere = std::dynamic_pointer_cast<GbSphere3D>(geoObjectPrototypeVector[i]);
      if (demCoProcessor->isSpheresIntersection(sphere->getX1Centroid(), sphere->getX2Centroid(), sphere->getX3Centroid(), sphere->getRadius()*2.0))
      {
         continue;
      }

      SPtr<GbObject3D> geoObject((GbObject3D*)(geoObjectPrototypeVector[i]->clone()));
      SPtr<MovableObjectInteractor> geoObjectInt = SPtr<MovableObjectInteractor>(new MovableObjectInteractor(geoObject, grid, velocityBcParticleAdapter, Interactor3D::SOLID, extrapolationReconstructor, State::UNPIN));

      //SetBcBlocksBlockVisitor setBcVisitor(geoObjectInt);
      //grid->accept(setBcVisitor);

      //std::vector< std::shared_ptr<Block3D> > blockVector;
      //blockVector.clear();
      //UbTupleInt3 blockNX=grid->getBlockNX();
      //grid->getAllBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*2.0, AABB[1]-(double)val<2>(blockNX)*2.0, AABB[2]-(double)val<3>(blockNX)*2.0, AABB[3]+(double)val<1>(blockNX)*2.0, AABB[4]+(double)val<2>(blockNX)*2.0, AABB[5]+(double)val<3>(blockNX)*2.0, blockVector);
      //for (std::shared_ptr<Block3D> block : blockVector)
      //   geoObjectInt->setBCBlock(block);


      ////UBLOG(logINFO, "grid->accept(setBcVisitor) time = "<<timer.stop());
      //geoObjectInt->initInteractor();
      //UBLOG(logINFO, "geoObjectInt->initInteractor() time = "<<timer.stop());
      demCoProcessor->addInteractor(geoObjectInt, demObjectMaterial, initalVelocity[i]);
      //UBLOG(logINFO, "demCoProcessor->addInteractor() time = "<<timer.stop());
   }
}

