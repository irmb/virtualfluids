#include "DemCoProcessor.h"

#include "GbSphere3D.h"
#include "MovableObjectInteractor.h"
#include "Communicator.h"
#include "ForceCalculator.h"
#include "Grid3D.h"
#include "UbScheduler.h"
#include "ILBMKernel.h"
#include "DistributionArray3D.h"
#include "BCProcessor.h"
#include "DataSet3D.h"
#include "SetBcBlocksBlockVisitor.h"

#include "PhysicsEngineMaterialAdapter.h"
#include "PhysicsEngineGeometryAdapter.h"
#include "PhysicsEngineSolverAdapter.h"
#include "PePhysicsEngineSolverAdapter.h"
#include "PePhysicsEngineGeometryAdapter.h"

#include "BoundaryConditions.h"
#include "Block3D.h"
#include "BCArray3D.h"
#include "MPICommunicator.h"
#include "BoundaryConditionsBlockVisitor.h"

#include "UbLogger.h"


#include <array>
#include <functional>

DemCoProcessor::DemCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Communicator> comm, std::shared_ptr<ForceCalculator> forceCalculator, std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps) :
   CoProcessor(grid, s), comm(comm), forceCalculator(forceCalculator), physicsEngineSolver(physicsEngineSolver), intermediateDemSteps(intermediatePeSteps)
{
#ifdef TIMING
   timer.resetAndStart();
#endif

   std::shared_ptr<walberla::blockforest::BlockForest> forest = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
   std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();


   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      walberla::pe::Storage* storage = blockIt->getData< walberla::pe::Storage >(*storageId.get());
      walberla::pe::BodyStorage* bodyStorage = &(*storage)[0];
      walberla::pe::BodyStorage* bodyStorageShadowCopies = &(*storage)[1];

      bodyStorage->registerAddCallback("DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::addPeGeo), this));
      //bodyStorage->registerRemoveCallback("DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::removePeGeo), this));

      bodyStorageShadowCopies->registerAddCallback("DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::addPeShadowGeo), this));
      bodyStorageShadowCopies->registerRemoveCallback("DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::removePeShadowGeo), this));
   }
}

DemCoProcessor::~DemCoProcessor()
{
   std::shared_ptr<walberla::blockforest::BlockForest> forest = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
   std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();

   for (auto& currentBlock : *forest)
   {
      walberla::pe::Storage * storage = currentBlock.getData< walberla::pe::Storage >(*storageId.get());
      walberla::pe::BodyStorage& localStorage = (*storage)[0];
      walberla::pe::BodyStorage& shadowStorage = (*storage)[1];

      localStorage.clearAddCallbacks();
      localStorage.clearRemoveCallbacks();

      shadowStorage.clearAddCallbacks();
      shadowStorage.clearRemoveCallbacks();

      //bodyStorage->deregisterAddCallback("DemCoProcessor");
      ////bodyStorage->deregisterRemoveCallback("DemCoProcessor");

      //if (&bodyStorage != &bodyStorageShadowCopies)
      //{
      //   bodyStorageShadowCopies->deregisterAddCallback("DemCoProcessor");
      //   bodyStorageShadowCopies->deregisterRemoveCallback("DemCoProcessor");
      //}
   }
}

void DemCoProcessor::addInteractor(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial, Vector3D initalVelocity)
{
   interactors.push_back(interactor);
   const int id = static_cast<int>(interactors.size()-1);

   //if (comm->getProcessID()==2) UBLOG(logINFO, "DemCoProcessor::addInteractor() id = "<<id);
   interactor->setID(id);
   const auto peGeometryAdapter = this->createPhysicsEngineGeometryAdapter(interactor, physicsEngineMaterial);
   if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometryAdapter)->isActive())
   {
      peGeometryAdapter->setLinearVelolocity(initalVelocity);
   }
   SetBcBlocksBlockVisitor setBcVisitor(interactor);
   grid->accept(setBcVisitor);

   //std::vector< std::shared_ptr<Block3D> > blockVector;
   //UbTupleInt3 blockNX=grid->getBlockNX();
   //SPtr<GbObject3D> geoObject(interactor->getGbObject3D());
   //double ext = 0.0;
   //std::array<double, 6> AABB ={ geoObject->getX1Minimum(),geoObject->getX2Minimum(),geoObject->getX3Minimum(),geoObject->getX1Maximum(),geoObject->getX2Maximum(),geoObject->getX3Maximum() };
   //grid->getBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*ext, AABB[1]-(double)val<2>(blockNX)*ext, AABB[2]-(double)val<3>(blockNX)*ext, AABB[3]+(double)val<1>(blockNX)*ext, AABB[4]+(double)val<2>(blockNX)*ext, AABB[5]+(double)val<3>(blockNX)*ext, blockVector);
   //for (std::shared_ptr<Block3D> block : blockVector)
   //{
   //   if (block->getKernel())
   //   {
   //      interactor->setBCBlock(block);
   //      //UBLOG(logINFO, "DemCoProcessor::addInteractor() rank = "<<comm->getProcessID());
   //   }
   //}

   interactor->initInteractor();

   physicsEngineGeometrieAdapters.push_back(peGeometryAdapter);
}


std::shared_ptr<PhysicsEngineGeometryAdapter> DemCoProcessor::createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const
{
   const int id = static_cast<int>(interactors.size()-1);
   SPtr<GbSphere3D> vfSphere = std::static_pointer_cast<GbSphere3D>(interactor->getGbObject3D());
   const Vector3D position(vfSphere->getX1Centroid(), vfSphere->getX2Centroid(), vfSphere->getX3Centroid());

   auto peGeometryAdapter = this->physicsEngineSolver->createPhysicsEngineGeometryAdapter(id, position, vfSphere->getRadius(), physicsEngineMaterial);
   //if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometry)->isActive())
   //{
   interactor->setPhysicsEngineGeometry(peGeometryAdapter);
   return peGeometryAdapter;
   //}
   //else
   //{
   //   return peGeometry;
   //}
}


void DemCoProcessor::process(double actualTimeStep)
{
#ifdef TIMING
   timer.resetAndStart();
#endif

   this->applyForcesOnGeometries();

#ifdef TIMING
   if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::process start step: " << actualTimeStep);
   if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::applyForcesOnGeometries() time = "<<timer.stop()<<" s");
#endif

   if (scheduler->isDue(actualTimeStep))
   {
      //UBLOG(logINFO, "DemCoProcessor::update - START - timestep = " << actualTimeStep);
      const double demTimeStepsPerIteration = scheduler->getMinStep();

      if (demTimeStepsPerIteration != 1)
         this->scaleForcesAndTorques(1.0 / demTimeStepsPerIteration);

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::scaleForcesAndTorques() time = "<<timer.stop()<<" s");
      if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::calculateDemTimeStep():");
#endif

      if (this->intermediateDemSteps == 1)
         this->calculateDemTimeStep(demTimeStepsPerIteration);

      //#ifdef TIMING
      //      if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::calculateDemTimeStep() time = "<<timer.stop()<<" s");
      //#endif
            //if ((int)actualTimeStep % 100 == 0)
            //{
            //    if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[0])->isActive())
            //    {
            //        //UBLOG(logINFO, "v: (x,y,z) " << physicsEngineGeometries[0]->getLinearVelocity() << " actualTimeStep = " << UbSystem::toString(actualTimeStep));
            //    }
            //}

            // during the intermediate time steps of the collision response, the currently acting forces
            // (interaction forces, gravitational force, ...) have to remain constant.
            // Since they are reset after the call to collision response, they have to be stored explicitly before.
            // Then they are set again after each intermediate step.

      this->moveVfGeoObjects();

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::moveVfGeoObject() time = "<<timer.stop()<<" s");
#endif

      grid->accept(*boundaryConditionsBlockVisitor.get());

#ifdef TIMING
      if (comm->isRoot()) UBLOG(logINFO, "grid->accept(*boundaryConditionsBlockVisitor.get()) time = "<<timer.stop()<<" s");
#endif

      //UBLOG(logINFO, "DemCoProcessor::update - END - timestep = " << actualTimeStep);
   }

#ifdef TIMING
   if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::process stop step: " << actualTimeStep);
#endif
}
//////////////////////////////////////////////////////////////////////////
std::shared_ptr<PhysicsEngineSolverAdapter> DemCoProcessor::getPhysicsEngineSolver()
{
   return physicsEngineSolver;
}

void DemCoProcessor::applyForcesOnGeometries()
{
   for (int i = 0; i < physicsEngineGeometrieAdapters.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         this->setForcesToObject(grid, interactors[i], physicsEngineGeometrieAdapters[i]);

         //physicsEngineGeometries[i]->setLinearVelolocity(Vector3D(-0.001, 0.0, 0.0));
         //physicsEngineGeometries[i]->setAngularVelocity(Vector3D(0.01, 0.01, 0.01));
         //UBLOG(logINFO, "v: (x,y,z) " << physicsEngineGeometries[i]->getLinearVelocity());
      }
   }
}

void DemCoProcessor::setForcesToObject(SPtr<Grid3D> grid, SPtr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry)
{
   for (BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap())
   {
      SPtr<Block3D> block = t.first;
      SPtr<ILBMKernel> kernel = block->getKernel();
      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();
      SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
      distributions->swap();

      std::set< std::vector<int> >& transNodeIndicesSet = t.second;
      for (std::vector<int> node : transNodeIndicesSet)
      {
         int x1 = node[0];
         int x2 = node[1];
         int x3 = node[2];

         if (kernel->isInsideOfDomain(x1, x2, x3) && bcArray->isFluid(x1, x2, x3))
         {
            //TODO: calculate assumed boundary position 

            const Vector3D worldCoordinates = grid->getNodeCoordinates(block, x1, x2, x3);
            const auto boundaryVelocity = physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);

            SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
            const Vector3D force = forceCalculator->getForces(x1, x2, x3, distributions, bc, boundaryVelocity);
            physicsEngineGeometry->addForceAtPosition(force, worldCoordinates);
         }
      }
      distributions->swap();
   }
}


void DemCoProcessor::scaleForcesAndTorques(double scalingFactor)
{
   for (int i = 0; i < physicsEngineGeometrieAdapters.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         const Vector3D force = physicsEngineGeometrieAdapters[i]->getForce() * scalingFactor;
         const Vector3D torque = physicsEngineGeometrieAdapters[i]->getTorque() * scalingFactor;

         physicsEngineGeometrieAdapters[i]->resetForceAndTorque();

         physicsEngineGeometrieAdapters[i]->setForce(force);
         physicsEngineGeometrieAdapters[i]->setTorque(torque);

         //UBLOG(logINFO, "F: (x,y,z) " << force);
         //UBLOG(logINFO, "T: (x,y,z) " << torque);
      }
   }
}


void DemCoProcessor::calculateDemTimeStep(double step)
{
   physicsEngineSolver->runTimestep(step);

#ifdef TIMING
   if (comm->isRoot()) UBLOG(logINFO, "  physicsEngineSolver->runTimestep() time = "<< timer.stop() <<" s");
#endif

   //for (int i = 0; i < physicsEngineGeometries.size(); i++)
   //{
   //   if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
   //   {
   //      interactors[i]->setPhysicsEngineGeometry(physicsEngineGeometries[i]);
   //   }
   //}

//#ifdef TIMING
//   if (comm->isRoot()) UBLOG(logINFO, "  physicsEngineSolver->updateGeometry() time = "<<timer.stop()<<" s");
//#endif
}

void DemCoProcessor::moveVfGeoObjects()
{
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         walberla::pe::RigidBody* peGeoObject = getPeGeoObject(std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->getSystemId());
         if (peGeoObject != nullptr)
         {
            std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->setGeometry(peGeoObject);
            interactors[i]->moveGbObjectTo(physicsEngineGeometrieAdapters[i]->getPosition());
         }
         else
         {
            std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->setInactive();
         }

         //UBLOG(logINFO, "DemCoProcessor::moveVfGeoObject() id = "<<std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->getId()<<"  position="<<physicsEngineGeometrieAdapters[i]->getPosition()<<" rank="<<comm->getProcessID());
      }
   }
}

bool  DemCoProcessor::isDemObjectInAABB(std::array<double, 6> AABB)
{
   bool result = false;
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         SPtr<GbObject3D> geoObject = interactors[i]->getGbObject3D();
         std::array <double, 2> minMax1;
         std::array <double, 2> minMax2;
         std::array <double, 2> minMax3;
         minMax1[0] = geoObject->getX1Minimum();
         minMax2[0] = geoObject->getX2Minimum();
         minMax3[0] = geoObject->getX3Minimum();
         minMax1[1] = geoObject->getX1Maximum();
         minMax2[1] = geoObject->getX2Maximum();
         minMax3[1] = geoObject->getX3Maximum();

         for (int x3 = 0; x3 < 2; x3++)
            for (int x2 = 0; x2 < 2; x2++)
               for (int x1 = 0; x1 < 2; x1++)
               {
                  result = result || (minMax1[x1] >= AABB[0] && minMax2[x2] >= AABB[1] && minMax3[x3] >= AABB[2] && minMax1[x1] <= AABB[3] && minMax2[x2] <= AABB[4] && minMax3[x3] <= AABB[5]);
               }
      }
   }

   std::vector<int> values;
   values.push_back((int)result);
   std::vector<int> rvalues = comm->gather(values);

   if (comm->isRoot())
   {
      for (int i = 0; i < (int)rvalues.size(); i++)
      {
         result = result || (bool)rvalues[i];
      }
   }
   int iresult = (int)result;
   comm->broadcast(iresult);
   result = (bool)iresult;

   return result;
}

void DemCoProcessor::addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles)
{
   for (int i = 0; i < interactors.size(); i++)
   {
      //UBLOG(logINFO, "DemCoProcessor::addSurfaceTriangleSet()1");
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         //UBLOG(logINFO, "DemCoProcessor::addSurfaceTriangleSet()2");
         interactors[i]->getGbObject3D()->addSurfaceTriangleSet(nodes, triangles);
      }
   }
}

void DemCoProcessor::getObjectsPropertiesVector(std::vector<double>& p)
{
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         p.push_back(interactors[i]->getGbObject3D()->getX1Centroid());
         p.push_back(interactors[i]->getGbObject3D()->getX2Centroid());
         p.push_back(interactors[i]->getGbObject3D()->getX3Centroid());
         Vector3D v = physicsEngineGeometrieAdapters[i]->getLinearVelocity();
         p.push_back(v[0]);
         p.push_back(v[1]);
         p.push_back(v[2]);
      }
   }
}

void DemCoProcessor::addPeGeo(walberla::pe::RigidBody * peGeo)
{



   if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;

   auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   geometry->counter++;
   //if (comm->getProcessID()==3)UBLOG(logINFO, "DemCoProcessor::addPeGeo() geoID = "<<peGeo->getID()<<" counter="<<geometry->counter<<" shadowCounter="<<geometry->shadowCounter);
   if (geometry->getSystemId() == peGeo->getSystemID())
   {
      geometry->setActive();
      geometry->setGeometry(peGeo);
      return;
   }
   else
      throw UbException(UB_EXARGS, "PeGeo ID="+UbSystem::toString(peGeo->getID())+" is not matching geometry ID="+UbSystem::toString(geometry->getId()));
}

void DemCoProcessor::removePeGeo(walberla::pe::RigidBody * peGeo)
{
   //if (comm->getProcessID()==2) UBLOG(logINFO, "DemCoProcessor::removePeGeo() rank = "<<comm->getProcessID());

   //if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;


   auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   geometry->counter--;

   if (geometry->getSystemId() == peGeo->getSystemID())
   {
      if (geometry->shadowCounter == 0)
      {
         geometry->setInactive();
         geometry->setGeometry(NULL);
      }
      //if (comm->getProcessID()==3)UBLOG(logINFO, "DemCoProcessor::removePeGeo() geoID = "<<peGeo->getID()<<" counter="<<geometry->counter<<" shadowCounter="<<geometry->shadowCounter);
      return;
   }
   else
      throw UbException(UB_EXARGS, "PeGeo ID="+UbSystem::toString(peGeo->getID())+" is not matching geometry ID="+UbSystem::toString(geometry->getId()));


   //if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;

   //std::shared_ptr<walberla::blockforest::BlockForest> forest = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
   //std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();

   //for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   //{
   //   for (auto bodyIt = walberla::pe::BodyIterator::begin(*blockIt, *storageId.get()); bodyIt != walberla::pe::BodyIterator::end(); ++bodyIt)
   //   {
   //      if (bodyIt->getID() == peGeo->getID())
   //      {
   //         return;
   //      }
   //   }
   //}

   //auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   //if (geometry->getId() == peGeo->getID())
   //{
   //   geometry->setInactive();
   //   geometry->setGeometry(NULL);
   //   return;
   //}
   //else
   //   throw UbException(UB_EXARGS, "PeGeo ID is not matching!");
}

void DemCoProcessor::addPeShadowGeo(walberla::pe::RigidBody * peGeo)
{
   if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;

   auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   geometry->shadowCounter++;
   //if (comm->getProcessID()==3) UBLOG(logINFO, "DemCoProcessor::addPeShadowGeo() geoID = "<<peGeo->getID()<<" counter="<<geometry->counter<<" shadowCounter="<<geometry->shadowCounter);
   if (geometry->getSystemId() == peGeo->getSystemID())
   {
      if (geometry->shadowCounter == 1)
      {
         geometry->setActive();
         geometry->setGeometry(peGeo);
      }

      return;
   }
   else
      throw UbException(UB_EXARGS, "PeGeo ID="+UbSystem::toString(peGeo->getID())+" is not matching geometry ID="+UbSystem::toString(geometry->getId()));
}

void DemCoProcessor::removePeShadowGeo(walberla::pe::RigidBody * peGeo)
{
   //if (comm->getProcessID()==2) UBLOG(logINFO, "DemCoProcessor::removePeShadowGeo() rank = "<<comm->getProcessID());

   //if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;

   //std::shared_ptr<walberla::blockforest::BlockForest> forest = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
   //std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();

   //for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   //{
   //   for (auto bodyIt = walberla::pe::BodyIterator::begin(*blockIt, *storageId.get()); bodyIt != walberla::pe::BodyIterator::end(); ++bodyIt)
   //   {
   //      if (bodyIt->getID() == peGeo->getID())
   //      {
   //         return;
   //      }
   //   }
   //}



   //std::vector< std::shared_ptr<Block3D> > blockVector;
   //UbTupleInt3 blockNX=grid->getBlockNX();
   //SPtr<GbObject3D> geoObject(interactors[peGeo->getID()]->getGbObject3D());
   //double ext = 0.0;
   //std::array<double, 6> AABB ={ geoObject->getX1Minimum(),geoObject->getX2Minimum(),geoObject->getX3Minimum(),geoObject->getX1Maximum(),geoObject->getX2Maximum(),geoObject->getX3Maximum() };
   //grid->getBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*ext, AABB[1]-(double)val<2>(blockNX)*ext, AABB[2]-(double)val<3>(blockNX)*ext, AABB[3]+(double)val<1>(blockNX)*ext, AABB[4]+(double)val<2>(blockNX)*ext, AABB[5]+(double)val<3>(blockNX)*ext, blockVector);
   //if (blockVector.size() == 0)
   //{
   //   auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   //   if (geometry->getId() == peGeo->getID())
   //   {
   //      geometry->setInactive();
   //      geometry->setGeometry(NULL);
   //      return;
   //   }
   //   else
   //      throw UbException(UB_EXARGS, "PeGeo ID is not matching!");
   //}


   //if (peGeo->getID() < 0 || peGeo->getID() >= physicsEngineGeometrieAdapters.size()) return;

   auto geometry = std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[peGeo->getID()]);
   geometry->shadowCounter--;

   if (geometry->getSystemId() == peGeo->getSystemID())
   {
      if (geometry->shadowCounter == 0)
      {
         geometry->setInactive();
         geometry->setGeometry(NULL);
      }
      else
      {
         geometry->setActive();
         geometry->setGeometry(peGeo);
      }
      //if (comm->getProcessID()==3)UBLOG(logINFO, "DemCoProcessor::removePeShadowGeo() geoID = "<<peGeo->getID()<<" counter="<<geometry->counter<<" shadowCounter="<<geometry->shadowCounter);
      return;
   }
   else
      throw UbException(UB_EXARGS, "PeGeo ID="+UbSystem::toString(peGeo->getID())+" is not matching geometry ID="+UbSystem::toString(geometry->getId()));
}

bool DemCoProcessor::isSpheresIntersection(double centerX1, double centerX2, double centerX3, double d)
{
   bool result = false;
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         SPtr<GbObject3D> sphere = interactors[i]->getGbObject3D();
         result = result || (sqrt(pow(sphere->getX1Centroid()-centerX1, 2)+pow(sphere->getX2Centroid()-centerX2, 2)+pow(sphere->getX3Centroid()-centerX3, 2)) < d);
      }
   }
   std::vector<int> values;
   values.push_back((int)result);
   std::vector<int> rvalues = comm->gather(values);

   if (comm->isRoot())
   {
      for (int i = 0; i < (int)rvalues.size(); i++)
      {
         result = result || (bool)rvalues[i];
      }
   }
   int iresult = (int)result;
   comm->broadcast(iresult);
   result = (bool)iresult;

   return result;
}

void DemCoProcessor::distributeIDs()
{
   std::vector<unsigned long long> peIDsSend;
   std::vector<int> vfIDsSend;

   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
      {
         peIDsSend.push_back(std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->getSystemId());
         vfIDsSend.push_back(interactors[i]->getID());
      }
   }

   std::vector<unsigned long long> peIDsRecv;
   std::vector<int> vfIDsRecv;

   comm->allGather(peIDsSend, peIDsRecv);
   comm->allGather(vfIDsSend, vfIDsRecv);

   std::map<int, unsigned long long> idMap;

   for (int i = 0; i < peIDsRecv.size(); i++)
   {
      idMap.insert(std::make_pair(vfIDsRecv[i], peIDsRecv[i]));
   }

   for (int i = 0; i < interactors.size(); i++)
   {
      std::map<int, unsigned long long>::const_iterator it;
      if ((it=idMap.find(interactors[i]->getID())) == idMap.end())
      {
         throw UbException(UB_EXARGS, "Interactor ID is invalid! The DEM object may be not in PE domain!");
      }


      std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->setSystemId(it->second);
      //std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->setId(idMap.find(interactors[i]->getID())->second);
   }

   //for (int i = 0; i < physicsEngineGeometrieAdapters.size(); i++)
   //{
   //   //physicsEngineSolver->updateGeometry(physicsEngineGeometries[i]);
   //   if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive())
   //   {
   //      interactors[i]->setPhysicsEngineGeometry(physicsEngineGeometrieAdapters[i]);
   //   }
   //}
}

void DemCoProcessor::setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor)
{
   this->boundaryConditionsBlockVisitor = boundaryConditionsBlockVisitor;
}
//////////////////////////////////////////////////////////////////////////
walberla::pe::RigidBody* DemCoProcessor::getPeGeoObject(walberla::id_t id)
{
   std::shared_ptr<walberla::blockforest::BlockForest> forest = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
   std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();

   //for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   //{

   //   for (auto bodyIt = walberla::pe::BodyIterator::begin(*blockIt, *storageId.get()); bodyIt != walberla::pe::BodyIterator::end(); ++bodyIt)
   //   {
   //      if (bodyIt->getID() == id)
   //      {
   //         walberla::pe::RigidBody* geo = *(bodyIt);
   //         return geo;
   //      }
   //   }
   //}

   std::shared_ptr<walberla::pe::BodyStorage> globalBodyStorage = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getGlobalBodyStorage();
   //walberla::pe::BodyID peGeo =  walberla::pe::getBody(*globalBodyStorage, *forest, *storageId, id);
   //if (peGeo != nullptr)
   //{
   //   return peGeo;
   //}
   //return NULL;


   return walberla::pe::getBody(*globalBodyStorage, *forest, *storageId, id);



   //walberla::pe::RigidBody* peBody = NULL;

   //std::shared_ptr<walberla::pe::BodyStorage> globalStorage = std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getGlobalBodyStorage();

   //auto bodyIt = globalStorage->find(id);
   //if (bodyIt != globalStorage->end())
   //{
   //   return *bodyIt;
   //}

   //for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   //{
   //   walberla::pe::Storage* storage = blockIt->getData< walberla::pe::Storage >(*storageId.get());
   //   walberla::pe::BodyStorage* bodyStorage = &(*storage)[0];
   //   walberla::pe::BodyStorage* bodystorageShadowCopies = &(*storage)[1];
   //   
   //   auto bodyIt = bodyStorage->find(id);
   //   if (bodyIt != bodyStorage->end())
   //   {
   //      return *bodyIt;
   //   }

   //   bodyIt = bodystorageShadowCopies->find(id);
   //   if (bodyIt != bodystorageShadowCopies->end())
   //   {
   //      return *bodyIt;
   //   }
   //}
   //return NULL;
}