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

#include "PhysicsEngineMaterialAdapter.h"
#include "PhysicsEngineGeometryAdapter.h"
#include "PhysicsEngineSolverAdapter.h"
#include "PePhysicsEngineGeometryAdapter.h"

#include "BoundaryConditions.h"
#include "Block3D.h"
#include "BCArray3D.h"
#include "MPICommunicator.h"
#include "BoundaryConditionsBlockVisitor.h"


#include <array>

DemCoProcessor::DemCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Communicator> comm, std::shared_ptr<ForceCalculator> forceCalculator, std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps) :
   CoProcessor(grid, s), comm(comm), forceCalculator(forceCalculator), physicsEngineSolver(physicsEngineSolver), intermediateDemSteps(intermediatePeSteps)
{

}

DemCoProcessor::~DemCoProcessor()
{

}

void DemCoProcessor::addInteractor(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial, Vector3D initalVelocity)
{
   interactors.push_back(interactor);
   const int id = static_cast<int>(interactors.size()) - 1;
   interactor->setID(id);
   const auto peGeometry = this->createPhysicsEngineGeometryAdapter(interactor, physicsEngineMaterial);
   if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometry)->isActive())
   {
      peGeometry->setLinearVelolocity(initalVelocity);
      physicsEngineGeometries.push_back(peGeometry);
   }
   else
   {
      physicsEngineGeometries.push_back(peGeometry);
   }
   //distributeIDs();
}


std::shared_ptr<PhysicsEngineGeometryAdapter> DemCoProcessor::createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const
{
   const int id = static_cast<int>(interactors.size()) - 1;
   SPtr<GbSphere3D> vfSphere = std::static_pointer_cast<GbSphere3D>(interactor->getGbObject3D());
   const Vector3D position(vfSphere->getX1Centroid(), vfSphere->getX2Centroid(), vfSphere->getX3Centroid());

   auto peGeometry = this->physicsEngineSolver->createPhysicsEngineGeometryAdapter(id, position, vfSphere->getRadius(), physicsEngineMaterial);
   if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometry)->isActive())
   {
      interactor->setPhysicsEngineGeometry(peGeometry);
      return peGeometry;
   }
   else
   {
      return peGeometry;
   }
}


void DemCoProcessor::process(double actualTimeStep)
{
   this->applyForcesOnGeometries();

   if (scheduler->isDue(actualTimeStep))
   {
      //UBLOG(logINFO, "DemCoProcessor::update - START - timestep = " << actualTimeStep);
      const double demTimeStepsPerIteration = scheduler->getMinStep();

      if (demTimeStepsPerIteration != 1)
         this->scaleForcesAndTorques(1.0 / demTimeStepsPerIteration);

      if (this->intermediateDemSteps == 1)
         this->calculateDemTimeStep(demTimeStepsPerIteration);
      
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

      this->moveVfGeoObject();

      grid->accept(*boundaryConditionsBlockVisitor.get());

      //UBLOG(logINFO, "DemCoProcessor::update - END - timestep = " << actualTimeStep);
   }
}
//////////////////////////////////////////////////////////////////////////
std::shared_ptr<PhysicsEngineSolverAdapter> DemCoProcessor::getPhysicsEngineSolver()
{
   return physicsEngineSolver;
}

void DemCoProcessor::applyForcesOnGeometries()
{
   for (int i = 0; i < physicsEngineGeometries.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         this->setForcesToObject(grid, interactors[i], physicsEngineGeometries[i]);

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
   for (int i = 0; i < physicsEngineGeometries.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         const Vector3D force = physicsEngineGeometries[i]->getForce() * scalingFactor;
         const Vector3D torque = physicsEngineGeometries[i]->getTorque() * scalingFactor;

         physicsEngineGeometries[i]->resetForceAndTorque();

         physicsEngineGeometries[i]->setForce(force);
         physicsEngineGeometries[i]->setTorque(torque);

         //UBLOG(logINFO, "F: (x,y,z) " << force);
         //UBLOG(logINFO, "T: (x,y,z) " << torque);
      }
   }
}


void DemCoProcessor::calculateDemTimeStep(double step) const
{
   physicsEngineSolver->runTimestep(step);

   for (int i = 0; i < physicsEngineGeometries.size(); i++)
   {
      physicsEngineSolver->updateGeometry(physicsEngineGeometries[i]);
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         interactors[i]->setPhysicsEngineGeometry(physicsEngineGeometries[i]);
      }
   }
}

void DemCoProcessor::moveVfGeoObject()
{
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         interactors[i]->moveGbObjectTo(physicsEngineGeometries[i]->getPosition());
      }
   }
}

bool  DemCoProcessor::isDemObjectInAABB(std::array<double, 6> AABB)
{
   bool result = false;
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
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
      for (int i = 0; i < (int)rvalues.size(); i ++)
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
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         interactors[i]->getGbObject3D()->addSurfaceTriangleSet(nodes, triangles);
      }
   }
}

void DemCoProcessor::getObjectsPropertiesVector(std::vector<double>& p)
{
   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         p.push_back(interactors[i]->getGbObject3D()->getX1Centroid());
         p.push_back(interactors[i]->getGbObject3D()->getX2Centroid());
         p.push_back(interactors[i]->getGbObject3D()->getX3Centroid());
         Vector3D v = physicsEngineGeometries[i]->getLinearVelocity();
         p.push_back(v[0]);
         p.push_back(v[1]);
         p.push_back(v[2]);
      }
   }
}

void DemCoProcessor::distributeIDs()
{
   std::vector<int> peIDsSend;
   std::vector<int> vfIDsSend;

   for (int i = 0; i < interactors.size(); i++)
   {
      if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
      {
         peIDsSend.push_back(std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->getId());
         vfIDsSend.push_back(interactors[i]->getID());
      }
   }

   std::vector<int> peIDsRecv;
   std::vector<int> vfIDsRecv;

   comm->allGather(peIDsSend, peIDsRecv);
   comm->allGather(vfIDsSend, vfIDsRecv);

   std::map<int, int> idMap;

   for (int i = 0; i < peIDsRecv.size(); i++)
   {
      idMap.insert(std::make_pair(vfIDsRecv[i], peIDsRecv[i]));
   }

   for (int i = 0; i < interactors.size(); i++)
   {
       std::map<int, int>::const_iterator it;
      if ((it=idMap.find(interactors[i]->getID())) == idMap.end())
      {
         throw UbException(UB_EXARGS, "not valid ID!");
      }
      

      std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->setId(it->second);
      //std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->setId(idMap.find(interactors[i]->getID())->second);
   }

   for (int i = 0; i < physicsEngineGeometries.size(); i++)
   {
         physicsEngineSolver->updateGeometry(physicsEngineGeometries[i]);
         if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[i])->isActive())
         {
            interactors[i]->setPhysicsEngineGeometry(physicsEngineGeometries[i]);
         }
   }
}

void DemCoProcessor::setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor)
{
   this->boundaryConditionsBlockVisitor = boundaryConditionsBlockVisitor;
}
