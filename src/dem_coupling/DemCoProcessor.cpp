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

#include <dem_coupling/physicsEngineAdapter/PhysicsEngineMaterialAdapter.h>
#include <dem_coupling/physicsEngineAdapter/PhysicsEngineGeometryAdapter.h>
#include <dem_coupling/physicsEngineAdapter/PhysicsEngineSolverAdapter.h>

#include "BoundaryConditions.h"
#include "Block3D.h"
#include "BCArray3D.h"

DemCoProcessor::DemCoProcessor( Grid3DPtr grid, UbSchedulerPtr s, CommunicatorPtr comm, std::shared_ptr<ForceCalculator> forceCalculator, std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps) :
    CoProcessor(grid, s), comm(comm), forceCalculator(forceCalculator), physicsEngineSolver(physicsEngineSolver), intermediateDemSteps(intermediatePeSteps)
{

}

DemCoProcessor::~DemCoProcessor()
{

}

void DemCoProcessor::addInteractor(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial, Vector3D initalVelocity)
{
    interactors.push_back(interactor);
    const auto peGeometry = this->createPhysicsEngineGeometryAdapter(interactor, physicsEngineMaterial);
    peGeometry->setLinearVelolocity(initalVelocity);
    physicsEngineGeometries.push_back(peGeometry);
}


std::shared_ptr<PhysicsEngineGeometryAdapter> DemCoProcessor::createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const
{
    const int id = static_cast<int>(interactors.size()) - 1;
    GbSphere3DPtr vfSphere = std::static_pointer_cast<GbSphere3D>(interactor->getGbObject3D());
    const Vector3D position(vfSphere->getX1Centroid(), vfSphere->getX2Centroid(), vfSphere->getX3Centroid());

    auto peGeometry = this->physicsEngineSolver->createPhysicsEngineGeometryAdapter(id, position, vfSphere->getRadius(), physicsEngineMaterial);
    interactor->setPhysicsEngineGeometry(peGeometry);
    return peGeometry;
}


void DemCoProcessor::process(double actualTimeStep)
{
    this->applyForcesOnObjects();

   if(scheduler->isDue(actualTimeStep) )
   {
       //UBLOG(logINFO, "DemCoProcessor::update - START" << step);
       const double timeStep = scheduler->getMinStep();

       if (timeStep != 1)
           this->scaleForcesAndTorques(1.0 / timeStep);

       if (this->intermediateDemSteps == 1)
           this->calculateDemTimeStep(timeStep);
           
       //const double subPeTimestep = timeStep / intermediateDemSteps;
       //
       //for (int i = 0; i < intermediateDemSteps; i++)
       //{
       //    // in the first set forces on local bodies are already set by force synchronization
       //    if (i != 0)
       //    {
       //        physicsEngineGeometries[i]->addForce();
       //        physicsEngineGeometries[i]->addTorque();
       //    }

       //    this->calculatePeTimeStep(subPeTimestep);
       //}

       this->moveVfGeoObject();
          
       //UBLOG(logINFO, "DemCoProcessor::update - END" << step);
   }
}

void DemCoProcessor::applyForcesOnObjects()
{
    for (int i = 0; i < physicsEngineGeometries.size(); i++)
    {
        this->setForcesToObject(grid, interactors[i], physicsEngineGeometries[i]);
        //physicsEngineGeometries[i]->setLinearVelolocity(Vector3D(0.05, 0.0, 0.0));
        //physicsEngineGeometries[i]->setAngularVelocity(Vector3D(0.01, 0.01, 0.01));
        UBLOG(logINFO, "F: (x,y,z) " << physicsEngineGeometries[i]->getLinearVelocity());

    }
}

void DemCoProcessor::setForcesToObject(Grid3DPtr grid, MovableObjectInteractorPtr interactor, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry)
{
    for (BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap())
    {
        Block3DPtr block = t.first;
        ILBMKernelPtr kernel = block->getKernel();
        BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();
        DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
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

                BoundaryConditionsPtr bc = bcArray->getBC(x1, x2, x3);
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
        const Vector3D force = physicsEngineGeometries[i]->getForce() * scalingFactor;
        const Vector3D torque = physicsEngineGeometries[i]->getTorque() * scalingFactor;

        physicsEngineGeometries[i]->resetForceAndTorque();

        physicsEngineGeometries[i]->setForce(force);
        physicsEngineGeometries[i]->setTorque(torque);

        //UBLOG(logINFO, "F: (x,y,z) " << force);
        //UBLOG(logINFO, "T: (x,y,z) " << torque);
    }
}


void DemCoProcessor::calculateDemTimeStep(double step) const
{
    physicsEngineSolver->runTimestep(step);

    for (int i = 0; i < physicsEngineGeometries.size(); i++)
        physicsEngineSolver->updateGeometry(physicsEngineGeometries[i]);
}

void DemCoProcessor::moveVfGeoObject()
{
    for (int i = 0; i < interactors.size(); i++)
        interactors[i]->moveGbObjectTo(physicsEngineGeometries[i]->getPosition());
}
