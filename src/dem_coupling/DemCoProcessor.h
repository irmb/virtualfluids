/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef DEM_CO_PROCESSOR_H
#define DEM_CO_PROCESSOR_H

#include <memory>
#include <vector>

#include "VirtualFluidsBasics/basics/utilities/Vector3D.h"

#include <boost/serialization/access.hpp>

#include <VirtualFluids/CoProcessors/CoProcessor.h>

class PhysicsEngineGeometryAdapter;
class PhysicsEngineSolverAdapter;
class PhysicsEngineMaterialAdapter;

class UbScheduler;
class Grid3D;
class ForceCalculator;
class Communicator;
class MovableObjectInteractor;


class DemCoProcessor : public CoProcessor
{
public:
    DemCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm, std::shared_ptr<ForceCalculator> forceCalculator, std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps = 1.0);
    virtual ~DemCoProcessor();

    void addInteractor(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial, Vector3D initalVelocity = Vector3D(0.0, 0.0, 0.0));

    void process(double step) override;
  
private:
    std::shared_ptr<PhysicsEngineGeometryAdapter> createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const;

    void applyForcesOnObjects();
    void setForcesToObject(Grid3DPtr grid, std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry);
    void scaleForcesAndTorques(double scalingFactor);

    void calculateDemTimeStep(double step) const;
    void moveVfGeoObject();

private:
    std::shared_ptr<Communicator> comm;
    std::vector<std::shared_ptr<MovableObjectInteractor> > interactors;
    std::shared_ptr<ForceCalculator> forceCalculator;

    std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver;
    std::vector<std::shared_ptr<PhysicsEngineGeometryAdapter> > physicsEngineGeometries;

    double intermediateDemSteps;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        //ar & boost::serialization::base_object<CoProcessor>(*this);
    }
};


#endif

