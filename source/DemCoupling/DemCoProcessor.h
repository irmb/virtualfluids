/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef DEM_CO_PROCESSOR_H
#define DEM_CO_PROCESSOR_H

#include <memory>
#include <vector>
#include <map>

#include "Vector3D.h"

#include "CoProcessor.h"
#include "UbTuple.h"

class PhysicsEngineGeometryAdapter;
class PhysicsEngineSolverAdapter;
class PhysicsEngineMaterialAdapter;

class UbScheduler;
class Grid3D;
class ForceCalculator;
class Communicator;
class MovableObjectInteractor;
class Communicator;
class BoundaryConditionsBlockVisitor;


class DemCoProcessor : public CoProcessor
{
public:
    DemCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm, std::shared_ptr<ForceCalculator> forceCalculator, std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps = 1.0);
    virtual ~DemCoProcessor();

    void addInteractor(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial, Vector3D initalVelocity = Vector3D(0.0, 0.0, 0.0));
    void process(double step) override;
    std::shared_ptr<PhysicsEngineSolverAdapter> getPhysicsEngineSolver();
    void distributeIDs();
    void setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> blockVisitor);
    bool isGeoObjectInAABB(std::array<double,6> AABB);
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles);
  
private:
    std::shared_ptr<PhysicsEngineGeometryAdapter> createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const;
    void applyForcesOnGeometries();
    void setForcesToObject(SPtr<Grid3D> grid, std::shared_ptr<MovableObjectInteractor> interactor, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry);
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
    SPtr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor;
};


#endif

