/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DUMMY_PHYSICS_ENGINE_SOLVER_ADAPTER_H
#define DUMMY_PHYSICS_ENGINE_SOLVER_ADAPTER_H

#include <memory>

#include "UbTuple.h"

#include "PhysicsEngineSolverAdapter.h"



class DummyPhysicsEngineSolverAdapter : public PhysicsEngineSolverAdapter
{
public:
    DummyPhysicsEngineSolverAdapter() {};
    virtual ~DummyPhysicsEngineSolverAdapter() {}

    std::shared_ptr<PhysicsEngineGeometryAdapter> createPhysicsEngineGeometryAdapter(int id, const Vector3D& position, double radius, std::shared_ptr<PhysicsEngineMaterialAdapter> material) const override;
    void runTimestep(double step) override;

};

#endif

