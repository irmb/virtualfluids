/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef PHYSICS_ENGINE_SOLVER_ADAPTER_H
#define PHYSICS_ENGINE_SOLVER_ADAPTER_H

#include "Vector3D.h"

class PhysicsEngineGeometryAdapter;
class PhysicsEngineMaterialAdapter;

class PhysicsEngineSolverAdapter
{
public:
    virtual ~PhysicsEngineSolverAdapter() {}

    virtual std::shared_ptr<PhysicsEngineGeometryAdapter>
    createPhysicsEngineGeometryAdapter(int id, const Vector3D &position, double radius,
                                       std::shared_ptr<PhysicsEngineMaterialAdapter> material) const = 0;
    virtual void runTimestep(double step)                                                            = 0;
};

#endif
