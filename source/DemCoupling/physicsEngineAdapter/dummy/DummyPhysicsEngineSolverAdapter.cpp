#include "DummyPhysicsEngineSolverAdapter.h"

#include "DummyPhysicsEngineGeometryAdapter.h"


std::shared_ptr<PhysicsEngineGeometryAdapter> DummyPhysicsEngineSolverAdapter::createPhysicsEngineGeometryAdapter(int id, const Vector3D& position, double radius, std::shared_ptr<PhysicsEngineMaterialAdapter> material) const
{
    return std::static_pointer_cast<PhysicsEngineGeometryAdapter>(std::make_shared<DummyPhysicsEngineGeometryAdapter>());
}

void DummyPhysicsEngineSolverAdapter::runTimestep(double step)
{

}
