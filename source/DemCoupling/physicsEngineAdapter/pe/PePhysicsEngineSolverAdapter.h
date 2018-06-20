/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef PE_PHYSICS_ENGINE_SOLVER_ADAPTER_H
#define PE_PHYSICS_ENGINE_SOLVER_ADAPTER_H

#include <memory>
#include <shared_mutex>

#include <pe/basic.h>
#include "UbTuple.h"

#include "PhysicsEngineSolverAdapter.h"
#include "PePhysicsEngineSolverAdapter.h"


class PePhysicsEngineMaterialAdapter;
class PhysicsEngineGeometryAdapter;

namespace walberla
{
    namespace domain_decomposition
    {
        class BlockDataID;
    }
    namespace blockforest
    {
        class BlockForest;
    }
    namespace pe
    {
        class BodyStorage;
        class RigidBody;
        namespace cr
        {
            class HardContactSemiImplicitTimesteppingSolvers;
        }
    }
}

struct PeParameter
{
   PeParameter(double relaxationParameter, int maxPeIterations, Vector3D globalLinearAcceleration, std::shared_ptr<PePhysicsEngineMaterialAdapter> planes, std::array<double, 6> simulationDomain, UbTupleInt3 numberOfBlocks, UbTupleBool3 isPeriodic, Vector3D minOffset, Vector3D maxOffset)
        : relaxationParameter(relaxationParameter), maxPeIterations(maxPeIterations), globalLinearAcceleration(globalLinearAcceleration), simulationDomain(simulationDomain), numberOfBlocks(numberOfBlocks), isPeriodic(isPeriodic), planes(planes), minOffset(minOffset), maxOffset(maxOffset)
    {
    }

    double relaxationParameter;
    int maxPeIterations;
    Vector3D globalLinearAcceleration;

    std::array<double, 6> simulationDomain;
    UbTupleInt3 numberOfBlocks;
    UbTupleBool3 isPeriodic;

    std::shared_ptr<PePhysicsEngineMaterialAdapter> planes;

    Vector3D minOffset;
    Vector3D maxOffset;

};

class PePhysicsEngineSolverAdapter : public PhysicsEngineSolverAdapter
{
public:
    PePhysicsEngineSolverAdapter(std::shared_ptr<PeParameter> peParameter);
    virtual ~PePhysicsEngineSolverAdapter() {}

    std::shared_ptr<PhysicsEngineGeometryAdapter> createPhysicsEngineGeometryAdapter(int id, const Vector3D& position, double radius, std::shared_ptr<PhysicsEngineMaterialAdapter> material) const override;
    void runTimestep(double step) override;
    walberla::pe::RigidBody* getPeGeoObject(walberla::id_t id);
    void updateGeometry(std::shared_ptr<PhysicsEngineGeometryAdapter>) override;
    std::shared_ptr< walberla::blockforest::BlockForest > getForest();

private:
    void initalizePeEnvironment();
    void initialPeBodyStorage();
    void initalPeBlockForest();
    void initalBlockData();

    void initalPeIntegrator();
    static void executePeBodyTypeTuple();
    void initalPeChannel() const;

private:
    std::shared_ptr<PeParameter> peParameter;

    std::shared_ptr<walberla::pe::BodyStorage> globalBodyStorage;
    std::shared_ptr< walberla::blockforest::BlockForest > forest;
    std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId;
    std::shared_ptr<walberla::pe::cr::HardContactSemiImplicitTimesteppingSolvers> cr;
};

#endif

