/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef D3Q27_MOVABLE_OBJECT_INTERACTOR_H
#define D3Q27_MOVABLE_OBJECT_INTERACTOR_H

#include <memory>
#include <vector>

#include <boost/serialization/base_object.hpp>

#include <VirtualFluids/Interactors/D3Q27Interactor.h>

#include <VirtualFluidsBasics/basics/utilities/Vector3D.h>
#include "physicsEngineAdapter/PhysicsEngineGeometryAdapter.h"

class MovableObjectInteractor;
typedef std::shared_ptr<MovableObjectInteractor> MovableObjectInteractorPtr;

class Grid3D;
class Block3D;
class BCArray3D;
class BCAdapter;
class GbObject3D;
class BoundaryConditionsBlockVisitor;

class PhysicsEngineGeometryAdapter;
class Reconstructor;

class MovableObjectInteractor : public D3Q27Interactor
{
public:
    MovableObjectInteractor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type, std::shared_ptr<Reconstructor> reconstructor, State isPinned);
    virtual ~MovableObjectInteractor();

    void setPhysicsEngineGeometry(std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry);
    void setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> blockVisitor);

    void moveGbObjectTo(const Vector3D& position);

private:
    void rearrangeGrid();
    void setSolidNodesToFluid();
    void reconstructSolidNodes();
    void setBcs() const;
    void setBcBlocks();

    void updateVelocityBc();
    void setGeometryVelocityToBoundaryCondition(std::vector<int> node, std::shared_ptr<Block3D> block, std::shared_ptr<BCArray3D> bcArray) const;

    std::shared_ptr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor;
    std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry;

    std::shared_ptr<Reconstructor> reconstructor;
    State state;

protected:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<D3Q27Interactor>(*this);
    }
};


#endif
