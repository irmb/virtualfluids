/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef D3Q27_MOVABLE_OBJECT_INTERACTOR_H
#define D3Q27_MOVABLE_OBJECT_INTERACTOR_H

#include <memory>
#include <vector>

#include "D3Q27Interactor.h"

#include "Vector3D.h"
#include "PhysicsEngineGeometryAdapter.h"


class Grid3D;
class Block3D;
class BCArray3D;
class BCAdapter;
class GbObject3D;

class PhysicsEngineGeometryAdapter;
class Reconstructor;

class MovableObjectInteractor : public D3Q27Interactor
{
public:
   typedef std::map<SPtr<Block3D>, std::set< std::array<int,3> > > InBcNodeIndicesMap;
   typedef std::map<SPtr<Block3D>, std::set< std::array<int,3> > > OutBcNodeIndicesMap;
public:
    MovableObjectInteractor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type, std::shared_ptr<Reconstructor> reconstructor, State isPinned);
    virtual ~MovableObjectInteractor();

    void setPhysicsEngineGeometry(std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry);

    void moveGbObjectTo(const Vector3D& position);

private:
    void rearrangeGrid();
    void updateNodeLists();
    void setSolidNodesToFluid();
    void setBcNodesToFluid();
    void reconstructDistributionOnSolidNodes();
    void setBcBlocks();

    void updateVelocityBc();
    void setGeometryVelocityToBoundaryCondition(std::vector<int> node, std::shared_ptr<Block3D> block, std::shared_ptr<BCArray3D> bcArray) const;

    std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry;

    std::shared_ptr<Reconstructor> reconstructor;
    State state;
    std::vector< std::shared_ptr<Block3D> > blockVector;
};


#endif
