#include "MovableObjectInteractor.h"

#include <basics/utilities/UbLogger.h>
#include <numerics/geometry3d/GbObject3D.h>
#include "VirtualFluidsBasics/basics/utilities/Vector3D.h"

#include "Block3D.h"
#include "Grid3D.h"
#include "BCArray3D.h"
#include "BCAdapter.h"
#include "BCProcessor.h"
#include "ILBMKernel.h"

#include "SetSolidBlockVisitor.h"
#include "BoundaryConditionsBlockVisitor.h"

#include "physicsEngineAdapter/PhysicsEngineGeometryAdapter.h"
#include "reconstructor/Reconstructor.h"


MovableObjectInteractor::MovableObjectInteractor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type, std::shared_ptr<Reconstructor> reconstructor, State state)
   : D3Q27Interactor(geoObject3D, grid, bcAdapter, type), reconstructor(reconstructor), state(state)
{

}

MovableObjectInteractor::~MovableObjectInteractor()
{

}

void MovableObjectInteractor::setPhysicsEngineGeometry(std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry)
{
    this->physicsEngineGeometry = physicsEngineGeometry;
    physicsEngineGeometry->changeState(this->state);
}

void MovableObjectInteractor::setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor)
{
    this->boundaryConditionsBlockVisitor = boundaryConditionsBlockVisitor;
}

void MovableObjectInteractor::moveGbObjectTo(const Vector3D& position)
{
    //UBLOG(logINFO, "new position: (x,y,z) " << val<1>(position) << ", " << val<2>(position) << ", " << val<3>(position));

    this->getGbObject3D()->setCenterCoordinates(UbTupleDouble3(position[0], position[1], position[2]));
    this->rearrangeGrid();
}

void MovableObjectInteractor::rearrangeGrid()
{
    this->reconstructSolidNodes();
    this->setSolidNodesToFluid();

    this->removeSolidBlocks();
    this->removeBcBlocks();

    this->setBcBlocks();

    this->initInteractor();

    this->setBcs();
    this->updateVelocityBc();
}

void MovableObjectInteractor::reconstructSolidNodes()
{
    for(SolidNodeIndicesMap::value_type t : solidNodeIndicesMap)
    {
        Block3DPtr block = t.first;
        std::set< UbTupleInt3 >& solidNodeIndices = t.second;

        
        ILBMKernelPtr kernel = block->getKernel();

        for (UbTupleInt3 node : solidNodeIndices)
        {
            const int x1 = val<1>(node);
            const int x2 = val<2>(node);
            const int x3 = val<3>(node);

            const Vector3D worldCoordinates = this->grid.lock()->getNodeCoordinates(block, x1, x2, x3);
    
            if (kernel->isInsideOfDomain(x1, x2, x3))
                reconstructor->reconstructNode(x1, x2, x3, worldCoordinates, physicsEngineGeometry, kernel);
        }
    }
}

void MovableObjectInteractor::setSolidNodesToFluid()
{
    for (SolidNodeIndicesMap::value_type t : solidNodeIndicesMap)
    {
        Block3DPtr block = t.first;
        std::set< UbTupleInt3 >& solidNodeIndices = t.second;

        ILBMKernelPtr kernel = block->getKernel();
        BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();

        for (UbTupleInt3 node : solidNodeIndices)
            bcArray->setFluid(val<1>(node), val<2>(node), val<3>(node));
    }
}

void MovableObjectInteractor::setBcBlocks()
{
    SetSolidBlockVisitor v(shared_from_this(), BlockType::BC);
    this->grid.lock()->accept(v);
}

void MovableObjectInteractor::setBcs() const
{
    this->grid.lock()->accept(*boundaryConditionsBlockVisitor.get());
}

void MovableObjectInteractor::updateVelocityBc()
{
    for(BcNodeIndicesMap::value_type t : this->getBcNodeIndicesMap())
    {
        Block3DPtr block = t.first;
        std::set< std::vector<int> >& bcNodeIndices = t.second;

        BCArray3DPtr bcArray = block->getKernel()->getBCProcessor()->getBCArray();

        for(std::vector<int> node : bcNodeIndices)
            setGeometryVelocityToBoundaryCondition(node, block, bcArray);
    }
}


void MovableObjectInteractor::setGeometryVelocityToBoundaryCondition(std::vector<int> node, Block3DPtr block, BCArray3DPtr bcArray) const
{
    const BoundaryConditionsPtr bc = bcArray->getBC(node[0], node[1], node[2]);
    const Vector3D worldCoordinates = this->grid.lock()->getNodeCoordinates(block, node[0], node[1], node[2]);
    const Vector3D velocity = this->physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);
    bc->setBoundaryVelocity(velocity);
}
