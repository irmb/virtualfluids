#include "MovableObjectInteractor.h"

#include "UbLogger.h"
#include "GbObject3D.h"
#include "Vector3D.h"

#include "Block3D.h"
#include "Grid3D.h"
#include "BCArray3D.h"
#include "BCAdapter.h"
#include "BCProcessor.h"
#include "ILBMKernel.h"

#include "SetBcBlocksBlockVisitor.h"
#include "BoundaryConditionsBlockVisitor.h"

#include "PhysicsEngineGeometryAdapter.h"
#include "Reconstructor.h"


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

void MovableObjectInteractor::moveGbObjectTo(const Vector3D& position)
{
    //UBLOG(logINFO, "new position: (x,y,z) " << val<1>(position) << ", " << val<2>(position) << ", " << val<3>(position));

    this->getGbObject3D()->setCenterCoordinates(UbTupleDouble3(position[0], position[1], position[2]));
    this->rearrangeGrid();
}

void MovableObjectInteractor::rearrangeGrid()
{
    this->reconstructDistributionOnSolidNodes();
    
    this->setSolidNodesToFluid();
    this->setBcNodesToFluid();

    this->removeSolidBlocks();
    this->removeBcBlocks();

    this->setBcBlocks();

    this->initInteractor();

    this->updateVelocityBc();
}

void MovableObjectInteractor::reconstructDistributionOnSolidNodes()
{
    for(SolidNodeIndicesMap::value_type t : solidNodeIndicesMap)
    {
        SPtr<Block3D> block = t.first;
        std::set< UbTupleInt3 >& solidNodeIndices = t.second;

        
        SPtr<ILBMKernel> kernel = block->getKernel();

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
        SPtr<Block3D> block = t.first;
        std::set< UbTupleInt3 >& solidNodeIndices = t.second;

        SPtr<ILBMKernel> kernel = block->getKernel();
        SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

        for (UbTupleInt3 node : solidNodeIndices)
            bcArray->setFluid(val<1>(node), val<2>(node), val<3>(node));
    }
}

void MovableObjectInteractor::setBcNodesToFluid()
{
   for (BcNodeIndicesMap::value_type t : bcNodeIndicesMap)
   {
      SPtr<Block3D> block = t.first;
      std::set< std::vector<int> >& bcNodeIndices = t.second;

      SPtr<ILBMKernel> kernel = block->getKernel();
      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

      for (std::vector<int> node : bcNodeIndices)
         bcArray->setFluid(node[0], node[1], node[2]);
   }
}

void MovableObjectInteractor::setBcBlocks()
{
    SetBcBlocksBlockVisitor v(shared_from_this());
    this->grid.lock()->accept(v);
}

void MovableObjectInteractor::updateVelocityBc()
{
    for(BcNodeIndicesMap::value_type t : this->getBcNodeIndicesMap())
    {
        SPtr<Block3D> block = t.first;
        std::set< std::vector<int> >& bcNodeIndices = t.second;

        SPtr<BCArray3D> bcArray = block->getKernel()->getBCProcessor()->getBCArray();

        for(std::vector<int> node : bcNodeIndices)
            setGeometryVelocityToBoundaryCondition(node, block, bcArray);
    }
}


void MovableObjectInteractor::setGeometryVelocityToBoundaryCondition(std::vector<int> node, SPtr<Block3D> block, SPtr<BCArray3D> bcArray) const
{
    const SPtr<BoundaryConditions> bc = bcArray->getBC(node[0], node[1], node[2]);
    const Vector3D worldCoordinates = this->grid.lock()->getNodeCoordinates(block, node[0], node[1], node[2]);
    const Vector3D velocity = this->physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);

    bc->setBoundaryVelocity(velocity);
}
