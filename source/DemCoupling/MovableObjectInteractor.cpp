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
#include "CoordinateTransformation3D.h"

#include "SetBcBlocksBlockVisitor.h"
#include "BoundaryConditionsBlockVisitor.h"

#include "PhysicsEngineGeometryAdapter.h"
#include "Reconstructor.h"

#include <array>

//#define TIMING

#ifdef TIMING
   #include "UbTiming.h"
#endif


MovableObjectInteractor::MovableObjectInteractor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type, std::shared_ptr<Reconstructor> reconstructor, State state)
   : D3Q27Interactor(geoObject3D, grid, bcAdapter, type), reconstructor(reconstructor), state(state)
{
   //grid->getBlocks(0, grid->getRank(), true, blockVector);
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
#ifdef TIMING
   UbTimer timer;
   timer.resetAndStart();
#endif

#ifdef TIMING
   UBLOG(logINFO, "MovableObjectInteractor::rearrangeGrid():start");
#endif

    this->reconstructDistributionOnSolidNodes();

#ifdef TIMING
    UBLOG(logINFO, "reconstructDistributionOnSolidNodes() time = "<<timer.stop()<<" s");
#endif

    this->setSolidNodesToFluid();

#ifdef TIMING
    UBLOG(logINFO, "setSolidNodesToFluid() time = "<<timer.stop()<<" s");
#endif

    this->setBcNodesToFluid();

#ifdef TIMING
    UBLOG(logINFO, "setBcNodesToFluid() time = "<<timer.stop()<<" s");
#endif

    this->removeSolidBlocks();

#ifdef TIMING
    UBLOG(logINFO, "removeSolidBlocks() time = "<<timer.stop()<<" s");
#endif

    this->removeBcBlocks();

#ifdef TIMING
    UBLOG(logINFO, "removeBcBlocks() time = "<<timer.stop()<<" s");
#endif

    this->setBcBlocks();

#ifdef TIMING
    UBLOG(logINFO, "setBcBlocks() time = "<<timer.stop()<<" s");
#endif

    this->initInteractor();

#ifdef TIMING
    UBLOG(logINFO, "initInteractor() time = "<<timer.stop()<<" s");
#endif

    this->updateVelocityBc();

#ifdef TIMING
    UBLOG(logINFO, "updateVelocityBc() time = "<<timer.stop()<<" s");
#endif
}

void MovableObjectInteractor::updateNodeLists()
{
   //for (BcNodeIndicesMap::value_type t : bcNodeIndicesMap)
   //{
   //   SPtr<Block3D> block = t.first;
   //   std::set< UbTupleInt3 >& bcNodeIndices = t.second;


   //   SPtr<ILBMKernel> kernel = block->getKernel();

   //   for (UbTupleInt3 node : bcNodeIndices)
   //   {

   //   }
   //}
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

    //////////////////////////////////////////////////////////////////////////
   //SPtr<GbObject3D> geoObject = this->getGbObject3D();
   //std::array<double, 6> AABB ={ geoObject->getX1Minimum(),geoObject->getX2Minimum(),geoObject->getX3Minimum(),geoObject->getX1Maximum(),geoObject->getX2Maximum(),geoObject->getX3Maximum() };
   //blockVector.clear();
   //UbTupleInt3 blockNX=grid.lock()->getBlockNX();
   //double ext = 0.0;
   //grid.lock()->getBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*ext, AABB[1]-(double)val<2>(blockNX)*ext, AABB[2]-(double)val<3>(blockNX)*ext, AABB[3]+(double)val<1>(blockNX)*ext, AABB[4]+(double)val<2>(blockNX)*ext, AABB[5]+(double)val<3>(blockNX)*ext, blockVector);

   //for(std::shared_ptr<Block3D> block : this->blockVector)
   //{
   //   if (block->getKernel())
   //   {
   //      setBCBlock(block);
   //   }
   //}
   //////////////////////////////////////////////////////////////////////////
   //SPtr<GbObject3D> geoObject = this->getGbObject3D();
   //std::array <double, 2> minMax1;
   //std::array <double, 2> minMax2;
   //std::array <double, 2> minMax3;
   //minMax1[0] = geoObject->getX1Minimum();
   //minMax2[0] = geoObject->getX2Minimum();
   //minMax3[0] = geoObject->getX3Minimum();
   //minMax1[1] = geoObject->getX1Maximum();
   //minMax2[1] = geoObject->getX2Maximum();
   //minMax3[1] = geoObject->getX3Maximum();

   //SPtr<CoordinateTransformation3D> trafo = grid.lock()->getCoordinateTransformator();

   //for (int x3 = 0; x3 < 2; x3++)
   //   for (int x2 = 0; x2 < 2; x2++)
   //      for (int x1 = 0; x1 < 2; x1++)
   //      {
   //         int ix1 = (int)trafo->transformForwardToX1Coordinate(minMax1[x1], minMax2[x2], minMax3[x3]);
   //         int ix2 = (int)trafo->transformForwardToX2Coordinate(minMax1[x1], minMax2[x2], minMax3[x3]);
   //         int ix3 = (int)trafo->transformForwardToX3Coordinate(minMax1[x1], minMax2[x2], minMax3[x3]);
   //         blockVector.push_back(grid.lock()->getBlock(ix1, ix2, ix3, 0));
   //      }
   //for(std::shared_ptr<Block3D> block : this->blockVector)
   //{
   //   if (block->getKernel())
   //   {
   //      setBCBlock(block);
   //   }
   //}
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
