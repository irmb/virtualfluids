#include "GenBlocksGridVisitor.h"
#include "Grid3DSystem.h"
#include "CoordinateTransformation3D.h"
#include "Block3D.h"
#include "Grid3D.h"

#include <numerics/geometry3d/GbObject3D.h>

GenBlocksGridVisitor::GenBlocksGridVisitor(GbObject3DPtr boundingBox) :
   boundingBox(boundingBox)
{

}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::visit(const Grid3DPtr grid)
{
    double orgX1 = boundingBox->getX1Minimum();
    double orgX2 = boundingBox->getX2Minimum();
    double orgX3 = boundingBox->getX3Minimum();

    double dx = grid->getDeltaX(0);

    UbTupleInt3 blockNX = grid->getBlockNX();

    double blockLentghX1 = (double)val<1>(blockNX)*dx;
    double blockLentghX2 = (double)val<2>(blockNX)*dx;
    double blockLentghX3 = (double)val<3>(blockNX)*dx;

    CoordinateTransformation3DPtr trafo(new CoordinateTransformation3D(orgX1, orgX2, orgX3, blockLentghX1, blockLentghX2, blockLentghX3));
    grid->setCoordinateTransformator(trafo);

    genBlocks(grid);
}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::fillExtentWithBlocks( Grid3DPtr grid )
{
   for(int x3 =  val<3>(minInd); x3 <  val<3>(maxInd); x3++)
   {
      for(int x2 =  val<2>(minInd); x2 <  val<2>(maxInd); x2++)
      {
         for(int x1 =  val<1>(minInd); x1 <  val<1>(maxInd); x1++)
         {
            Block3DPtr block( new Block3D(x1,x2,x3,0) );
            grid->addBlock(block);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::genBlocks(Grid3DPtr grid)
{
    minInd = grid->getBlockIndexes(boundingBox->getX1Minimum(), boundingBox->getX2Minimum(), boundingBox->getX3Minimum());
    double geoMaxX1 = boundingBox->getX1Maximum();
    double geoMaxX2 = boundingBox->getX2Maximum();
    double geoMaxX3 = boundingBox->getX3Maximum();
    maxInd = grid->getBlockIndexes(geoMaxX1, geoMaxX2, geoMaxX3);
    UbTupleDouble3 blockCoord = grid->getBlockWorldCoordinates(static_cast<int>(val<1>(maxInd)), static_cast<int>(val<2>(maxInd)), static_cast<int>(val<3>(maxInd)), 0);
    if (geoMaxX1 > val<1>(blockCoord))
        val<1>(maxInd) += 1;
    if (geoMaxX2 > val<2>(blockCoord))
        val<2>(maxInd) += 1;
    if (geoMaxX3 > val<3>(blockCoord))
        val<3>(maxInd) += 1;

    this->fillExtentWithBlocks(grid);

    grid->setNX1(val<1>(maxInd));
    grid->setNX2(val<2>(maxInd));
    grid->setNX3(val<3>(maxInd));
}
