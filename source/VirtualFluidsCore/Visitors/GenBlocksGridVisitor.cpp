#include "GenBlocksGridVisitor.h"
#include "Grid3DSystem.h"
#include <boost/foreach.hpp>

GenBlocksGridVisitor::GenBlocksGridVisitor(GbObject3DPtr boundingBox) :
   boundingBox(boundingBox),
   nx1(0),
   nx2(0),
   nx3(0),
   withDeltaX(true)
{

}
GenBlocksGridVisitor::GenBlocksGridVisitor(int nx1, int nx2, int nx3) :
   nx1(nx1),
   nx2(nx2),
   nx3(nx3),
   withDeltaX(false)
{

}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::visit( const Grid3DPtr grid )
{
   findOrigin(grid);
   UbTupleInt3 blockNX = grid->getBlockNX();
   double blockLentghX1, blockLentghX2, blockLentghX3; 
   double dx;
   double geoMaxX1 = boundingBox->getX1Maximum();
   double geoMaxX2 = boundingBox->getX2Maximum();
   double geoMaxX3 = boundingBox->getX3Maximum();

   if (withDeltaX)
   {
      dx = grid->getDeltaX(0);
      blockLentghX1 = (double)val<1>(blockNX)*dx;
      blockLentghX2 = (double)val<2>(blockNX)*dx;
      blockLentghX3 = (double)val<3>(blockNX)*dx;
   } 
   else
   {
      int gNX1 = grid->getNX1();
      dx = boundingBox->getLengthX1()/double(val<1>(blockNX)*gNX1);
      grid->setDeltaX(dx);
      blockLentghX1 = val<1>(blockNX)*dx;
      blockLentghX2 = val<2>(blockNX)*dx;
      blockLentghX3 = val<3>(blockNX)*dx;
   }

   CoordinateTransformation3DPtr trafo(new CoordinateTransformation3D(orgX1,orgX2,orgX3,blockLentghX1,blockLentghX2,blockLentghX3));
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
   if(geoMaxX1 >  val<1>(blockCoord))
      val<1>(maxInd) += 1;
   if(geoMaxX2 >  val<2>(blockCoord))
      val<2>(maxInd) += 1;
   if(geoMaxX3 >  val<3>(blockCoord))
      val<3>(maxInd) += 1;

   this->fillExtentWithBlocks(grid);

   grid->setNX1(val<1>(maxInd));
   grid->setNX2(val<2>(maxInd));
   grid->setNX3(val<3>(maxInd));
}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::findOrigin( Grid3DPtr grid )
{
   orgX1 = boundingBox->getX1Minimum();
   orgX2 = boundingBox->getX2Minimum();
   orgX3 = boundingBox->getX3Minimum();

   //double minX1, minX2, minX3;

   //minX1 = boundingBox->getX1Minimum();
   //minX2 = boundingBox->getX2Minimum();
   //minX3 = boundingBox->getX3Minimum();

   //if(minX1 <= orgX1) orgX1 = minX1;
   //if(minX2 <= orgX2) orgX2 = minX2;
   //if(minX3 <= orgX3) orgX3 = minX3;
}
