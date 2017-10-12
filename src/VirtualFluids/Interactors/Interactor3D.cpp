#include "Interactor3D.h"

#include <boost/foreach.hpp>

#include <fstream>
#include <numerics/geometry3d/GbCuboid3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbFileOutput.h>

using namespace std;

const int Interactor3D::SOLID	           = (1<<0); //1
const int Interactor3D::INVERSESOLID       = (1<<1); //2
const int Interactor3D::TIMEDEPENDENT      = (1<<2); //4   //zeitlich
const int Interactor3D::FLUID              = (1<<3); //8
const int Interactor3D::MOVEABLE           = (1<<4); //16  // geometrisch
const int Interactor3D::CHANGENOTNECESSARY = (1<<5); //32

//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D()
  : type(SOLID)
{

}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(Grid3DPtr grid, int type)
   :   type(type)
     , grid(grid)
{
}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(GbObject3DPtr geoObject3D, Grid3DPtr grid, int type)
   :   geoObject3D(geoObject3D)
     , grid(grid)
     , type(type)
     , accuracy(SIMPLE)
{
}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(GbObject3DPtr geoObject3D, Grid3DPtr grid, int type, Interactor3D::Accuracy a)
   :   geoObject3D(geoObject3D)
   , grid(grid)
   , type(type)
   , accuracy(a)
{
}
//////////////////////////////////////////////////////////////////////////
Interactor3D::~Interactor3D()
{
}
//////////////////////////////////////////////////////////////////////////
//void Interactor3D::deleteSolidBlocks(int level)
//{
//   //hier werden die Bloecke aktiv oder nicht aktiv gesetzt
//   double minX1,minX2,minX3,maxX1,maxX2,maxX3,x1,x2,x3;
//   int gridRank = grid.lock()->getRank();
//
//   vector<Block3DPtr> blockVector;
//   bool activ = true;
//   grid.lock()->getBlocks(level, gridRank, activ, blockVector);
//   BOOST_FOREACH(Block3DPtr block, blockVector)
//   {
//      double deltaX = grid.lock()->getDeltaX(block);
//      UbTupleDouble3 blockLengths  = grid.lock()->getBlockLengths(block);
//
//      //Koords bestimmen
//      UbTupleDouble3 org = grid.lock()->getBlockWorldCoordinates(block);
//
//      x1 = val<1>(org);
//      x2 = val<2>(org);
//      x3 = val<3>(org);
//
//      minX1 = x1;
//      minX2 = x2;
//      minX3 = x3;
//      maxX1 = x1 + val<1>(blockLengths);
//      maxX2 = x2 + val<2>(blockLengths);
//      maxX3 = x3 + val<3>(blockLengths);
//
//      if(this->isInverseSolid())
//      {
//         switch (accuracy)
//         {
//         //simple duff
//         case SIMPLE:
//            if(!this->geoObject3D->isCellInsideOrCuttingGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3))
//               block->setActive(false);
//            break;
//         //test only edges
//         case EDGES:
//            if(arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, deltaX) &&
//
//               arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, deltaX) &&
//               arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, deltaX) &&
//
//               arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, deltaX))   
//                  block->setActive(false);
//            break;
//         //test only faces
//         case FACES:
//            if(arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, deltaX))
//                  block->setActive(false);
//            break;
//         //test all points
//         case POINTS:
//            if(arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
//               block->setActive(false);
//            break;
//         default:
//            UB_THROW( UbException(UB_EXARGS, "Accuracy isn't correct") );
//            break;
//         }
//      }
//      else //solid 
//      {
//         switch (accuracy)
//         {
//         //simple duff
//         case SIMPLE:
//            if(this->geoObject3D->isCellInsideGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3))
//               block->setActive(false);
//            break;
//         //test only edges
//         case EDGES:
//            if(arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, deltaX) &&
//
//               arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, deltaX) &&
//               arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, deltaX) &&
//
//               arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, deltaX))   
//               block->setActive(false);
//            break;
//         //test only faces
//         case FACES:
//            if(arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, deltaX) &&
//               arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, deltaX))
//               block->setActive(false);
//            break;
//         //test all points
//         case POINTS:
//            if(arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
//               block->setActive(false);
//            break;
//         default:
//            UB_THROW( UbException(UB_EXARGS, "Accuracy isn't correct") );
//            break;
//         }
//      }
//   }
//}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsInsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   bool result = true;
   for (double ix3=minX3; ix3<=maxX3; ix3+=delta)
      for (double ix2=minX2; ix2<=maxX2; ix2+=delta)
         for (double ix1=minX1; ix1<=maxX1; ix1+=delta)
            result = result && this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3);

   return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsOutsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   bool result = true;
   for (double ix3=minX3; ix3<=maxX3; ix3+=delta)
      for (double ix2=minX2; ix2<=maxX2; ix2+=delta)
         for (double ix1=minX1; ix1<=maxX1; ix1+=delta)
            result = result && (!this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3));

   return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsCuttingGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   bool result = true;
   for (double ix3=minX3; ix3<=maxX3; ix3+=delta)
      for (double ix2=minX2; ix2<=maxX2; ix2+=delta)
         for (double ix1=minX1; ix1<=maxX1; ix1+=delta)
            result = result || this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3);

   return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockOutsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   switch (accuracy)
   {
      //simple duff
   case SIMPLE:
      return !this->geoObject3D->isCellInsideOrCuttingGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3);
      //test only edges
   case EDGES:
      return arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) &&
             arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
             arePointsOutsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) &&
             
             arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) &&
             arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
             arePointsOutsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) &&
             
             arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
             arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) &&
             arePointsOutsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);   
      //test only faces
   case FACES:
      return arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
             arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) &&
             arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
      //test all points
   case POINTS:
      return arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
   default:
      UB_THROW( UbException(UB_EXARGS, "Accuracy isn't correct") );
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockInsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   switch (accuracy)
   {
      //simple duff
   case SIMPLE:
      return this->geoObject3D->isCellInsideGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3);
      //test only edges
   case EDGES:
      return arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) &&
             arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
             arePointsInsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) &&
             
             arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) &&
             arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
             arePointsInsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) &&
             
             arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
             arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) &&
             arePointsInsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);   
      //test only faces
   case FACES:
      return arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
             arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) &&
             arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) &&
             arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
             arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
      //test all points
   case POINTS:
      return arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
   default:
      UB_THROW( UbException(UB_EXARGS, "Accuracy isn't correct") );
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockCuttingGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta)
{
   switch (accuracy)
   {
      //simple duff
   case SIMPLE:
      return this->geoObject3D->isCellCuttingGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3);
      //test only edges
   case EDGES:
      return arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) ||
             arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) ||
             arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) ||
             arePointsCuttingGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) ||
                                                                             
             arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) ||
             arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) ||
             arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) ||
             arePointsCuttingGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) ||
                                                                             
             arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) ||
             arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) ||
             arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) ||
             arePointsCuttingGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);   
      //test only faceCutting
   case FACES:        
      return arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) ||
             arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) ||
             arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) ||
             arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) ||
             arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) ||
             arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
      //test all pointCutting
   case POINTS:       
      return arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
   default:
      UB_THROW( UbException(UB_EXARGS, "Accuracy isn't correct") );
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setSolidBlock(Block3DPtr block)
{
   double minX1,minX2,minX3,maxX1,maxX2,maxX3;

   double deltaX = grid.lock()->getDeltaX(block);
   UbTupleDouble3 blockLengths  = grid.lock()->getBlockLengths(block);
   UbTupleDouble3 org = grid.lock()->getBlockWorldCoordinates(block);
   UbTupleDouble3 nodeOffset = grid.lock()->getNodeOffset(block);

   //coordinates of block without ghost layer
   minX1 = val<1>(org) + val<1>(nodeOffset);
   minX2 = val<2>(org) + val<2>(nodeOffset);
   minX3 = val<3>(org) + val<3>(nodeOffset);
   maxX1 = val<1>(org) + val<1>(blockLengths) - val<1>(nodeOffset);
   maxX2 = val<2>(org) + val<2>(blockLengths) - val<2>(nodeOffset);
   maxX3 = val<3>(org) + val<3>(blockLengths) - val<3>(nodeOffset);

   if(this->isInverseSolid())
   {
      if(isBlockOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
      {
         block->setActive(false);
         this->solidBlocks.push_back(block);
      }
   }
   else //solid 
   {
      if(isBlockInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
      {
         block->setActive(false);
         this->solidBlocks.push_back(block);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setTransBlock(Block3DPtr block)
{
   double minX1,minX2,minX3,maxX1,maxX2,maxX3;

   double deltaX = grid.lock()->getDeltaX(block);
   UbTupleDouble3 blockLengths  = grid.lock()->getBlockLengths(block);
   UbTupleDouble3 org = grid.lock()->getBlockWorldCoordinates(block);
   UbTupleDouble3 nodeOffset = grid.lock()->getNodeOffset(block);

   //coordinates of block with ghost layer
   minX1 = val<1>(org) - val<1>(nodeOffset);
   minX2 = val<2>(org) - val<2>(nodeOffset);
   minX3 = val<3>(org) - val<3>(nodeOffset);
   maxX1 = val<1>(org) + val<1>(blockLengths) + val<1>(nodeOffset);
   maxX2 = val<2>(org) + val<2>(blockLengths) + val<2>(nodeOffset);
   maxX3 = val<3>(org) + val<3>(blockLengths) + val<3>(nodeOffset);

   if(isBlockCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
      this->transBlocks.push_back(block);
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::initInteractor(const double& timeStep)
{
   //UBLOG(logINFO, "transBlocks.size = "<<transBlocks.size());

   BOOST_FOREACH(Block3DPtr block, transBlocks)
   {
      this->setDifferencesToGbObject3D(block);
   }
}
//////////////////////////////////////////////////////////////////////////
//SOLID:
//bloecke werden nicht activ gesetzt, wenn der Block vollstaendig in der Geometrie
//blocke werden in transBlock hinzugefuegt, wenn block+delta schnittpunkt oder beruehrungspunkt mit geo hat
//fuer jeden transblock wird "setDifferencesToGbObject3D aufgerufen
//void Interactor3D::initInteractor(const double& timeStep)
//{
   //this->removeTransBlocks();
   //this->removeSolidBlocks();

   ////hier werden die Bloecke aktiv oder nicht aktiv gesetzt
   //double minX1,minX2,minX3,maxX1,maxX2,maxX3,x1,x2,x3;
   //int gridRank = grid.lock()->getRank();

   //int minInitLevel = this->grid.lock()->getCoarsestInitializedLevel();
   //int maxInitLevel = this->grid.lock()->getFinestInitializedLevel();

   //for(int level = minInitLevel; level<=maxInitLevel;level++)
   //{
   //   vector<Block3DPtr> blockVector;
   //   grid.lock()->getBlocks(level, gridRank, blockVector);
   //   BOOST_FOREACH(Block3DPtr block, blockVector)
   //   {
   //      //Koords bestimmen
   //      UbTupleDouble3 org = grid.lock()->getBlockWorldCoordinates(block);
   //      UbTupleDouble3 blockLengths  = grid.lock()->getBlockLengths(block);
   //      double dx = grid.lock()->getDeltaX(block);
   //      UbTupleDouble3 orgDelta = grid.lock()->getNodeOffset(block);
   //      UbTupleDouble3 coords = grid.lock()->getNodeCoordinates(block, 0, 0, 0);

   //      //x1 = val<1>(org);
   //      //x2 = val<2>(org);
   //      //x3 = val<3>(org);
   //      x1 = val<1>(coords);
   //      x2 = val<2>(coords);
   //      x3 = val<3>(coords);

   //      minX1 = x1;
   //      minX2 = x2;
   //      minX3 = x3;
   //      //maxX1 = x1 + val<1>(blockLengths);
   //      //maxX2 = x2 + val<2>(blockLengths);
   //      //maxX3 = x3 + val<3>(blockLengths);
   //      maxX1 = val<1>(coords);
   //      maxX2 = val<2>(coords);
   //      maxX3 = val<3>(coords);

   //      if(this->isInverseSolid())
   //      {
   //         if(   UbMath::lessEqual(minX1,geoObject3D->getX1Minimum()) 
   //            && UbMath::lessEqual(minX2,geoObject3D->getX2Minimum())
   //            && UbMath::lessEqual(minX3,geoObject3D->getX2Minimum())
   //            && UbMath::greaterEqual(maxX1,geoObject3D->getX1Maximum())
   //            && UbMath::greaterEqual(maxX2,geoObject3D->getX2Maximum()) 
   //            && UbMath::greaterEqual(maxX3,geoObject3D->getX2Maximum()))
   //         {
   //            this->transBlocks.push_back(block);
   //            this->setDifferencesToGbObject3D(block/*, x1, x2, x3, val<1>(blockLengths), val<2>(blockLengths), val<3>(blockLengths), timeStep*/);
   //         }
   //         else if(this->geoObject3D->isCellCuttingGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3))
   //         {
   //            this->transBlocks.push_back(block);
   //            this->setDifferencesToGbObject3D(block/*, x1, x2, x3, val<1>(blockLengths), val<2>(blockLengths), val<3>(blockLengths), timeStep*/);
   //         }
   //      }
   //      else //solid 
   //      {
   //         //1. Fall: block komlett in geo ist in deleteSolidBlocks() erledigt
   //         //2. Fall:  Celle umhuellt Geo oder Cell schneidet Geo
   //         if( this->geoObject3D->isCellCuttingGbObject3D(minX1,minX2,minX3,maxX1,maxX2,maxX3) )
   //         {
   //            this->transBlocks.push_back(block);
   //            this->setDifferencesToGbObject3D(block/*,x1,x2,x3,val<1>(blockLengths),val<2>(blockLengths),val<3>(blockLengths),timeStep*/);
   //         }
   //         //3. Fall: cell umhuellt geo:
   //         else if(   UbMath::lessEqual(    minX1, geoObject3D->getX1Minimum() )
   //            && UbMath::lessEqual(    minX2, geoObject3D->getX2Minimum() )
   //            && UbMath::lessEqual(    minX3, geoObject3D->getX3Minimum() )
   //            && UbMath::greaterEqual( maxX1, geoObject3D->getX1Maximum() )
   //            && UbMath::greaterEqual( maxX2, geoObject3D->getX2Maximum() )
   //            && UbMath::greaterEqual( maxX3, geoObject3D->getX3Maximum() ) ) //block umhuellt geo
   //         {
   //            throw UbException(UB_EXARGS,"//3. Fall: cell umhuellt geo sollte mit Fall 2 abgedeckt sein!!!");

   //            this->transBlocks.push_back(block);
   //            this->setDifferencesToGbObject3D(block/*,x1,x2,x3,val<1>(blockLengths),val<2>(blockLengths),val<3>(blockLengths),timeStep*/);
   //         }
   //      }
   //   }
   //}
//}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::updateInteractor(const double& timeStep)
{
   UB_THROW( UbException("Interactor3D::updateInteractor - toDo") );
}
//////////////////////////////////////////////////////////////////////////

