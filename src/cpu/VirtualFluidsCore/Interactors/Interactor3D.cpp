#include "Interactor3D.h"



#include <fstream>
#include <geometry3d/GbCuboid3D.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbFileOutput.h>
#include "UbException.h"

#include "Grid3D.h"
#include "Block3D.h"
#include "GbObject3D.h"


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
Interactor3D::Interactor3D(SPtr<Grid3D> grid, int type)
   :   type(type)
     , grid(grid)
{

}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type)
   :   geoObject3D(geoObject3D)
     , grid(grid)
     , type(type)
     , accuracy(SIMPLE)
{

}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type, Interactor3D::Accuracy a)
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
void Interactor3D::setSolidBlock(SPtr<Block3D> block)
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
void Interactor3D::setBCBlock(SPtr<Block3D> block)
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
      this->bcBlocks.push_back(block);
}

UbTupleDouble3 Interactor3D::getForces()
{
    UB_THROW( UbException("UbTupleDouble3 getForces() - gehoert in die abgeleitete klasse") );
}
void Interactor3D::setID(int id)
{
   this->id = id;
}
//////////////////////////////////////////////////////////////////////////
int Interactor3D::getID()
{
   return id;
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setActive()
{
   active = true;
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setInactive()
{
   active = false;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isActive()
{
   return active;
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::initInteractor(const double& timeStep)
{
   //UBLOG(logINFO, "transBlocks.size = "<<transBlocks.size());

   for(SPtr<Block3D> block : bcBlocks)
   {
      this->setDifferencesToGbObject3D(block);
   }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::updateInteractor(const double& timeStep)
{
   UB_THROW( UbException("Interactor3D::updateInteractor - toDo") );
}
//////////////////////////////////////////////////////////////////////////

