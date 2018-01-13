#ifndef INTERACTOR3D_H
#define INTERACTOR3D_H

#include <vector>
#include <PointerDefinitions.h>

#include "UbSystem.h"
#include "UbTuple.h"

class Block3D;
class Grid3D;
class UbFileInput;
class UbFileOutput;
class GbObject3D;
class Block3D;


class Interactor3D : public enableSharedFromThis<Interactor3D>
{
public:
   enum Accuracy{SIMPLE, EDGES, FACES, POINTS};
   Interactor3D();
   Interactor3D(SPtr<Grid3D> grid, int type=Interactor3D::SOLID);
   Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type);
   //! constructor
   //! \param a set accuracy for arePointsInObject() and arePointsNotInObject()
   Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type, Interactor3D::Accuracy a);
   
   virtual ~Interactor3D();
   virtual void initInteractor(const double& timestep=0); 
   virtual void updateInteractor(const double& timestep=0)=0;

   void setSolidBlock(SPtr<Block3D> block);
   void setBCBlock(SPtr<Block3D> block);

    virtual UbTupleDouble3 getForces();

   void setSolid()        { UbSystem::setBit(this->type, SOLID   ); }
   void setMoveable()     { UbSystem::setBit(this->type, MOVEABLE); }
   
   bool isSolid()         { return UbSystem::bitCheck(this->type, SOLID        ); }
   bool isInverseSolid()  { return UbSystem::bitCheck(this->type, INVERSESOLID ); }
   bool isTimeDependent() { return UbSystem::bitCheck(this->type, TIMEDEPENDENT); }
   bool isMoveable()      { return UbSystem::bitCheck(this->type, MOVEABLE     ); }
   
   SPtr<Grid3D> getGrid3D()  const { return grid.lock();   }
   void setGrid3D(SPtr<Grid3D> grid) { this->grid = grid; }
   virtual SPtr<GbObject3D>  getGbObject3D() const { return geoObject3D; }
   virtual bool setDifferencesToGbObject3D(const SPtr<Block3D> block/*, const double& x1, const double& x2, const double& x3, const double& blockLengthX1, const double& blockLengthX2, const double& blockLengthX3, const double& timestep=0*/)
   {
      return false;  
   }

   virtual std::vector<SPtr<Block3D> >& getBcBlocks() { return this->bcBlocks; }
   virtual void removeBcBlocks() { this->bcBlocks.clear(); }
   virtual std::vector<SPtr<Block3D> >& getSolidBlockSet() { return this->solidBlocks; }
   virtual void removeSolidBlocks() { this->solidBlocks.clear(); }

protected:
   void setTimeDependent()   { UbSystem::setBit(this->type  , TIMEDEPENDENT); }
   void unsetTimeDependent() { UbSystem::unsetBit(this->type, TIMEDEPENDENT); }
   
   //! detect that points are inside object
   //! \param min/max coordinates of bounding box
   //! \param delta is delta x
   bool arePointsInsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);
   
   //! detect that points aren't inside object
   //! \param min/max coordinates of bounding box
   //! \param delta is delta x
   bool arePointsOutsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);

   //! detect that points are cutting object
   //! \param min/max coordinates of bounding box
   //! \param delta is delta x
   bool arePointsCuttingGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);
   
   bool isBlockOutsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);
   bool isBlockInsideGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);
   bool isBlockCuttingGeoObject(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3, double delta);

   int type;
   
   WPtr<Grid3D> grid;
   SPtr<GbObject3D> geoObject3D;

   std::vector<SPtr<Block3D> > bcBlocks;
   std::vector<SPtr<Block3D> > solidBlocks;
   int accuracy;

public:
   static const int SOLID	            ;//= (1<<0); //1
   static const int INVERSESOLID       ;//= (1<<1); //2
   static const int TIMEDEPENDENT      ;//= (1<<2); //4   //zeitlich
   static const int FLUID              ;//= (1<<3); //8
   static const int MOVEABLE           ;//= (1<<4); //16  // geometrisch
   static const int CHANGENOTNECESSARY ;//= (1<<5); //32

private:

};



#endif
