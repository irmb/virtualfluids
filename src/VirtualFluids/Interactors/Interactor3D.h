#ifndef INTERACTOR3D_H
#define INTERACTOR3D_H

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <cmath>

class UbFileInput;
class UbFileOutput;
class GbObject3D;
class Block3D;

#include <boost/serialization/shared_ptr.hpp>
class Interactor3D;
typedef boost::shared_ptr<Interactor3D> Interactor3DPtr;

#include "UbException.h"
#include "UbTuple.h"
#include "ObObject.h"
#include "GbObject3D.h"
#include "Grid3D.h"

#include <boost/serialization/serialization.hpp>
#include <boost/enable_shared_from_this.hpp>

class Interactor3D 
{
public:
   enum Accuracy{SIMPLE, EDGES, FACES, POINTS};
   Interactor3D();
   Interactor3D(Grid3DPtr grid, int type=Interactor3D::SOLID);
   Interactor3D(GbObject3DPtr geoObject3D, Grid3DPtr grid, int type);
   //! constructor
   //! \param a set accuracy for arePointsInObject() and arePointsNotInObject()
   Interactor3D(GbObject3DPtr geoObject3D, Grid3DPtr grid, int type, Interactor3D::Accuracy a);
   
   virtual ~Interactor3D();
   virtual void initInteractor(const double& timestep=0); 
   virtual void updateInteractor(const double& timestep=0)=0;
   //virtual void deleteSolidBlocks(int level);

   void setSolidBlock(Block3DPtr block);
   void setTransBlock(Block3DPtr block);
      
   virtual UbTupleDouble3 getForces() { UB_THROW( UbException("UbTupleDouble3 getForces() - gehoert in die abgeleitete klasse") ); }

   void setSolid()        { UbSystem::setBit(this->type, SOLID   ); }
   void setMoveable()     { UbSystem::setBit(this->type, MOVEABLE); }
   
   bool isSolid()         { return UbSystem::bitCheck(this->type, SOLID        ); }
   bool isInverseSolid()  { return UbSystem::bitCheck(this->type, INVERSESOLID ); }
   bool isTimeDependent() { return UbSystem::bitCheck(this->type, TIMEDEPENDENT); }
   bool isMoveable()      { return UbSystem::bitCheck(this->type, MOVEABLE     ); }
   
   Grid3DPtr getGrid3D()  const { return grid.lock();   }
   void setGrid3D(Grid3DPtr grid) { this->grid = grid; }
   virtual GbObject3DPtr  getGbObject3D() const { return geoObject3D; }
   virtual bool setDifferencesToGbObject3D(const Block3DPtr block/*, const double& x1, const double& x2, const double& x3, const double& blockLengthX1, const double& blockLengthX2, const double& blockLengthX3, const double& timestep=0*/)
   {
      return false;  
   }

   virtual std::vector<Block3DPtr>& getTransBlockSet() { return this->transBlocks; }
   virtual void removeTransBlocks() { this->transBlocks.clear(); }
   virtual std::vector<Block3DPtr>& getSolidBlockSet() { return this->solidBlocks; }
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
   
   Grid3DWeakPtr grid;
   GbObject3DPtr geoObject3D;

   std::vector<Block3DPtr> transBlocks;
   std::vector<Block3DPtr> solidBlocks;
   int accuracy;

public:
   static const int SOLID	            ;//= (1<<0); //1
   static const int INVERSESOLID       ;//= (1<<1); //2
   static const int TIMEDEPENDENT      ;//= (1<<2); //4   //zeitlich
   static const int FLUID              ;//= (1<<3); //8
   static const int MOVEABLE           ;//= (1<<4); //16  // geometrisch
   static const int CHANGENOTNECESSARY ;//= (1<<5); //32

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & transBlocks;
      ar & solidBlocks;
   }

};



#endif
