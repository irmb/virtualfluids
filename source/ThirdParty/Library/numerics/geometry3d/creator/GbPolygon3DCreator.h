#ifndef GBPOLYGON3DCREATOR_H
#define GBPOLYGON3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbPoint3D.h>

class GbPolygon3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbPolygon3DCreator* getInstance()
   {
      static GbPolygon3DCreator instance;
      return &instance;
   }

   GbPolygon3D* createGbObject3D() { return new GbPolygon3D(); }

   std::string getGbObject3DTypeID() { return "GbPolygon3D";        }
   std::string toString()            { return "GbPolygon3DCreator"; }

private:
   GbPolygon3DCreator( const GbPolygon3DCreator& );                  //no copy allowed 
   const GbPolygon3DCreator& operator=( const GbPolygon3DCreator& ); //no copy allowed

   GbPolygon3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbPolygon3DCreator::getInstance()), CAB_GbPolygon3DCreator);
#endif

#endif
