#ifndef GBTRIANGLE3DCREATOR_H
#define GBTRIANGLE3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbTriangle3D.h>

class GbTriangle3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbTriangle3DCreator* getInstance()
   {
      static GbTriangle3DCreator instance;
      return &instance;
   }

   GbTriangle3D* createGbObject3D() { return new GbTriangle3D(); }
   
   std::string getGbObject3DTypeID(){ return "GbTriangle3D";        }
   std::string toString()           { return "GbTriangle3DCreator"; }

private:
   GbTriangle3DCreator( const GbTriangle3DCreator& ); //no copy allowed 
   const GbTriangle3DCreator& operator=( const GbTriangle3DCreator& ); //no copy allowed
   GbTriangle3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbTriangle3DCreator::getInstance()), CAB_GbTriangle3DCreator);
#endif

#endif
