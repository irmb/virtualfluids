#ifndef GBLINE3DCREATOR_H
#define GBLINE3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbLine3D.h>

class GbLine3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbLine3DCreator* getInstance()
   {
      static GbLine3DCreator instance;
      return &instance;
   }

   GbLine3D* createGbObject3D() { return new GbLine3D(); }

   std::string getGbObject3DTypeID(){ return "GbLine3D";        }
   std::string toString()           { return "GbLine3DCreator"; }

private:
   GbLine3DCreator( const GbLine3DCreator& );                  //no copy allowed 
   const GbLine3DCreator& operator=( const GbLine3DCreator& ); //no copy allowed
   GbLine3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbLine3DCreator::getInstance()), CAB_GbLine3DCreator);
#endif

#endif
