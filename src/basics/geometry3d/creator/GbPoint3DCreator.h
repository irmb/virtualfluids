#ifndef GBPOINT3DCREATOR_H
#define GBPOINT3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>
#include <geometry3d/GbPoint3D.h>

class GbPoint3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbPoint3DCreator* getInstance()
   {
      static GbPoint3DCreator instance;
      return &instance;
   }
   
   GbPoint3D* createGbObject3D() { return new GbPoint3D(); }
   
   std::string getGbObject3DTypeID() { return "GbPoint3D";        }
   std::string toString()            { return "GbPoint3DCreator"; }

private:
   GbPoint3DCreator( const GbPoint3DCreator& );                  //no copy allowed 
   const GbPoint3DCreator& operator=( const GbPoint3DCreator& ); //no copy allowed
   GbPoint3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbPoint3DCreator::getInstance()), CAB_GbPoint3DCreator);
#endif

#endif
