#ifndef GBOBJECT3DFACTORY_H
#define GBOBJECT3DFACTORY_H

#include <string>
#include <sstream>
#include <map>

#include <basics/objects/ObObjectFactory.h>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

class GbObject3D;
class GbObject3DCreator;

class GbObject3DFactory : public ObObjectFactory
{
private:
   GbObject3DFactory();
   GbObject3DFactory( const GbObject3DFactory& );                  //no copy allowed 
   const GbObject3DFactory& operator=( const GbObject3DFactory& ); //no copy allowed
public:
   static GbObject3DFactory* getInstance();
   
   GbObject3D* createGbObject3D(UbFileInput* in);

   //void addGbObject3DCreator(GbObject3DCreator* creator);
   //void deleteGbObject3DCreator(GbObject3DCreator* creator);
   //std::map<std::string, GbObject3DCreator*>* getCreatorSet() { return &creatorSet;   }

   std::string toString();
   GbObject3D* createEmptyGbObject3D(std::string objectType);
   //GbObject3DCreator* getCreator(std::string objectTypeID);

private:
   
   
   //std::map<std::string, GbObject3DCreator*> creatorSet;
   //typedef std::map<std::string, GbObject3DCreator*>::iterator CreatorIterator;
};
/*=========================================================================*/
#endif

