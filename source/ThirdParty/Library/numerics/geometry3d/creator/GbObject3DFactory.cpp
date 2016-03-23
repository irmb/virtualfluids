#include <numerics/geometry3d/creator/GbObject3DFactory.h>
#include <numerics/geometry3d/creator/GbObject3DCreator.h>
//#include <numerics/geometry3d/creator/GbPoint3DCreator.h>
// #include <numerics/geometry3d/creator/GbCuboid3DCreator.h>
// #include <numerics/geometry3d/creator/GbSphere3DCreator.h>
// #include <numerics/geometry3d/creator/GbCylinder3DCreator.h>
// #include <numerics/geometry3d/creator/GbLine3DCreator.h>
// #include <numerics/geometry3d/creator/GbPolygon3DCreator.h>
// #include <numerics/geometry3d/creator/GbTriangle3DCreator.h>
// #include <numerics/geometry3d/creator/GbTriangularMesh3DCreator.h>

using namespace std;

//GbObject3DFactory* GbObject3DFactory::instance = NULL;

/*======================================================================*/
GbObject3DFactory::GbObject3DFactory() 
   : ObObjectFactory()
{
}
/*======================================================================*/  
GbObject3DFactory* GbObject3DFactory::getInstance()
{
   static GbObject3DFactory instance;
   return &instance;
}

///*======================================================================*/
//void GbObject3DFactory::addGbObject3DCreator(GbObject3DCreator* creator)
//{
//   //cout<<"Meth:"<<creator->toString()<<" Meth-ID:"<<creator->getGbObject3DTypeID()<<endl;
//   creatorSet.insert(pair<string, GbObject3DCreator*>(creator->getGbObject3DTypeID(), creator));
//}
//
//void GbObject3DFactory::deleteGbObject3DCreator(GbObject3DCreator* creator)
//{
//   throw UbException(UB_EXARGS,"GbObject3DFactory::deleteGbObject3DCreator not yet implemented");
//   // this.creatorSet.delete(creator);
//}

/*======================================================================*/
GbObject3D* GbObject3DFactory::createGbObject3D(UbFileInput *in) 
{
   string str = in->readString();
   //cout<<"GbObject3DFactory::createGbObject3D:"<<str<<endl;

   GbObject3D *gbObject3D = createEmptyGbObject3D(str);

   if(!gbObject3D)
      throw UbException(UB_EXARGS,"creator for type available");
   
   gbObject3D->read(in);

   return gbObject3D;
}
/*======================================================================*/
GbObject3D* GbObject3DFactory::createEmptyGbObject3D(string objectType)
{
   typedef std::map<string, ObObjectCreator*>::iterator CreatorIterator;
   std::map<string, ObObjectCreator*>* creatorSet = this->getCreatorSet();
   CreatorIterator creatorIterator = creatorSet->find(objectType);

   if(creatorIterator == creatorSet->end()) 
      throw UbException(UB_EXARGS,"factory has no creator for "+objectType);

   GbObject3DCreator *creator = dynamic_cast<GbObject3DCreator*>(creatorIterator->second);

   if(!creator) 
      throw UbException(UB_EXARGS,"factory has no creator for "+objectType);

   return creator->createGbObject3D();
}
/*======================================================================*/
//GbObject3DCreator* GbObject3DFactory::getCreator(string objectType)
//{
//   CreatorIterator creatorIterator = creatorSet.find(objectType);
//   if(creatorIterator == creatorSet.end()) throw UbException(UB_EXARGS,"factory has no creator for "+objectType);
//   GbObject3DCreator *creator = creatorIterator->second;
//   if(!creator) throw UbException(UB_EXARGS,"factory has no creator for "+objectType);
//   return(creator);
//}
/*======================================================================*/
string GbObject3DFactory::toString() 
{
   stringstream ss;
   ss<<"GbObject2DFactory";
   int a=1;
   std::map<std::string, ObObjectCreator*>::iterator creatorIterator; 
   std::map<std::string, ObObjectCreator*>* tmp = this->getCreatorSet();
   for(creatorIterator=tmp->begin(); creatorIterator!=tmp->end(); creatorIterator++)
   {
      ss<<(a++)<<". ";
      ss<<creatorIterator->second->getTypeID();
      ss<<endl;
   }
   return(ss.str());
}
