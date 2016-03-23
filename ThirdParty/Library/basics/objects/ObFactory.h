//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBFACTORY_H
#define OBFACTORY_H


#include <string>
#include <map>
#include <sstream>
#include <iomanip>
#include <typeinfo>

#include <basics/objects/ObCreator.h>

/*=========================================================================*/
/*  ObFactory                                                            */
/*                                                                         */
/**
generic factory
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 14.06.07
@version 1.1 - 12.04.08
*/ 

/*
usage:  T       = zu erzeugende Klasse
        Creator = Erzeugerklasse
//////////////////////////////////////////////////////////////////////////
//example
//////////////////////////////////////////////////////////////////////////
//  class Base{ 
//  public:
//        OBCREATOR_EXT(Base)
//  };
//  //automatisches registrieren:
//  UB_AUTO_RUN_NAMED(ObFactory<Base>::getInstance()->addObCreator(ObCreatorImpl<Base,Base>::getInstance()), CAB_Base);
//  class Derived : public Base 
//  {
//   public:
//        OBCREATOR_EXT(Derived)
//};
//  //automatisches registrieren:
//  UB_AUTO_RUN_NAMED(ObFactory<Base>::getInstance()->addObCreator(ObCreatorImpl<Base,Derived>::getInstance()), CAB_Derived);
////////////////////////////////////////////////////////////////////////////
//  int main()
//  {
//       //Alternativ zu UB_AUTO_RUN_NAMED: haendisches registrieren
//       ObFactory<Base>::getInstance()->addObCreator(ObCreatorImpl<Base>::getInstance());
//       ObFactory<Base>::getInstance()->addObCreator(ObCreatorImpl<Derived,Base>::getInstance());
// 
//       //create objects - method1
//       Base* test1 = ObFactory<Base>::getInstance()->createObject<Base>();
//       Base* test2 = ObFactory<Base>::getInstance()->createObject<Derived>();
// 
//       //create objects - method2
//       Base* test1 = ObFactory<Base>::getInstance()->createObject(Base::getStaticClassObjectTypeID()    );
//       Base* test2 = ObFactory<Base>::getInstance()->createObject(Derived::getStaticClassObjectTypeID() );
//   //...
// }
*/


template<class  T, typename Creator = ObCreator< T > >
class ObFactory
{
   typedef std::map<  std::string, Creator* > CreatorMap;
   typedef typename CreatorMap::iterator      CreatorMapIt;
   typedef std::pair< std::string, Creator* > CreatorMapElement;

protected:
   ObFactory() {}  //so ist vererbung gewahrleistet

private:
   ObFactory( const ObFactory< T, Creator >& );    //no copy allowed 
   const ObFactory& operator=( const ObFactory& ); //no copy allowed


public:
   virtual ~ObFactory() {}

   static ObFactory< T, Creator >* getInstance() 
   {
      static ObFactory< T, Creator > instance;
      return &instance;
   }

   bool addObCreator(Creator* creator);
   bool removeObCreator(Creator* creator);

   T* createObject(const std::string& objectTypeID);
   
   template< typename T2 > 
   T* createObject() { return this->createObject( T2::getStaticClassObjectTypeID() ); }
   
   Creator* getCreator(const std::string& objectTypeID);

   virtual std::string toString();
  
private:
   CreatorMap creatorMap;
};

//////////////////////////////////////////////////////////////////////////
//Implementation
template<class  T, typename Creator >
bool ObFactory< T, Creator >::addObCreator(Creator* creator)
{
	if(creatorMap.insert( CreatorMapElement(creator->getObjectTypeID(), creator) ).second )
   {
      //insert succeeded
      return true;
   }
   //insert fails
   return false;
}
/*======================================================================*/
template<class  T, typename Creator >
bool ObFactory< T, Creator >::removeObCreator(Creator* creator)
{
   if(creator && creatorMap->erase( creator->getClassObjectTypeID() ) ) 
      return true;

   return false;
}
/*======================================================================*/
template<class  T, typename Creator >
Creator* ObFactory< T, Creator >::getCreator(const std::string& obtypeID)
{
   CreatorMapIt it = creatorMap.find(obtypeID);
   if(it == creatorMap.end()) return NULL;

   Creator* creator = it->second;
   if(!creator) return NULL;

   return creator;
}
/*======================================================================*/
 template<class  T, typename Creator >
 T* ObFactory< T, Creator >::createObject(const std::string& objectTypeID)
 {
    Creator* creator = this->getCreator(objectTypeID);
    
    if(!creator) 
    {
       UB_THROW( UbException(UB_EXARGS,"no creator avaivlable for ID="+objectTypeID ) );
    }
 
    return creator->createObject();
 }
/*======================================================================*/
template<class  T, typename Creator >
std::string ObFactory< T, Creator >::toString() 
{
   std::size_t maxL = 6;
   for(CreatorMapIt it=creatorMap.begin(); it!=creatorMap.end(); ++it)
      if( it->first.size() > maxL ) 
         maxL = it->first.size();
   
   std::stringstream os;
   os<<(std::string)typeid(*this).name()<<" - info:"<<std::endl;
   os<<"   "<<std::left<<std::setw(maxL)<<"object"<<" <-> "<<"creator "<<std::endl;
   for(CreatorMapIt it=creatorMap.begin(); it!=creatorMap.end(); ++it)
      os<< " - " << std::setw(maxL) << it->first << " <-> " << (std::string)typeid(*it->second).name() << std::endl;

   return os.str();
}
/*======================================================================*/

#endif //OBFACTORY_H
