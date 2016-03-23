//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBCREATOR_H
#define OBCREATOR_H

#include <string>

/*=========================================================================*/
/*  ObCreator / ObCreatorImpl                                              */
/*                                                                         */
/**
generic factory
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 14.06.07
@version 1.1 - 12.04.08
*/ 

/*
usage: see bottom of file "./ObFactory.h"
*/

//////////////////////////////////////////////////////////////////////////
// ObCreator
// Um in der Factory verschiedene Typen von Creaors in einer 
// std::map<std::string,ObCreator<BaseT>*> halten zu koennen
// muss eine gemeinsame Basisklasse existieren
//////////////////////////////////////////////////////////////////////////
template< class BaseT >
class ObCreator
{
public:
   virtual std::string  getObjectTypeID()=0;
   virtual BaseT*       createObject() = 0;

   virtual ~ObCreator() {  }

protected:
   ObCreator() {}
private:
   ObCreator( const ObCreator< BaseT >& );         //no copy allowed 
   const ObCreator& operator=( const ObCreator& ); //no copy allowed
};

//////////////////////////////////////////////////////////////////////////
// ObCreatorImpl
// Implementierung des speziellen Creators 
//////////////////////////////////////////////////////////////////////////
template< class T, class BaseT=T >
class ObCreatorImpl : public ObCreator< BaseT >
{
public:
   static ObCreator<BaseT >* getInstance()
   {
      static ObCreatorImpl< T, BaseT > instance;
      return &instance;
   }

public:
   ~ObCreatorImpl() {}

   //aus portabilitaetsgruenden kann man nicht typeinfo nehmen, da diese compilerabhaengig ist
   std::string getObjectTypeID()  { return T::getStaticClassObjectTypeID();  } 
   
   virtual T*  createObject() { return new T(); }

protected:
	ObCreatorImpl() {}
private:
	ObCreatorImpl( const ObCreatorImpl< T, BaseT >& );      //no copy allowed 
   const ObCreatorImpl& operator=( const ObCreatorImpl& ); //no copy allowed
};

//////////////////////////////////////////////////////////////////////////
// ObCreatorImpl
// Implementierung des speziellen Creators fuer Singletons
//////////////////////////////////////////////////////////////////////////
template< class T, class BaseT=T >
class ObSingletonCreatorImpl : public ObCreator< BaseT >
{
public:
   static ObCreator<BaseT >* getInstance()
   {
      static ObSingletonCreatorImpl< T, BaseT > instance;
      return &instance;
   }
public:
   ~ObSingletonCreatorImpl() {}

   //aus portabilitaetsgruenden kann man nicht typeinfo nehmen, da diese compilerabhaengig ist
   std::string getObjectTypeID()  { return T::getStaticClassObjectTypeID();  } 

   virtual T* createObject() { return T::getInstance(); }

protected:
   ObSingletonCreatorImpl() {}
private:
   ObSingletonCreatorImpl( const ObSingletonCreatorImpl< T, BaseT >& );      //no copy allowed 
   const ObSingletonCreatorImpl& operator=( const ObSingletonCreatorImpl& ); //no copy allowed
};

//workaround for the not perfect C++ world. typeinfo::name is not usable for this purpose!
//see Andrei Alexandrescu, "Modern C++ Design: Generic Programming and Design Patterns Applied", Chapter 8.5
#define OBCREATOR_EXT( ClassObject ) \
   static  std::string  getStaticClassObjectTypeID() { return #ClassObject;                 } \
   virtual std::string  getClassObjectTypeID()       { return getStaticClassObjectTypeID(); } 

#endif //OBCREATOR_H
