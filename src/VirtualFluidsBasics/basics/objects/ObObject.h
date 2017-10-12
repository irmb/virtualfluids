//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBOBJECT_H
#define OBOBJECT_H

#include <string>

#include <basics/objects/ObObjectCreator.h>
#include <basics/utilities/UbObservable.h>

#ifdef CAB_RCF
#include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif


class ObObjectCreator;

class ObObject : public UbObservable
{
public:
   ObObject() : name("") { }
   ObObject(const std::string& name) : name(name) { }

   virtual ~ObObject() { }

   virtual ObObject*   clone()=0;
   virtual std::string getTypeID()=0;

   virtual std::string getName()  { return name; }
   void setName(std::string name) { this->name=name; }

   virtual std::string toString()=0;

   virtual ObObjectCreator* getCreator()=0;

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar) 
   {
      //SF::SF_SERIALIZE_PARENT<UbObservable>(ar, *this);
      SF_SERIALIZE_PARENT<UbObservable>(ar, *this);
      ar & name;
   }
#endif //CAB_RCF

private:
   std::string name;
};

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
SF_NO_CTOR(ObObject);
UB_AUTO_RUN_NAMED( ( SF::registerType<ObObject>("ObObject") ),                SF_ObObject     );
UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived<UbObservable, ObObject >() ), SF_ObObject_BD1 );
#endif //RCF_USE_SF_SERIALIZATION

#endif
