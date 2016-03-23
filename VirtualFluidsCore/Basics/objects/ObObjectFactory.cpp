#include <basics/objects/ObObjectFactory.h>

/**** Eigene ****/
#include <basics/objects/ObObjectCreator.h>
#include <basics/utilities/UbException.h>

using namespace std;

//ObObjectFactory::ObObjectFactory()
//{
//}
//
//ObObjectFactory::~ObObjectFactory()
//{
//}
/*======================================================================*/  
//ObObjectFactory* ObObjectFactory::getInstance()
//{
//	static ObObjectFactory instance;
//	return &instance;
//}
/*======================================================================*/
void ObObjectFactory::addObObjectCreator(ObObjectCreator *creator)
{
	//cout<<"Meth:"<<creator->toString()<<" Meth-ID:"<<creator->getTypeID()<<endl;
	creatorSet.insert(std::pair<string, ObObjectCreator*>(creator->getTypeID(), creator));
}
/*======================================================================*/
void ObObjectFactory::removeObObjectCreator(ObObjectCreator *creator)
{
	UB_THROW( UbException(UB_EXARGS,"not implemented") );
}
/*======================================================================*/
ObObjectCreator* ObObjectFactory::getCreator(string objectType) 
{
	std::map<string, ObObjectCreator*>::iterator creatorIterator = creatorSet.find(objectType);
	if(creatorIterator == creatorSet.end()) UB_THROW( UbException(UB_EXARGS,"factory has no creator for "+objectType) );
	ObObjectCreator *creator = creatorIterator->second;
	if(!creator) UB_THROW( UbException(UB_EXARGS,"no time series creator for type available") );
	return creator;
}
/*======================================================================*/
string ObObjectFactory::toString() 
{
   stringstream text;

   std::map<string, ObObjectCreator*>::iterator creatorIterator;
   std::map<string, ObObjectCreator*>* creatorSet = this->getCreatorSet();

   for(creatorIterator = creatorSet->begin(); creatorIterator!=creatorSet->end(); ++creatorIterator)
      text<<"   - "<<(creatorIterator->second)->toString()<<" for "<<(creatorIterator->first)<<endl;

   return text.str();
}
