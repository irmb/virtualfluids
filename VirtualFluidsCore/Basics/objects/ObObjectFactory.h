//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBOBJECTFACTORY_H
#define OBOBJECTFACTORY_H

#include <string>
#include <map>

class ObObjectCreator; 

class ObObjectFactory
{
public:
   ObObjectFactory() {}
   virtual ~ObObjectFactory() {}

   //static geht nicht, da abgeleitete Factories existieren ...
   //static ObObjectFactory* getInstance();
   //virtual ObObjectFactory* getInstance()=0;

   ObObjectCreator* getCreator(std::string objectType);

	void addObObjectCreator(ObObjectCreator* creator);
	void removeObObjectCreator(ObObjectCreator* creator);

   std::map<std::string, ObObjectCreator*>* getCreatorSet() { return &creatorSet;  }

   virtual std::string toString();
	
private:
   ObObjectFactory( const ObObjectFactory& );                  //no copy allowed 
   const ObObjectFactory& operator=( const ObObjectFactory& ); //no copy allowed

   std::map<std::string, ObObjectCreator*> creatorSet;
};


#endif //OBOBJECTFACTORY_H
