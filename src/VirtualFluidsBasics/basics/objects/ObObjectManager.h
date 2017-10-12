//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBOBJECTMANAGER_H
#define OBOBJECTMANAGER_H

#include <string>
#include <sstream>
#include <vector>

#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbObserver.h>

class UbException;
class UbTableModel;
class ObObjectManager;
class ObObjectFactory;
class ObObject;


class ObObjectEntry
{
   friend class ObObjectManager;
public:
   ObObjectManager* getParent() { return parent; }
   ObObject*        getObject() { return object; }
   
   ObObjectEntry(ObObjectManager* parent, ObObject* object);
   virtual ~ObObjectEntry() {  }

protected:
   ObObjectManager* parent;
   ObObject* object;
};


class ObObjectManager : public UbObservable, public UbObserver
{
public:
	ObObjectManager();
	~ObObjectManager();
	
   //virtual bool addObObject(ObObject* object);   
   virtual bool addObObjectEntry(ObObjectEntry* objectEntry);

   virtual ObObjectEntry* createNewObObjectEntry(ObObject* obj) { return new ObObjectEntry(this, obj); }

	bool removeObObject(ObObject* object);
	bool removeObObject(int index);
	void removeAllObObjects();
	bool selectObObject(int index);
	bool selectObObject(ObObject* object);
	ObObject* getSelectedObObject();
	int getSelectedIndex();

	int getNumberOfObObjects();                 
   std::vector<ObObject*>* getAllObObjects();
   std::vector<ObObjectEntry*>* getAllObObjectEntries();
	ObObject* getObObject(int index);
   ObObjectEntry* getObObjectEntry(int index);

	std::string toString();

	virtual void objectChanged(UbObservable* observable);
	virtual void objectWillBeDeleted(UbObservable* observable);

	UbTableModel* getTableModel();
   virtual ObObjectFactory* getObObjectFactory()=0;

protected:
	 std::vector<ObObjectEntry*> objectList;
	 ObObject* selectedObject;
	 UbTableModel* tableModel;
};

#endif //OBOBJECTMANAGER_H
