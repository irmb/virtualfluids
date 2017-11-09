//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBOBJECTCREATOR_H
#define OBOBJECTCREATOR_H

#include <string>

class ObObject;
class ObObjectManager;

class Presentator;
class QViewer;

#ifdef CAB_QT 
class QObObjectSpecificInstrument;
class QWidget;
class QActionGroup;
#endif

class ObObjectCreator
{
public:
   virtual ~ObObjectCreator() {}

	virtual ObObject* createObObject()=0;

	virtual std::string getTypeID()	{ return "ObObject"; }
	virtual std::string toString()	{ return "ObObjectCreator"; }
   
#ifdef CAB_QT 
   //virtual Presentator* createObjectPresentator(ObObject *object)=0;
   virtual Presentator* createObjectPresentator(ObObject *object) { return NULL; }
   virtual QActionGroup* getSpecificPresentatorGroup(ObObject* object, QViewer *viewer, QWidget* parent) { return NULL; }
   virtual QActionGroup* getSpecificActionGroup(ObObjectManager* manager, ObObject* object, QWidget* parent) 
   { 
      return NULL; 
   }

   virtual ObObject* createObObjectWithQt() { return NULL; }
   virtual void showSpecificInstrument(ObObject* object, QWidget* parent=0) {}
   virtual QObObjectSpecificInstrument* getSpecificInstrument() { return NULL; }
   
   //virtual QActionGroup *getSpecificContextMenuActionGroup() { return NULL; }
#endif

protected:
	ObObjectCreator() {}

private:
   ObObjectCreator( const ObObjectCreator& );                  //no copy allowed 
   const ObObjectCreator& operator=( const ObObjectCreator& ); //no copy allowed

};
#endif
