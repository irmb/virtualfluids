#include <basics/objects/ObObjectManager.h>
#include <basics/objects/ObObject.h>
#include <basics/objects/ObObjectCreator.h>
#include <basics/utilities/UbTableModel.h>
#include <basics/utilities/UbException.h>

using namespace std;

ObObjectEntry::ObObjectEntry(ObObjectManager *parent, ObObject *object)
{
   this->parent = parent;
   this->object = object;
}
/*======================================================*/
ObObjectManager::ObObjectManager()
{
	this->selectedObject = NULL;
	this->tableModel = NULL;
}

/*======================================================*/
ObObjectManager::~ObObjectManager()
{
	//cerr<<"NEIN, notifyObserversObjectWillBeDeleted wird AUSSCHLIESSLICH von BasisKlasse aufgerufen!!!"<<endl;
 //  cerr<<"das muss so sein, denn ansonsten duerfte diese funktion nur in der speziellen klasse stehen, da\n";
 //  cerr<<"virtuelle destruktoren sich rekursiv vom speziellen ins allg. aufrufen --> notify.. wuerde\n";
 //  cerr<<"oefters aufgerufen werden...\n";

	this->objectList.clear();
	if(this->tableModel) delete this->tableModel;
}
/*======================================================*/
UbTableModel* ObObjectManager::getTableModel()
{ 
	return tableModel; 
}
/*======================================================*/
//bool ObObjectManager::addObObject(ObObject *object)
//{
//   cout<<"ObObjectManager::addObObject "<<object->toString()<<endl;
//	for(int pos=0; pos<(int)this->objectList.size(); pos++)
//		if(this->objectList[pos]->object==object) 
//			return false;
//
//	this->objectList.push_back(new ObObjectEntry(this,object));
//	//object->addObserver(this);
//	this->selectObObject(object);
//	return true;
//}
/*======================================================*/
bool ObObjectManager::addObObjectEntry(ObObjectEntry* objectEntry)
{
   for(int pos=0; pos<(int)this->objectList.size(); pos++)
      if(this->objectList[pos]->object==objectEntry->object) 
         return false;

   this->objectList.push_back(objectEntry);
//   objectEntry->getObject()->addObserver(this);
   this->selectObObject(objectEntry->object);
   return true;
}
/*======================================================*/
bool ObObjectManager::removeObObject(ObObject* object)
{
	if (this->selectedObject == object) this->selectedObject=NULL;
	for(int pos=0; pos<(int)this->objectList.size(); pos++)
	{

		if(this->objectList[pos]->object==object) 
		{
         return this->removeObObject(pos);
//			this->objectList.erase(objectList.begin()+pos);
//			//this->removeObserver(this);
//			return true;
		}
	}
	return false;
}
/*======================================================*/
bool ObObjectManager::removeObObject(int index)
{
	try
	{
		if ( objectList[index]->object == this->selectedObject ) this->selectedObject=NULL;
      //den entry loeschen ... das object im Entry ??? erstmal ausserhalb ...
      delete objectList[index]; 
		objectList.erase(objectList.begin()+index);
   	this->notifyObserversObjectChanged();
   	return true;
	}
	catch(const std::exception& e)  {  cerr<<e.what()<<endl;    }
   catch(...)                      {  cerr<<"Fehler in ObObjectManager::removeObObject"<<endl; }
   return false;
}
/*======================================================*/
void ObObjectManager::removeAllObObjects() 
{  
	//TODO: implementieren!!
	//foreach grid:
	//grid->removeObserver(this);
	//vector<ObObject*>::iterator it;
	//for(it=objectList.begin();  it!=objectList.end(); it++)
	//{
	//	it->removeObserver(this);
	//}
// 	for(int i=0; i<(int)objectList.size(); i++)
// 	{
// 		delete objectList[i]->object->removeObserver(this);
// 	} 
	this->objectList.clear();
	this->selectedObject = NULL;
	this->notifyObserversObjectChanged();
}
/*======================================================*/
int ObObjectManager::getNumberOfObObjects()
{ 
	return (int)this->objectList.size();
}
/*======================================================*/
vector<ObObject*>* ObObjectManager::getAllObObjects()  
{ 
   UB_THROW( UbException(UB_EXARGS,"hier muss noch was getan werden") );
//	return this->objectList;  
}
vector<ObObjectEntry*>* ObObjectManager::getAllObObjectEntries()
{
   return &this->objectList;  
}
/*======================================================*/
ObObject* ObObjectManager::getObObject(int index)
{
	if(index <  0)                            return NULL;
	if(index >= (int)this->objectList.size()) return NULL;

	return(this->objectList[index]->object);
}
/*======================================================*/
ObObjectEntry* ObObjectManager::getObObjectEntry(int index)
{
   if(index <  0)                            return NULL;
   if(index >= (int)this->objectList.size()) return NULL;

   return(this->objectList[index]);
}
/*====================================================*/
string ObObjectManager::toString()
{
	stringstream ss; ss<<endl;

	for(int pos=0; pos<(int)this->objectList.size(); pos++)          
	{
		ObObject* object = this->objectList[pos]->object;
		ss<<(pos+1)<<". "<<object->toString()<<endl;
	}
	return ss.str();
}
/*======================================================*/
void ObObjectManager::objectChanged(UbObservable* observable)
{
   //cout<<"ObObjectManager::objectChanged ??";
	this->notifyObserversObjectChanged();
}
/*======================================================*/
void ObObjectManager::objectWillBeDeleted(UbObservable* observable)
{
   cout<<"ObObjectManager::objectWillBeDeleted ??";
	//observable->removeObserver(this);
}
/*======================================================*/
bool ObObjectManager::selectObObject(int index)
{
   if((int)this->objectList.size()==0) 
   {
      this->selectedObject = NULL; return false; 
   }
	if (index > (int)this->objectList.size()-1 || index < 0) return false; 
	if ( this->selectedObject == this->getObObject(index) ) return true;
   
	this->selectedObject = this->getObObject(index);
   //cout<<this->getObserverList()->size()<<endl;

	this->notifyObserversObjectChanged();
	return true;
}
/*======================================================*/
bool ObObjectManager::selectObObject(ObObject* object)
{
   if((int)this->objectList.size()==0) { this->selectedObject = NULL; return false; }
	for(int pos=0; pos<(int)this->objectList.size(); pos++)
	{
		if(this->objectList[pos]->object==object) 
		{
			return this->selectObObject(pos);
		}
	}
	return false;
}
/*======================================================*/
ObObject* ObObjectManager::getSelectedObObject()
{
	return this->selectedObject;
}
/*======================================================*/
int ObObjectManager::getSelectedIndex()
{
	for(int pos=0; pos<(int)this->objectList.size(); pos++)
	{
		if(this->objectList[pos]->object==this->selectedObject) 
		{
			return pos;
		}
	}
	return -1;
}
/*======================================================*/

