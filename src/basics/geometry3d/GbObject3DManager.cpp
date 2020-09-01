#include <geometry3d/GbObject3DManager.h>
#include <geometry3d/GbObject3D.h>
#include <geometry3d/creator/GbObject3DFactory.h>

using namespace std;

GbObject3DEntry::GbObject3DEntry(GbObject3DManager *parent, GbObject3D *geoObject, bool active, string name):ObObjectEntry(parent, geoObject)
{
   //this->parent = parent;
   //this->geoObject = geoObject;
   this->active = active;
   this->name = name;
}

/*======================================================*/
/*==     nun halt der Manager                       ====*/
/*======================================================*/

GbObject3DManager::GbObject3DManager(): ObObjectManager()
{
   this->tableModel = new GbObject3DTableModel(this);
}
/*======================================================*/
GbObject3DManager::~GbObject3DManager()
{
//   this->gbObject3DList.clear();
}
/*======================================================*/
ObObjectFactory* GbObject3DManager::getObObjectFactory()
{
   return GbObject3DFactory::getInstance();
}

/*======================================================*/
ObObjectEntry* GbObject3DManager::createNewObObjectEntry(ObObject *obj)
{ 
   GbObject3D *geoobject = dynamic_cast<GbObject3D*>(obj);
   return new GbObject3DEntry(this, geoobject, true, geoobject->getTypeID()); 
}

/*======================================================*/
//UbTableModel* GbObject3DManager::getTableModel()
//{
//   return this->tablemodel;
//}

/*======================================================*/
//bool GbObject3DManager::addGbObject3D(GbObject3D *geoObject)
//{
//   return this->addGbObject3D(geoObject, true, "GeoObject");
//}
//
///*======================================================*/
bool GbObject3DManager::addGbObject3D(GbObject3D *geoObject, string name)
{
   GbObject3DEntry *entry = new GbObject3DEntry (this, geoObject, true, name);
   return ObObjectManager::addObObjectEntry(entry);
}
//bool GbObject3DManager::addGbObject3D(GbObject3D *geoObject, bool active, string name)  
//{
//   GbObject3DEntry *entry = new GbObject3DEntry (this, geoObject, true, name);
//   return ObObjectManager::addObObjectEntry(entry);
//}
//
/*======================================================*/
bool GbObject3DManager::removeGbObject3D(GbObject3D *geoObject)
{
   return ObObjectManager::removeObObject(geoObject);
}
/*======================================================*/
bool GbObject3DManager::removeGbObject3D(int index)
{
   return ObObjectManager::removeObObject(index);
}
/*======================================================*/
//void GbObject3DManager::removeAllGbObject3Ds() 
//{  
//    this->gbObject3DList.clear();
//}
/*======================================================*/
int GbObject3DManager::getNumberOfGbObject3Ds()
{ 
   return GbObject3DManager::getNumberOfObObjects();
}
/*======================================================*/
vector<GbObject3D*>* GbObject3DManager::getAllGbObject3Ds()  
{ 
   throw UbException(UB_EXARGS,"not implemented");
   //vector<GbObject3D*> *geoVektor = new vector<GbObject3D*>;
   //for(int u=0; u<(int)this->gbObject3DList.size();u++)
   //{
   //   GbObject3D* geoObject = dynamic_cast<GbObject3D*>((gbObject3DList)[u]->getObject());
   //   geoVektor->push_back(geoObject);
   //}
   //return geoVektor;  
}
/*======================================================*/
GbObject3D* GbObject3DManager::getGbObject3D(int index)
{
   if(index <  0)                            return NULL;
   if(index >= this->getNumberOfObObjects()) return NULL;

   GbObject3D* geoObject = dynamic_cast<GbObject3D*>(this->getObObject(index));
   return(geoObject);
}
/*======================================================*/
//GbObject3DEntry* GbObject3DManager::getGbObject3DEntry(int index)
//{
//   if(index <  0)                                 return NULL;
//   if(index >= (int)this->gbObject3DList.size())  return NULL;
//
//   return((this->gbObject3DList)[index]);
//}
/*====================================================*/
void GbObject3DManager::write(UbFileOutput *out) 
{                    
   int size = this->getNumberOfObObjects();
   out->writeInteger(size);
   out->writeString("// #GbObjects");

   GbObject3D *geoObject;
   for(int pos=0; pos<size; pos++)          
   {
      out->writeLine();
      geoObject = dynamic_cast<GbObject3D*>(this->getObObject(pos));
      cout<<pos<<".:"<<geoObject->toString()<<endl;
      geoObject->write(out);
   }
}
/*======================================================================*/
void GbObject3DManager::read(UbFileInput *in) 
{
   this->removeAllObObjects();
   
   int n = in->readInteger();                          
   
   cout<<"GbObject3DManager::read "<<n<<" GbObject3Ds\n";
   GbObject3D *geoObject;
   for(int pos=1; pos<=n; pos++)
   {
      in->readLine();
      cout<<" - GbObject3D "<<pos<<" ...";
      geoObject = GbObject3DFactory::getInstance()->createGbObject3D(in);
      
      GbObject3DEntry *entry = new GbObject3DEntry(this, geoObject, true, "GeoObject");
      this->addObObjectEntry(entry);
      cout<<"done\n";
   }
}
/*======================================================*/
string GbObject3DManager::toString()
{
   stringstream ss; ss<<endl;
   
   int size = this->getNumberOfObObjects();
   for(int pos=0; pos<size; pos++)          
   {
      ObObject* geoObject = this->getObObject(pos);
      ss<<(pos+1)<<". "<<geoObject->toString()<<endl;
   }
   string back = ss.str();
   return back;
}

/*======================================================*/
/*======================================================*/
/*======================================================*/

GbObject3DTableModel::GbObject3DTableModel(GbObject3DManager* manager)
{
   this->objectManager = manager;
}

/*======================================================*/
GbObject3DTableModel::~GbObject3DTableModel(void)
{
}

/*======================================================*/
//Gibt die Anzahl der Spalten zurueck.
int GbObject3DTableModel::getColumnNumber()
{
   return 3;
}

/*======================================================*/
std::string GbObject3DTableModel::getColumnLabel(int column)
{
   switch(column)
   {
   case COL_NAME: return "Name";
   case COL_TYPE: return "Type";
   case COL_ACTIVE: return "Active";
   default: throw UbException(UB_EXARGS,"falscher Spaltenindex");
   }
}

/*======================================================*/
int GbObject3DTableModel::getRowNumber()
{
   return this->objectManager->getNumberOfGbObject3Ds();
}
/*======================================================*/
int GbObject3DTableModel::getSelectedRowIndex()
{
   return 0;
   //	return this->objectManager->getSelectedIndex();
}

/*======================================================*/
int GbObject3DTableModel::getColumnType(int column)
{
   switch(column) {
   case COL_NAME :
      return UbTableModel::COL_TYPE_STRING;
      break;
   case COL_TYPE :
      return UbTableModel::COL_TYPE_STRING;
      break;
   case COL_ACTIVE :
      return UbTableModel::COL_TYPE_BOOL;
      break;
   }
   return -1;
}

//Gibt den Eintag der Tabelle an der angegebenen Spalten- und 
//Zeilenposition in Form eines String Werts zurueck.
std::string GbObject3DTableModel::getStringValue(int row, int col)
{
   GbObject3DEntry* gbObjEntry = dynamic_cast<GbObject3DEntry*>(this->objectManager->getObObjectEntry(row));
   switch(col) {
   case COL_NAME:
      return gbObjEntry->name;
      break;
   case COL_TYPE:
      return gbObjEntry->getObject()->getTypeID();
      break;
   case COL_ACTIVE:
      if ( gbObjEntry->active ) return "True";
      return "False";
      break;
   }
   return "Fehler";
}

/*======================================================*/
//bool GbObject3DManager::selectGbObject3D(int index)
//{
//   cout<<"GbObject3DManager::selectGbObject3D(int index):"<<index<<endl;
//   if (index > (int)this->gbObject3DList.size()-1 || index < 0) return false; 
//   if ( this->selectedObject == this->getGbObject3D(index) ) return true;
//   this->selectedObject = this->getGbObject3D(index);
//   this->notifyObserversObjectChanged();
//   return true;
//}
///*======================================================*/
//bool GbObject3DManager::selectGbObject3D(GbObject3D* geoObject)
//{
//   for(int pos=0; pos<(int)this->gbObject3DList.size(); pos++)
//   {
//      if(this->gbObject3DList[pos]->geoObject==geoObject) 
//      {
//         return this->selectGbObject3D(pos);
//      }
//   }
//   return false;
//}
///*======================================================*/
//GbObject3D* GbObject3DManager::getSelectedGbObject3D()
//{
//   return this->selectedObject;
//}
///*======================================================*/
//int GbObject3DManager::getSelectedIndex()
//{
//   for(int pos=0; pos<(int)this->gbObject3DList.size(); pos++)
//   {
//      if(this->gbObject3DList[pos]->geoObject==this->selectedObject) 
//      {
//         return pos;
//      }
//   }
//   return -1;
//}
