//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBOBJECT3DMANAGER_H
#define GBOBJECT3DMANAGER_H

#include <string>
#include <sstream>
#include <vector>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>
#include <basics/utilities/UbFileOutput.h>
#include <basics/utilities/UbTableModel.h>
#include <basics/objects/ObObjectManager.h>

#include <basics/memory/MbSharedPointerDefines.h>
class GbObject3DManager;
typedef VFSharedPtr<GbObject3DManager> GbObject3DManagerPtr;

                                                  
class GbObject3D;   
class GbObject3DManager;
class GbObject3DTableModel;

class GbObject3DEntry  : public ObObjectEntry
{
   friend class GbObject3DManager;
   friend class GbObject3DTableModel;
public:
   std::string getName() { return this->name;}
private:
   //GbObject3DManager *parent;
   //GbObject3D        *geoObject;
   bool        active;
   std::string name;
   

   GbObject3DEntry(GbObject3DManager* parent, GbObject3D* geoObject, bool active, std::string name);
 
};


class GbObject3DManager  : public ObObjectManager 
{                                           
public:
   GbObject3DManager();
   ~GbObject3DManager();

   ObObjectEntry* createNewObObjectEntry(ObObject *obj);
   ObObjectFactory* getObObjectFactory();


   //bool addGbObject3D(GbObject3D *geoObject3D);  
   bool addGbObject3D(GbObject3D *geoObject3D, std::string name);   
   //bool addGbObject3D(GbObject3D *geoObject3D, bool active, std::string name);   

   bool removeGbObject3D(GbObject3D *geoObject3D);
   bool removeGbObject3D(int index);

   int getNumberOfGbObject3Ds();                 
   std::vector<GbObject3D*>* getAllGbObject3Ds();
   GbObject3D* getGbObject3D(int index);

   //keine Definition dafuer da ...
   //void writeValidationAVSFile(string filename);
   //void writeSurfaceAVSFile(string filename);

   //public final OctConstructionDescriptor[] getAllActiveConstructions()
   //public final OctSpecificConstructionInstrument getSpecificConstructionInstrumentInstance(int index)
   //public final boolean isConstructionActive(int index)
   //public final boolean isConstructionVisible(int index)
   //public final void activateConstruction(int index, boolean active)
   //public final void visibleConstruction(int index, boolean visible)
  // UbTableModel* getTableModel();
   //void objectChanged(UbObservable* observable);
   //void objectWillBeDeleted(UbObservable* observable);


   void read(UbFileInput *in);
   void write(UbFileOutput *out); 

   std::string toString();

private:
   //GbObject3DTableModel* tablemodel;
  // GbObject3D* selectedObject;
   //vector<GbObject3DEntry*> gbObject3DList;
};

class GbObject3DTableModel : public UbTableModel
{
public:

   static const int COL_NAME   = 0;
   static const int COL_TYPE   = 1;
   static const int COL_ACTIVE = 2;

   GbObject3DTableModel(GbObject3DManager* manager);
   ~GbObject3DTableModel(void);

   //////////////////////////////////////////////////////////////////////////
   //Geerbt von CabTable
   int getColumnNumber(void);
   int getRowNumber();

   std::string getColumnLabel(int column);

   int getColumnType(int);
   std::string getStringValue(int row, int col);
   int getSelectedRowIndex();
   //bool GetBoolValue(int row, int col);
   void setStringValue(int row, int col, std::string str) { throw UbException(UB_EXARGS,"not implemented"); }

protected:
   GbObject3DManager* objectManager;
};

#endif
