#ifndef GBOBJECT3DCREATOR_H
#define GBOBJECT3DCREATOR_H

#include <string>

#include <basics/objects/ObObjectCreator.h>
#include <basics/utilities/UbAutoRun.hpp>

#include <numerics/geometry3d/GbObject3D.h>

#ifdef CAB_QT 
#include <qdialog.h>
#endif

#ifdef CAB_VTK
#include <userinterface/presentation/vtkPoElement3D.h>
#endif

#ifdef CAB_PARAVIEW 
#include "vtkPVSource.h"
#endif          

class GbObject3DCreator : public ObObjectCreator                           
{                                       
protected:
   GbObject3DCreator() {}
private:
   GbObject3DCreator( const GbObject3DCreator& );                  //no copy allowed !!!
   const GbObject3DCreator& operator=( const GbObject3DCreator& ); //no copy allowed
public:
   virtual ~GbObject3DCreator(){}

   virtual std::string getTypeID() { return getGbObject3DTypeID();}
   virtual ObObject* createObObject()
   {
      return this->createGbObject3D();
   }


   virtual GbObject3D* createGbObject3D()=0;
   virtual std::string getGbObject3DTypeID()=0;                       
   virtual std::string toString() { return "GbObject3DCreator"; }     

#ifdef CAB_QT 
   virtual GbObject3D* createGbObject3DwithQt(QWidget* parent=0, Qt::WFlags flags=0) { throw UbException(UB_EXARGS,"Not implemented..."); }
   virtual void editGbObject3DwithQt(GbObject3D* gbObj, QWidget* parent=0, Qt::WFlags flags=0)  { throw UbException(UB_EXARGS,"Not implemented..."); }
#endif
   //die Teile von ObObjectCreator ...
#ifdef CAB_QT 
   void showSpecificInstrument(ObObject* object,QWidget* parent=0)
   {
      GbObject3D* geoObj = dynamic_cast<GbObject3D*>(object);
      this->editGbObject3DwithQt(geoObj, parent);
   }
   virtual ObObject* createObObjectWithQt() { return this->createGbObject3DwithQt();}
   virtual QObObjectSpecificInstrument* getSpecificInstrument() { throw UbException(UB_EXARGS,"not implemented"); }

#endif
#ifdef CAB_VTK 
   virtual Presentator* createObjectPresentator(ObObject *object) { return NULL; }
#endif


#ifdef CAB_PARAVIEW 
   virtual vtkPVSource* createPVSource(vtkPVWindow *Window) {  throw UbException(UB_EXARGS,"vtkPVSource* createPVSource"); }
#endif


};

#include <numerics/geometry3d/creator/GbObject3DFactory.h>

/*=========================================================================*/
#endif



