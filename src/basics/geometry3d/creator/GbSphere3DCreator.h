#ifndef GBSPHERE3DCREATOR_H
#define GBSPHERE3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>
#include <geometry3d/GbSphere3D.h>

#ifdef CAB_QT 
#include <geometry3d/presentation/QGbSphere3DInstrument.h>
#include <QtGui/QWidget>
#include <QtGui/QDialog>
#endif

#ifdef CAB_VTK
#include <geometry3d/presentation/vtkGbSphere3D.h>
#endif

#ifdef CAB_PARAVIEW 
#include "vtkSMSourceProxy.h"
#include "vtkSMProperty.h"
#include "vtkSMDoubleVectorProperty.h"
#endif

class GbSphere3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbSphere3DCreator* getInstance()
   {
      static GbSphere3DCreator instance;
      return &instance;
   }

   GbSphere3D* createGbObject3D() { return new GbSphere3D(); }

   std::string getGbObject3DTypeID() { return "GbSphere3D"; };
   std::string toString()            { return "GbSphere3DCreator"; }

private:
   GbSphere3DCreator( const GbSphere3DCreator& );                  //no copy allowed 
   const GbSphere3DCreator& operator=( const GbSphere3DCreator& ); //no copy allowed
   GbSphere3DCreator() : GbObject3DCreator() {}

#ifdef CAB_QT 
public:

   GbSphere3D* createGbObject3DwithQt(QWidget* parent=0, Qt::WFlags flags=0)
   { 
      GbSphere3D* sphere = this->createGbObject3D();
      sphere->setRadius(3.0);
      sphere->setCenterX1Coordinate(6.0);

      QGbSphere3DInstrument instrument(parent, flags);
      instrument.setGbSphere3D(sphere);
      if (instrument.exec()) { return sphere; }
      delete sphere;
      return NULL;
   }
   QDialog* getSpecificInstrument(QWidget* parent=0, Qt::WFlags flags=0) { return new QGbSphere3DInstrument(parent, flags); }

   void editGbObject3DwithQt(GbObject3D* gbObj, QWidget* parent=0, Qt::WFlags flags=0)
   { 
      GbSphere3D* sphere = dynamic_cast<GbSphere3D*>(gbObj);
      if(!sphere) throw UbException(UB_EXARGS,"selected object to edit is no GbSphere3D");

      QGbSphere3DInstrument instrument(parent, flags);
      instrument.setGbSphere3D(sphere);
      instrument.exec();
   }
#endif
#ifdef CAB_VTK
public:
   Presentator* createObjectPresentator(ObObject *object) { return new vtkGbSphere3D(dynamic_cast<GbSphere3D*>(object)); }
#endif
  

#ifdef CAB_PARAVIEW
   vtkPVSource* createPVSource(vtkPVWindow *Window);
#endif
};

#ifdef CAB_PARAVIEW                  
vtkPVSource* GbSphere3DCreator::createPVSource(vtkPVWindow *Window)
{
   GbSphere3D *mysphere = this->createGbObject3D();
   mysphere->setCenterX1Coordinate(2.0);
   mysphere->setCenterX2Coordinate(1.0);
   mysphere->setCenterX3Coordinate(3.0);
   mysphere->setRadius(0.3);

   vtkPVSource* pvs = Window->CreatePVSource("SphereSource");
   pvs->CreateProperties();
   if (pvs)
   {
      pvs->IsPermanentOn();
      pvs->Accept(1, 1);
      //pvs->SetTraceReferenceObject(this->GetWindow());
      pvs->SetTraceReferenceObject(Window);
   }
   //vtkPVDisplayGUI *settingsGUI= pvs->GetPVOutput();

   vtkSMSourceProxy* proxy = pvs->GetProxy();
   vtkSMProperty *prop = proxy->GetProperty("Center");

   vtkSMDoubleVectorProperty *doubleprop = vtkSMDoubleVectorProperty::SafeDownCast(proxy->GetProperty("Center"));
   doubleprop->SetElement(0, mysphere->getX1Centroid());
   doubleprop->SetElement(1, mysphere->getX2Centroid());
   doubleprop->SetElement(2, mysphere->getX3Centroid());
   pvs->SetLabel("Kugel");

   doubleprop = vtkSMDoubleVectorProperty::SafeDownCast(proxy->GetProperty("Radius"));
   doubleprop->SetElement(0, mysphere->getRadius());

   pvs->GetPVWidget("Center")->ResetInternal();
   pvs->GetPVWidget("Radius")->ResetInternal();

   pvs->SetVisibility(TRUE);
   pvs->AcceptCallback();
   pvs->Update();
   return pvs;
}
#endif

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbSphere3DCreator::getInstance()), CAB_GbSphere3DCreator);
#endif

#endif
