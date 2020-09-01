#ifndef GBCYLINDER3DCREATOR_H
#define GBCYLINDER3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>
#include <geometry3d/GbCylinder3D.h>

#ifdef CAB_QT 
#include <geometry3d/presentation/QGbCylinder3DInstrument.h>
#include <QtGui/QWidget>
#include <QtGui/QDialog>
#endif

#ifdef CAB_VTK
#include <geometry3d/presentation/vtkGbCylinder3D.h>
#endif


class GbCylinder3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbCylinder3DCreator* getInstance()
   {
      static GbCylinder3DCreator instance;
      return &instance;
   }

   GbCylinder3D* createGbObject3D() { return new GbCylinder3D(); }
   
   std::string getGbObject3DTypeID(){ return "GbCylinder3D";        }
   std::string toString()           { return "GbCylinder3DCreator"; }

private:
   GbCylinder3DCreator( const GbCylinder3DCreator& );                  //no copy allowed 
   const GbCylinder3DCreator& operator=( const GbCylinder3DCreator& ); //no copy allowed
GbCylinder3DCreator() : GbObject3DCreator() {}

#ifdef CAB_QT
public:

   GbCylinder3D* createGbObject3DwithQt(QWidget* parent=0, Qt::WFlags flags=0)
   {                                                              
      GbCylinder3D* cylinder = this->createGbObject3D();
      cylinder->setRadius(2.0);
      cylinder->setPoint1(0.0, 0.0, 0.0);
      cylinder->setPoint2(0.0, 5.0, 0.0);

      QGbCylinder3DInstrument instrument(parent, flags);
      instrument.setGbCylinder3D(cylinder);
      if (instrument.exec()){ return cylinder; }
      delete cylinder;

      return NULL;
   }

   QDialog* getSpecificInstrument(QWidget* parent=0, Qt::WFlags flags=0)
   { 
      return new QGbCylinder3DInstrument(parent, flags);
   }

   void editGbObject3DwithQt(GbObject3D* gbObj, QWidget* parent=0, Qt::WFlags flags=0)
   { 
      GbCylinder3D* cylinder = dynamic_cast<GbCylinder3D*>(gbObj);
      if(!cylinder) throw UbException(UB_EXARGS,"selected object to edit is no GbCylinder3D!");

      QGbCylinder3DInstrument instrument(parent, flags);
      instrument.setGbCylinder3D(cylinder);
      instrument.exec();
   }
#endif
#ifdef CAB_VTK
public:
   Presentator* createObjectPresentator(ObObject *object) { return new vtkGbCylinder3D(dynamic_cast<GbCylinder3D*>(object)); }
#endif

};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbCylinder3DCreator::getInstance()), CAB_GbCylinder3DCreator);
#endif

#endif   
