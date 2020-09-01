#ifndef GBCUBOID3DCREATOR_H
#define GBCUBOID3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>
#include <geometry3d/GbCuboid3D.h>

#ifdef CAB_QT 
#include <geometry3d/presentation/QGbCuboid3DInstrument.h>
#include <QtGui/QWidget>
#include <QtGui/QDialog>
#endif

#ifdef CAB_VTK 
#include <geometry3d/presentation/vtkGbCuboid3D.h>
#endif

class GbCuboid3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbCuboid3DCreator* getInstance()
   {
      static GbCuboid3DCreator instance;
      return &instance;
   }

   GbCuboid3D* createGbObject3D() { return new GbCuboid3D(); }          

   std::string getGbObject3DTypeID() { return "GbCuboid3D"; };
   std::string toString()            { return "GbCuboid3DCreator"; }

private:
   GbCuboid3DCreator() : GbObject3DCreator() {}

   GbCuboid3DCreator( const GbCuboid3DCreator& );                  //no copy allowed 
   const GbCuboid3DCreator& operator=( const GbCuboid3DCreator& ); //no copy allowed

#ifdef CAB_QT 
public:
   GbCuboid3D* createGbObject3DwithQt(QWidget* parent=0, Qt::WFlags flags=0)
   {                                                              
      GbCuboid3D* cuboid = this->createGbObject3D();
      cuboid->getPoint2()->setX1(2.0);
      cuboid->getPoint2()->setX2(2.0);
      cuboid->getPoint2()->setX3(2.0);

      QGbCuboid3DInstrument instrument(parent, flags);
      instrument.setGbCuboid3D(cuboid);
      if (instrument.exec()) { return cuboid; }
      delete cuboid;

      return NULL;
   }

   QDialog* getSpecificInstrument(QWidget* parent=0, Qt::WFlags flags=0) { return new QGbCuboid3DInstrument(parent, flags); }

   void editGbObject3DwithQt(GbObject3D* gbObj, QWidget* parent=0, Qt::WFlags flags=0)
   { 
      GbCuboid3D* cuboid = dynamic_cast<GbCuboid3D*>(gbObj);
      if(!cuboid) throw UbException(UB_EXARGS,"selected object to edit is no GbCuboid3D!");

      QGbCuboid3DInstrument instrument(parent, flags);
      instrument.setGbCuboid3D(cuboid);
      instrument.exec();
   }
#endif
#ifdef CAB_VTK
public:
   Presentator* createObjectPresentator(ObObject *object) { return new vtkGbCuboid3D(dynamic_cast<GbCuboid3D*>(object)); }
#endif

};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbCuboid3DCreator::getInstance()), CAB_GbCuboid3DCreator);
#endif

#endif   
