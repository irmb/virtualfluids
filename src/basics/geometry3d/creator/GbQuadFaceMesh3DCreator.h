#ifndef GBQUADFACEMESH3DCREATOR_H
#define GBQUADFACEMESH3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>

class GbQuadFaceMesh3D;

#ifdef CAB_QT 

#endif

#ifdef CAB_VTK
//#include <geometry3d/presentation/vtkGbQuadangularMesh3D.h>
#endif

class GbQuadFaceMesh3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbQuadFaceMesh3DCreator* getInstance()
   {
      static GbQuadFaceMesh3DCreator instance;
      return &instance;
   }
   static GbQuadFaceMesh3D *createQuadMesh3D(int nodesX1, int nodesX2, float startX1, float startX2, double knotenabstandX1, double knotenabstandX2, float nullNiveau, std::string name);

   GbObject3D* createGbObject3D();
   
   std::string getGbObject3DTypeID() { return "GbQuadFaceMesh3D";        }
   std::string toString()            { return "GbQuadFaceMesh3DCreator"; }

#ifdef CAB_QT 

   GbQuadFaceMesh3D* createGbObject3DwithQt();
   //QDialog* getSpecificInstrument()  {  return 0;}
   void editGbObject3DwithQt(GbObject3D* gbObj);
#endif
#ifdef CAB_VTK
   Presentator* createObjectPresentator(ObObject *object);
#endif


private:
   GbQuadFaceMesh3DCreator( const GbQuadFaceMesh3DCreator& );                  //no copy allowed 
   const GbQuadFaceMesh3DCreator& operator=( const GbQuadFaceMesh3DCreator& ); //no copy allowed
   GbQuadFaceMesh3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbQuadFaceMesh3DCreator::getInstance()), CAB_GbQuadFaceMesh3DCreator);
#endif

#endif
