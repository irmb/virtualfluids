#ifndef FETRIFACEMESH3DCREATOR_H
#define FETRIFACEMESH3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/fem/FeTriFaceMesh3D.h>
#include <basics/utilities/UbFileInputASCII.h>

#ifdef CAB_QT 
#include <qfiledialog.h>    
#endif

#ifdef CAB_VTK
//#include <numerics/geometry3d/presentation/vtkGbTriangularMesh3D.h>
#endif

class FeTriFaceMesh3DCreator : public GbObject3DCreator              
{                                       
public:
   static FeTriFaceMesh3DCreator* getInstance()
   {
      static FeTriFaceMesh3DCreator instance;
      return &instance;
   }
   static FeTriFaceMesh3D* readMeshFromFile(std::string filename, std::string meshName, bool removeRedundantNodes=true);

   static FeTriFaceMesh3D* readMeshFromPLYFile(std::string filename, std::string meshName, bool removeRedundantNodes=true);
   static FeTriFaceMesh3D* readMeshFromPLYFile(UbFileInput* in, std::string meshName, bool removeRedundantNodes=true);

   static FeTriFaceMesh3D* readMeshFromSTLFile(std::string filename, std::string meshName, bool removeRedundantNodes=true); 
   static FeTriFaceMesh3D* readMeshFromSTLFile(UbFileInput* in, std::string meshName, bool removeRedundantNodes=true);

   static FeTriFaceMesh3D* readMeshFromAVSFile(std::string filename, std::string meshName, bool removeRedundantNodes=true); 
   static FeTriFaceMesh3D* readMeshFromAVSFile(UbFileInput* in, std::string meshName, bool removeRedundantNodes=true);

   static FeTriFaceMesh3D* readMeshFromVTKASCIIFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true); 
   static FeTriFaceMesh3D* readMeshFromVTKASCIIFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);


   FeTriFaceMesh3D* createGbObject3D() { return new FeTriFaceMesh3D(); }
   
   std::string getGbObject3DTypeID(){ return "FeTriFaceMesh3D"; };
   std::string toString()           { return "FeTriFaceMesh3DCreator"; }

#ifdef CAB_QT 


   FeTriFaceMesh3D* createGbObject3DwithQt()
   {
	   //QString s = QFileDialog::getOpenFileName(NULL,NULL,NULL,"open file dialog","Choose a STL file" );
	   QString s = QFileDialog::getOpenFileName(NULL, "Choose a STL file", "/home", "*.stl");
      //QFileDialog* fd = new QFileDialog( NULL );
      //fd->setIconText(QString("Hallo"));
      //fd->show();
      //TODO: Open File Dialog einbauen.		
      UbFileInputASCII in( s.toAscii().data() );
      stringstream stream;
      stream <<"TriangularMesh3D ";//<<_objCount++;
      FeTriFaceMesh3D *mesh = NULL;//FeTriFaceMesh3DCreator::readMeshFromSTLFile(&in, stream.str() );
      return mesh;
   }
   //QDialog* getSpecificInstrument()  {  return 0;}
   void editGbObject3DwithQt(GbObject3D* gbObj)
   { 
   }
#endif
#ifdef CAB_VTK
public:
   Presentator* createObjectPresentator(ObObject *object) { return new vtkGbTriangularMesh3D((GbTriangularMesh3D*)object); }
#endif


private:
   FeTriFaceMesh3DCreator( const FeTriFaceMesh3DCreator& );                  //no copy allowed 
   const FeTriFaceMesh3DCreator& operator=( const FeTriFaceMesh3DCreator& ); //no copy allowed
   FeTriFaceMesh3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(FeTriFaceMesh3DCreator::getInstance()), CAB_FeTriFaceMesh3DCreator);
#endif

#endif
