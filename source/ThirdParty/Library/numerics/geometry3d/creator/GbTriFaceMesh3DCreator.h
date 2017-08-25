#ifndef GBTRIFACEMESH3DCREATOR_H
#define GBTRIFACEMESH3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbTriFaceMesh3D.h>
#include <basics/utilities/UbFileInputASCII.h>

#ifdef CAB_QT 
#include <qfiledialog.h>    
#endif

#ifdef CAB_VTK
#include <numerics/geometry3d/presentation/vtkGbTriangularMesh3D.h>
#endif

class GbTriFaceMesh3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbTriFaceMesh3DCreator* getInstance()
   {
      static GbTriFaceMesh3DCreator instance;
      return &instance;
   }
   static GbTriFaceMesh3D* readMeshFromFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);

   static GbTriFaceMesh3D* readMeshFromMeshFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);
   static GbTriFaceMesh3D* readMeshFromMeshFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);

   static GbTriFaceMesh3D* readMeshFromPLYFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);
   static GbTriFaceMesh3D* readMeshFromPLYFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);

   static GbTriFaceMesh3D* readMeshFromSTLFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true); 
   static GbTriFaceMesh3D* readMeshFromSTLFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);
   static GbTriFaceMesh3D* readMeshFromSTLFile2(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true,  bool isBinaryFormat=true);

   static GbTriFaceMesh3D* readMeshFromAVSFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true); 
   static GbTriFaceMesh3D* readMeshFromAVSFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);

   static GbTriFaceMesh3D* readMeshFromVTKASCIIFile(std::string filename, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true); 
   static GbTriFaceMesh3D* readMeshFromVTKASCIIFile(UbFileInput* in, std::string meshName="", GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg = GbTriFaceMesh3D::KDTREE_SAHPLIT, bool removeRedundantNodes=true);

   GbTriFaceMesh3D* createGbObject3D() { return new GbTriFaceMesh3D(); }
   
   std::string getGbObject3DTypeID(){ return "GbTriFaceMesh3D"; };
   std::string toString()           { return "GbTriFaceMesh3DCreator"; }

#ifdef CAB_QT 


   GbTriFaceMesh3D* createGbObject3DwithQt()
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
      GbTriFaceMesh3D *mesh = NULL;//GbTriFaceMesh3DCreator::readMeshFromSTLFile(&in, stream.str() );
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
   GbTriFaceMesh3DCreator( const GbTriFaceMesh3DCreator& );                  //no copy allowed 
   const GbTriFaceMesh3DCreator& operator=( const GbTriFaceMesh3DCreator& ); //no copy allowed
   GbTriFaceMesh3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbTriFaceMesh3DCreator::getInstance()), CAB_GbTriFaceMesh3DCreator);
#endif

#endif
