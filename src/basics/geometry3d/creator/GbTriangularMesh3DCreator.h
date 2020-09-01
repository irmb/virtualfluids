#ifndef GBTRIANGULARMESH3DCREATOR_H
#define GBTRIANGULARMESH3DCREATOR_H

#include <geometry3d/creator/GbObject3DCreator.h>
#include <geometry3d/GbTriangularMesh3D.h>
#include <basics/utilities/UbFileInputASCII.h>

#ifdef CAB_QT 
#include <qfiledialog.h>    
#endif

#ifdef CAB_VTK
#include <geometry3d/presentation/vtkGbTriangularMesh3D.h>
#endif

class GbTriangularMesh3DCreator : public GbObject3DCreator              
{                                       
public:
   static GbTriangularMesh3DCreator* getInstance()
   {
      static GbTriangularMesh3DCreator instance;
      return &instance;
   }
   
   static GbTriangularMesh3D* readMeshFromFile(std::string filename, std::string meshName="");

   static GbTriangularMesh3D* readMeshFromSTLFile(std::string filename, std::string meshName);
   static GbTriangularMesh3D* readMeshFromGTSFile(std::string filename, std::string meshName);     
   static GbTriangularMesh3D* readMeshFromPLYFile(std::string filename, std::string meshName);
   //static GbTriangularMesh3D* readMeshFromAVSFile(std::string filename, std::string meshName);
   static GbTriangularMesh3D* readMeshFromRAWFile(std::string inputFile, std::string meshName);

   static GbTriangularMesh3D* readMeshFromSTLFile(UbFileInput* infile, std::string meshName);
   static GbTriangularMesh3D* readMeshFromGTSFile(UbFileInput* infile, std::string meshName);     
   static GbTriangularMesh3D* readMeshFromPLYFile(UbFileInput* infile, std::string meshName);
   //static GbTriangularMesh3D* readMeshFromAVSFile(UbFileInput* infile, std::string meshName);
   static GbTriangularMesh3D* readMeshFromRAWFile(UbFileInput* infile, std::string meshName);
	
   
	GbTriangularMesh3D* createGbObject3D() { return new GbTriangularMesh3D(); }
   
   std::string getGbObject3DTypeID(){ return "GbTriangularMesh3D"; };
   std::string toString()           { return "GbTriangularMesh3DCreator"; }

#ifdef CAB_QT 

   GbTriangularMesh3D* createGbObject3DwithQt(QWidget* parent=0, Qt::WFlags flags=0)
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
      GbTriangularMesh3D *mesh = GbTriangularMesh3DCreator::readMeshFromSTLFile(&in, stream.str() );
      mesh->deleteRedundantNodes();
      return mesh;
   }
   //QDialog* getSpecificInstrument()  {  return 0;}
   void editGbObject3DwithQt(GbObject3D* gbObj, QWidget* parent=0, Qt::WFlags flags=0)
   { 
   }
#endif
#ifdef CAB_VTK
public:
   Presentator* createObjectPresentator(ObObject *object) { return new vtkGbTriangularMesh3D((GbTriangularMesh3D*)object); }
#endif


private:
   GbTriangularMesh3DCreator( const GbTriangularMesh3DCreator& );                  //no copy allowed 
   const GbTriangularMesh3DCreator& operator=( const GbTriangularMesh3DCreator& ); //no copy allowed
   GbTriangularMesh3DCreator() : GbObject3DCreator() {}
};

#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbTriangularMesh3DCreator::getInstance()), CAB_GbTriangularMesh3DCreator);
#endif

#endif
