#include <numerics/geometry3d/creator/GbQuadFaceMesh3DCreator.h>
#include <numerics/geometry3d/GbQuadFaceMesh3D.h>
#include <basics/utilities/UbFileInputASCII.h>

using namespace std;

/***************************************************************************/
GbObject3D* GbQuadFaceMesh3DCreator::createGbObject3D() 
{ 
   return new GbQuadFaceMesh3D(); 
}
/***************************************************************************/
GbQuadFaceMesh3D* GbQuadFaceMesh3DCreator::createQuadMesh3D(int nodesX1, int nodesX2, float startX1, float startX2, double knotenabstandX1, double knotenabstandX2, float nullNiveau, string name)
{
   vector<GbQuadFaceMesh3D::Vertex> *vertices = new vector<GbQuadFaceMesh3D::Vertex>;
   vector<GbQuadFaceMesh3D::QuadFace> *quads = new vector<GbQuadFaceMesh3D::QuadFace>;
   for(int x1=0;x1<nodesX1;x1++)
   {
      for(int x2=0;x2<nodesX2;x2++)
      {
         vertices->push_back(GbQuadFaceMesh3D::Vertex((float)(x1*knotenabstandX1+startX1), (float)(x2*knotenabstandX2+startX2), nullNiveau));
      }
   }
   for(int x1=0;x1<nodesX1-1;x1++)
   {
      for(int x2=0;x2<nodesX2-1;x2++)
      {
         int index = x1*nodesX2+x2;
         quads->push_back(GbQuadFaceMesh3D::QuadFace(index, index+nodesX2, index+nodesX2+1, index+1));
      }
   }
   
   return (new GbQuadFaceMesh3D(name, vertices, quads));
}

/*============================================================*/

#ifdef CAB_QT 

GbQuadFaceMesh3D* GbQuadFaceMesh3DCreator::createGbObject3DwithQt()
{
   //QString s = QFileDialog::getOpenFileName(NULL,NULL,NULL,"open file dialog","Choose a STL file" );
   //QString s = QFileDialog::getOpenFileName(NULL, "Choose a STL file", "/home", "*.stl");
   //QFileDialog* fd = new QFileDialog( NULL );
   //fd->setIconText(QString("Hallo"));
   //fd->show();
   //TODO: Open File Dialog einbauen.		
   //UbFileInputASCII in( s.toAscii().data() );
   //stringstream stream;
   //stream <<"TriangularMesh3D ";//<<_objCount++;
   //GbQuadFaceMesh3D *mesh = NULL;//GbQuadFaceMesh3DCreator::readMeshFromSTLFile(&in, stream.str() );
   //return mesh;
   return NULL;
}
//QDialog* getSpecificInstrument()  {  return 0;}
void GbQuadFaceMesh3DCreator::editGbObject3DwithQt(GbObject3D* gbObj)
{ 

}
#endif
#ifdef CAB_VTK
Presentator* GbQuadFaceMesh3DCreator::createObjectPresentator(ObObject *object) 
{
   return NULL;
}
#endif
