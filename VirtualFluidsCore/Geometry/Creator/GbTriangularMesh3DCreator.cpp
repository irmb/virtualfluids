#include <numerics/geometry3d/creator/GbTriangularMesh3DCreator.h>
#include <algorithm>
#include <basics/utilities/UbLogger.h>

using namespace std;
                                               
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromFile(string filename, string meshName)
{
   if(meshName.empty())
   {
      size_t pos=filename.rfind("/");
      if(pos!=string::npos) meshName = filename.substr(pos+1);
      else                  meshName = filename; 
   }
   
   UbFileInputASCII stlfile(filename);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   string ext=stlfile.getFileExtension();

   //in "kleinbuchstaben" umwandeln
   transform(ext.begin(), ext.end(), ext.begin(), (int(*)(int))tolower);

   if     ( !ext.compare("ply") ) return GbTriangularMesh3DCreator::readMeshFromPLYFile(filename, meshName);
   else if( !ext.compare("stl") ) return GbTriangularMesh3DCreator::readMeshFromSTLFile(filename, meshName);
   else if( !ext.compare("gts") ) return GbTriangularMesh3DCreator::readMeshFromGTSFile(filename, meshName);
   else if( !ext.compare("raw") ) return GbTriangularMesh3DCreator::readMeshFromRAWFile(filename, meshName);
   else throw UbException(UB_EXARGS,"unrecognized fileformat "+ext);

   return NULL;
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromSTLFile(string filename, string meshName) 
{
   UbFileInputASCII stlfile(filename);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriangularMesh3DCreator::readMeshFromSTLFile(&stlfile,meshName);
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromSTLFile(UbFileInput* in, string meshName) 
{
   UBLOG(logINFO,"GbTriangularMesh3DFile.readMeshFromSTLFile !!! Dieses Format hat leider redundante Knoten ...");
   vector<GbPoint3D*>     *nodes     = new vector<GbPoint3D*>;
   vector<GbTriangle3D*>  *triangles = new vector<GbTriangle3D*>;
   nodes->resize(0, NULL);
   triangles->resize(0, NULL);
   double x, y, z;
   //int nr=0;
   string dummy;
   GbPoint3D     *node1    = NULL;                      
   GbPoint3D     *node2    = NULL;
   GbPoint3D     *node3    = NULL;
   GbTriangle3D *triangle = NULL;
   in->readLine();
   while(dummy!="endsolid")
   {		
      in->readLine();	
      in->readLine();	
      dummy = in->readString(); if(dummy!="vertex") throw UbException(UB_EXARGS,"no vertex format");
      x=in->readDouble(); 
      y=in->readDouble();           
      z=in->readDouble(); 
      node1 = new GbPoint3D(x,y,z); nodes->push_back(node1);
      in->readLine();	
      in->readString();	
      x=in->readDouble();
      y=in->readDouble(); 
      z=in->readDouble();	
      node2 = new GbPoint3D(x,y,z); nodes->push_back(node2);
      in->readLine();	
      in->readString();	
      x=in->readDouble();
      y=in->readDouble(); 
      z=in->readDouble();	
      node3 = new GbPoint3D(x,y,z); nodes->push_back(node3); 
      triangle = new GbTriangle3D(node1, node2, node3); triangles->push_back(triangle);
      in->readLine();
      in->readLine();
      in->readLine();
      dummy = in->readString();		
   }
   return new GbTriangularMesh3D(meshName, nodes, triangles);
}                                     
/*======================================================================*/
/**
* Returns a triangular mesh created from the specified TICAD source ASCII stream (system.dat format).
* @param in the input stream
* @param meshName the name of the created mesh
* @return a triangular mesh created from the specified TICAD source ASCII stream
* @exception IOException if any error occurs in creating the triangular mesh
*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromGTSFile(string inputfile, string meshName) 
{
   UbFileInputASCII gtlfile(inputfile);
   if(!gtlfile) throw UbException(UB_EXARGS,"cannot open file "+inputfile);
   return GbTriangularMesh3DCreator::readMeshFromGTSFile(&gtlfile,meshName);
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromGTSFile(UbFileInput *in, string meshName) 
{
   UBLOG(logINFO,"GbTriangularMesh3DFile.readMeshFromGTSFile !!! ");
   vector<GbPoint3D*>     *nodes     = new vector<GbPoint3D*>;
   vector<GbLine3D*>      *edges     = new vector<GbLine3D*>;
   vector<GbTriangle3D*>  *triangles = new vector<GbTriangle3D*>;
   nodes->resize(0, NULL);
   edges->resize(0, NULL);
   triangles->resize(0, NULL);
   double x, y, z;
   int point1, point2, point3;
   //int nr = 0;
   //in->readLine();
   int nodesize     = in->readInteger();
   int edgesize     = in->readInteger();
   int trianglesize = in->readInteger();
   UBLOG(logINFO,"node-/edge-/trianglesize: "<<nodesize<<" / "<<edgesize<<" / "<<trianglesize);

   for(int i=0; i<nodesize;i++)
   {		
      in->readLine();	
      x=in->readDouble(); 
      y=in->readDouble();  
      z=in->readDouble(); 
      nodes->push_back(new GbPoint3D(x,y,z));
   }
   for(int i=0; i<edgesize;i++)
   {		
      in->readLine();	
      point1=in->readInteger()-1; 
      point2=in->readInteger()-1; 
      edges->push_back(new GbLine3D((*nodes)[point1],(*nodes)[point2]));
   }
   for(int i=0; i<trianglesize;i++)
   {		
      in->readLine();	
      point1=in->readInteger();                
      point2=in->readInteger(); 
      point3=in->readInteger(); 
      //triangles->push_back(new GbTriangle3D((*nodes)[point1-1],(*nodes)[point2-1],(*nodes)[point3-1]));
      triangles->push_back(new GbTriangle3D((GbPoint3D*)(*edges)[point1-1]->getPoint1(),(GbPoint3D*)(*edges)[point2-1]->getPoint1(),(GbPoint3D*)(*edges)[point3-1]->getPoint1()));
   }
   return(new GbTriangularMesh3D(meshName, nodes, edges, triangles));
}                                  
/*======================================================================*/
/**
* Returns a triangular mesh created from the specified TICAD source ASCII stream (system.dat format).
* @param in the input stream
* @param meshName the name of the created mesh
* @return a triangular mesh created from the specified TICAD source ASCII stream
* @exception IOException if any error occurs in creating the triangular mesh
*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromPLYFile(string inputfile, string meshName) 
{
   UbFileInputASCII plyfile(inputfile);
   if(!plyfile) throw UbException(UB_EXARGS,"cannot open file "+inputfile);
   return GbTriangularMesh3DCreator::readMeshFromPLYFile(&plyfile,meshName);
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromPLYFile(UbFileInput *in, string meshName) 
{
   //cout<<"GbTriangularMesh3DFile.readMeshFromPLYFile !!! Dieses Format hat leider redundante Knoten ..."<<endl;
   vector<GbPoint3D*>     *nodes     = new vector<GbPoint3D*>;
   vector<GbTriangle3D*>  *triangles = new vector<GbTriangle3D*>;
   nodes->resize(0, NULL);
   triangles->resize(0, NULL);
   double x, y, z;
   int nr=0;
   string dummy;
   int numVertices, numFaces;
   GbPoint3D     *node     = NULL;
   GbPoint3D     *node1    = NULL;
   GbPoint3D     *node2    = NULL;
   GbPoint3D     *node3    = NULL;
   GbTriangle3D *triangle = NULL;
   in->readLine();
   in->readLine();
   in->readString(); in->readString(); numVertices = in->readInteger();
   in->readLine();
   in->readLine();
   in->readLine();
   in->readLine();
   in->readLine();
   in->readLine();
   in->readLine();
   in->readString(); in->readString(); numFaces = in->readInteger(); in->readLine();
   in->readLine();
   in->readLine();
   UBLOG(logINFO,"Number of vertices "<<numVertices);
   UBLOG(logINFO,"Number of faces    "<<numFaces);
   for (int i=0; i<numVertices; i++)
   {
      x = in->readDouble();
      y = in->readDouble();
      z = in->readDouble();
      //       cout<<x<<y<<z;
      //       in->readString(); in->readString(); in->readString();
      in->readLine();
      node = new GbPoint3D(x,y,z); nodes->push_back(node); 
   }
   nr=0;

   for (int i=0; i<numFaces; i++)
   {
      in->readString();
      int j,k,l;
      j = in->readInteger(); k = in->readInteger(); l = in->readInteger();
      node1 = (*nodes)[j];
      node2 = (*nodes)[k];
      node3 = (*nodes)[l];
      in->readLine();
      nr++;
      triangle = new GbTriangle3D(node1, node2, node3); triangles->push_back(triangle); 
   }

   return(new GbTriangularMesh3D(meshName, nodes, triangles));
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromRAWFile(string inputfile, string meshName) 
{
   UbFileInputASCII stlfile(inputfile);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+inputfile);
   return GbTriangularMesh3DCreator::readMeshFromRAWFile(&stlfile,meshName);
}
/*======================================================================*/
GbTriangularMesh3D* GbTriangularMesh3DCreator::readMeshFromRAWFile(UbFileInput *in, string meshName) 
{
   UBLOG(logINFO,"GbTriangularMesh3DFile.readMeshFromGTSFile !!! ");
   vector<GbPoint3D*>     *nodes     = new vector<GbPoint3D*>;
   vector<GbLine3D*>      *edges     = new vector<GbLine3D*>;
   vector<GbTriangle3D*>  *triangles = new vector<GbTriangle3D*>;
   nodes->resize(0, NULL);
   edges->resize(0, NULL);
   triangles->resize(0, NULL);
   double x, y, z;
   int point1, point2, point3;
   //int nr = 0;
   //in->readLine();
   int nodesize     = in->readInteger();
   int trianglesize = in->readInteger();
   int edgesize = 0;
   UBLOG(logINFO,"node-/edge-/trianglesize "<<nodesize<<" / "<<edgesize<<" / "<<trianglesize);

   for(int i=0; i<nodesize;i++)
   {		
      in->readLine();	
      x=in->readDouble(); 
      y=in->readDouble();  
      z=in->readDouble(); 
      nodes->push_back(new GbPoint3D(x,y,z));
   }
   for(int i=0; i<trianglesize;i++)
   {		
      in->readLine();	
      point1=in->readInteger();                
      point2=in->readInteger(); 
      point3=in->readInteger(); 
      triangles->push_back(new GbTriangle3D((*nodes)[point1],(*nodes)[point2],(*nodes)[point3]));
   }
   return(new GbTriangularMesh3D(meshName, nodes, edges, triangles));
}                                  
