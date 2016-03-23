#include <numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbTiming.h>

using namespace std;

GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
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
   transform(ext.begin(), ext.end(), ext.begin(), (int(*)(int))tolower); //(int(*)(int)) ist irgendso ein fieser cast, weil tolower ne alte c-methode ist

   //UBLOG(logINFO, "GbTriFaceMesh3DCreator::readMeshFromFile - read " <<filename );

   if     ( !ext.compare("ply" ) ) return GbTriFaceMesh3DCreator::readMeshFromPLYFile(filename, meshName,splitAlg , removeRedundantNodes);
   else if( !ext.compare("stl" ) ) return GbTriFaceMesh3DCreator::readMeshFromSTLFile(filename, meshName,splitAlg , removeRedundantNodes);
   else if( !ext.compare("inp" ) ) return GbTriFaceMesh3DCreator::readMeshFromAVSFile(filename, meshName,splitAlg , removeRedundantNodes);
   else if( !ext.compare("mesh") ) return GbTriFaceMesh3DCreator::readMeshFromMeshFile(filename, meshName,splitAlg , removeRedundantNodes);
   //else if( !ext.compare("raw") ) return GbTriFaceMesh3DCreator::readMeshFromRAWFile(filename,splitAlg , meshName);
   else throw UbException(UB_EXARGS,"fileformat "+ext);

   return NULL;
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromMeshFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UbFileInputASCII meshfile(filename);
   if(!meshfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriFaceMesh3DCreator::readMeshFromMeshFile(&meshfile,meshName,splitAlg,removeRedundantNodes);
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromMeshFile(UbFileInput* in, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   vector<GbTriFaceMesh3D::Vertex>    *nodes     = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *triangles = new vector<GbTriFaceMesh3D::TriFace>;

   while( !in->eof() )
   {
      std::string line = in->readStringLine();
      if( line.find("Vertices") !=std::string::npos )
         break;
   }
   int numVertices = in->readInteger();

   UBLOG(logDEBUG1,"Number of vertices "<<numVertices);

   nodes->resize(numVertices);
   
   float x, y, z;
   for (int i=0; i<numVertices; i++)
   {
      x = in->readFloat();
      y = in->readFloat();
      z = in->readFloat();
      in->readLine();
      (*nodes)[i] = GbTriFaceMesh3D::Vertex(x,y,z);
   }
   UBLOG(logDEBUG1," - read vertices (#"<<numVertices<<") done");

   while( !in->eof() )
   {
      std::string line = in->readStringLine();
      if( line.find("Triangles") !=std::string::npos )
         break;
   }
   int numFaces = in->readInteger();
   triangles->reserve(numFaces);
   UBLOG(logDEBUG1,"Number of faces    "<<numFaces);

   int j,k,l;
   for(int i=0; i<numFaces; i++)
   {
      j = in->readInteger()-1;
      k = in->readInteger()-1;
      l = in->readInteger()-1;
      in->readLine();

      (*triangles).push_back(GbTriFaceMesh3D::TriFace(j,k,l));
   }
   UBLOG(logDEBUG1," - read faces (#"<<(int)triangles->size()<<", #nonTriangles="<<(int)triangles->size()-numFaces<<") done");

   GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D(meshName, nodes, triangles, splitAlg, removeRedundantNodes );

   return mesh;
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromPLYFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UbFileInputASCII plyfile(filename);
   if(!plyfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriFaceMesh3DCreator::readMeshFromPLYFile(&plyfile,meshName,splitAlg,removeRedundantNodes);
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromPLYFile(UbFileInput* in, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   //cout<<"GbTriangularMesh3DFile.readMeshFromPLYFile !!! Dieses Format hat leider redundante Knoten ..."<<endl;
   vector<GbTriFaceMesh3D::Vertex>    *nodes     = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *triangles = new vector<GbTriFaceMesh3D::TriFace>;
   
   float x, y, z;
   string dummy;

   int numVertices = in->readIntegerAfterString("element vertex");
   int numFaces    = in->readIntegerAfterString("element face");
   in->setPosAfterLineWithString("end_header");
   
   UBLOG(logDEBUG1,"Number of vertices "<<numVertices);
   UBLOG(logDEBUG1,"Number of faces    "<<numFaces);
   
   nodes->resize(numVertices);
   triangles->reserve(numFaces);

   int onePercent = (int)UbMath::max(1,UbMath::integerRounding(numVertices*0.01));
   for (int i=0; i<numVertices; i++)
   {
      if( i%onePercent==0 )
         cout<<" - read vertices (#"<<numVertices<<") "<<UbMath::integerRounding(i/(double)numVertices*100.0)<<"% "<<"\r"<<flush;
      x = in->readFloat();
      y = in->readFloat();
      z = in->readFloat();
      in->readLine();
      (*nodes)[i] = GbTriFaceMesh3D::Vertex(x,y,z);
   }
   UBLOG(logDEBUG1," - read vertices (#"<<numVertices<<") done");

   int p,j,k,l,n;
   onePercent = (int)UbMath::max(1,UbMath::integerRounding(numFaces*0.01));
   for(int i=0; i<numFaces; i++)
   {
      if( i%onePercent==0 ) cout<<" - read faces (#"<<numFaces<<") "<<UbMath::integerRounding(i/(double)numFaces*100.0)<<"% "<<"\r"<<flush;
      
      p = in->readInteger();
      if(p==3)  //Dreieck, alles andere wird stumpf ingnoriert
      {
         j = in->readInteger();
         k = in->readInteger();
         l = in->readInteger();

         if(   !UbMath::inClosedInterval(j,0,numVertices-1) 
            || !UbMath::inClosedInterval(k,0,numVertices-1) 
            || !UbMath::inClosedInterval(l,0,numVertices-1) ) 
         {         
            throw UbException(UB_EXARGS,"dreiecksindex ist groesser als max Knotenindex oder kleiner 0");
         }
         triangles->push_back(GbTriFaceMesh3D::TriFace(j,k,l));
      }
      else if(p==4)  //Viereck --> wird zu zwei Dreiecken!
      {
         j = in->readInteger();
         k = in->readInteger();
         l = in->readInteger();
         n = in->readInteger();
         numFaces++;
         i++;

         if(   !UbMath::inClosedInterval(j,0,numVertices-1) 
            || !UbMath::inClosedInterval(k,0,numVertices-1) 
            || !UbMath::inClosedInterval(l,0,numVertices-1) 
            || !UbMath::inClosedInterval(n,0,numVertices-1) 
            ) 
         {         
            throw UbException(UB_EXARGS,"vierecksindex ist groesser als max Knotenindex oder kleiner 0");
         }
         triangles->push_back(GbTriFaceMesh3D::TriFace(j,k,l));
         triangles->push_back(GbTriFaceMesh3D::TriFace(l,n,j));
      }

      in->readLine();

   }
   UBLOG(logDEBUG1," - read faces (#"<<(int)triangles->size()<<", #nonTriangles="<<(int)triangles->size()-numFaces<<") done");

   GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D(meshName, nodes, triangles, splitAlg, removeRedundantNodes);
   
   return mesh;
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromSTLFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UbFileInputASCII stlfile(filename);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriFaceMesh3DCreator::readMeshFromSTLFile(&stlfile,meshName,splitAlg,removeRedundantNodes);
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromSTLFile(UbFileInput *in, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UBLOG(logDEBUG1,"GbTriFaceMesh3DCreator::readMeshFromSTLFile !!! Dieses Format hat leider redundante Knoten ...");

   vector<GbTriFaceMesh3D::Vertex>    *nodes     = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *triangles = new vector<GbTriFaceMesh3D::TriFace>;
   string dummy;

   double x, y, z;
   int nr=0;

   in->readLine();
   while(dummy!="endsolid")
   {
      in->readLine();
      in->readLine();
      dummy = in->readString();
      if(dummy!="vertex") throw UbException(UB_EXARGS,"no vertex format");
      x=in->readDouble();
      y=in->readDouble();
      z=in->readDouble();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      in->readLine();
      in->readString();
      x=in->readDouble();
      y=in->readDouble();
      z=in->readDouble();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      in->readLine();
      in->readString();
      x=in->readDouble();
      y=in->readDouble();
      z=in->readDouble();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      triangles->push_back(GbTriFaceMesh3D::TriFace(nr,nr+1,nr+2));
      in->readLine();
      in->readLine();
      in->readLine();
      dummy = in->readString();
      nr+=3;
      //std::cout<<"read mesh "<< nr <<" \n";
   }

   GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D(meshName, nodes, triangles, splitAlg, removeRedundantNodes);
   
   return mesh;
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromAVSFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UbFileInputASCII stlfile(filename);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriFaceMesh3DCreator::readMeshFromAVSFile(&stlfile,meshName,splitAlg,removeRedundantNodes);
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromAVSFile(UbFileInput *in, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg , bool removeRedundantNodes)
{
   UBLOG(logDEBUG1,"GbTriFaceMesh3DCreator.readMeshFromAVSFile !!! Dieses Format hat leider redundante Knoten ...");

   vector<GbTriFaceMesh3D::Vertex>    *nodes     = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *triangles = new vector<GbTriFaceMesh3D::TriFace>;
   string dummy;

   in->readLine();
   int numberNodes = in->readInteger();
   int numberTris  = in->readInteger();
   in->readLine();

   double x,y,z;
   for(int u=0;u<numberNodes;u++)
   {
      in->readInteger();
      x=in->readDouble();
      y=in->readDouble();
      z=in->readDouble();
      in->readLine();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
   }
   int id1,id2,id3;
   for(int u=0;u<numberTris;u++)
   {
      in->readInteger();
      in->readInteger();
      in->readString();
      id1 = in->readInteger();
      id2 = in->readInteger();
      id3 = in->readInteger();
      triangles->push_back(GbTriFaceMesh3D::TriFace(id1-1,id2-1,id3-1));
   }

   GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D(meshName, nodes, triangles, splitAlg, removeRedundantNodes);
   
   return mesh;
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromVTKASCIIFile(string filename, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UbFileInputASCII stlfile(filename);
   if(!stlfile) throw UbException(UB_EXARGS,"cannot open file "+filename);
   return GbTriFaceMesh3DCreator::readMeshFromVTKASCIIFile(&stlfile,meshName,splitAlg,removeRedundantNodes);
}
/*======================================================================*/
GbTriFaceMesh3D* GbTriFaceMesh3DCreator::readMeshFromVTKASCIIFile(UbFileInput *in, string meshName, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
{
   UBLOG(logDEBUG1,"GbTriFaceMesh3DCreator.readMeshFromVTKASCIIFile !!! Dieses Format hat leider redundante Knoten ...");

   vector<GbTriFaceMesh3D::Vertex>    *nodes     = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *triangles = new vector<GbTriFaceMesh3D::TriFace>;
   string dummy;

   in->readLine();
   in->readLine();
   in->readLine();
   in->readLine();
   
   in->readString();
   int numberNodes = in->readInteger();
   in->readLine();

   double x,y,z;
   for(int u=0;u<numberNodes;u++)
   {
      x=in->readDouble();
      y=in->readDouble();
      z=in->readDouble();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      //cout<<u<<" - x,y,z:"<<x<<","<<y<<","<<z<<endl;
      //x=in->readDouble();
      //y=in->readDouble();
      //z=in->readDouble();
      //nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      //x=in->readDouble();
      //y=in->readDouble();
      //z=in->readDouble();
      //nodes->push_back(GbTriFaceMesh3D::Vertex((float)x,(float)y,(float)z));
      //in->readLine();
   }
   in->readLine();
   in->readString();
   int numberTris  = in->readInteger();
   in->readLine();
   UBLOG(logDEBUG1,"numberTris:"<<numberTris);

   int id1,id2,id3;
   for(int u=0;u<numberTris;u++)
   {
      in->readInteger();
      id1 = in->readInteger();
      id2 = in->readInteger();
      id3 = in->readInteger();
      triangles->push_back(GbTriFaceMesh3D::TriFace(id1,id2,id3));
      //cout<<u<<" - id1,id2,id3:"<<id1<<","<<id2<<","<<id3<<endl;
   }
   UBLOG(logDEBUG1,"Tris gelesen");

   GbTriFaceMesh3D* mesh = new GbTriFaceMesh3D(meshName, nodes, triangles, splitAlg, removeRedundantNodes);
   UBLOG(logDEBUG1,"mesh erzeugt (with remove redundant nodes = "<< boolalpha <<removeRedundantNodes<<")");


   return mesh;
}
