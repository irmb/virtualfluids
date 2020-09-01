#include <geometry3d/fem/FeTriFaceMesh3D.h>

#include <geometry3d/fem/creator/FeTriFaceMesh3DCreator.h>
#include <geometry3d/GbTriangle3D.h>

#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>

using namespace std;

FeTriFaceMesh3D::FeTriFaceMesh3D():GbTriFaceMesh3D()
{
   this->attributes = new vector<VertexAttributes>;
//   this->createVertexTriFaceMap();
}
/*====================================================*/
FeTriFaceMesh3D::FeTriFaceMesh3D(std::string name, std::vector<Vertex>* nodes, std::vector<TriFace>* triangles):GbTriFaceMesh3D(name,nodes,triangles)
{
   this->attributes = new vector<VertexAttributes>;
   this->attributes->resize(nodes->size());
//   this->createVertexTriFaceMap();
}
/*======================================================================*/
ObObjectCreator* FeTriFaceMesh3D::getCreator()
{
   return FeTriFaceMesh3DCreator::getInstance();
}
/*====================================================*/
void FeTriFaceMesh3D::resizeAttributes()
{
   this->attributes->resize(nodes->size());
}
/*====================================================*/
//void FeTriFaceMesh3D::createVertexTriFaceMap()
//{
//   vertexTriFaceMap.clear();
//   vector<TriFace>& tris = *this->triangles;
//   vector<Vertex>&  pts  = *this->nodes;
//
//   for(size_t t=0; t<tris.size(); t++)
//   {
//      TriFace& tri = tris[t];
//      Vertex& vert1 = pts[tri.v1];
//      Vertex& vert2 = pts[tri.v2];
//      Vertex& vert3 = pts[tri.v3];
//      vertexTriFaceMap.setVertexTriFaceRelation(&vert1,&tri);
//      vertexTriFaceMap.setVertexTriFaceRelation(&vert2,&tri);
//      vertexTriFaceMap.setVertexTriFaceRelation(&vert3,&tri);
//   }
//}
/*====================================================*/
FeTriFaceMesh3D* FeTriFaceMesh3D::createMeshByTriangles(std::string name, std::vector<GbTriangle3D*> *triangles)
{
   vector<GbTriFaceMesh3D::Vertex>    *nodes = new vector<GbTriFaceMesh3D::Vertex>;
   vector<GbTriFaceMesh3D::TriFace>   *tris  = new vector<GbTriFaceMesh3D::TriFace>;
   int nr=0;
   for(int u=0;u<(int)triangles->size();u++)
   {
      if(UbMath::zero((*triangles)[u]->getArea())) continue;

      GbPoint3D* p1 = (*triangles)[u]->getPoint1();
      GbPoint3D* p2 = (*triangles)[u]->getPoint2();
      GbPoint3D* p3 = (*triangles)[u]->getPoint3();
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)p1->getX1Coordinate(),(float)p1->getX2Coordinate(),(float)p1->getX3Coordinate()));
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)p2->getX1Coordinate(),(float)p2->getX2Coordinate(),(float)p2->getX3Coordinate()));
      nodes->push_back(GbTriFaceMesh3D::Vertex((float)p3->getX1Coordinate(),(float)p3->getX2Coordinate(),(float)p3->getX3Coordinate()));
      tris->push_back(GbTriFaceMesh3D::TriFace(nr,nr+1,nr+2));
      nr+=3;
   }
   FeTriFaceMesh3D* triMesh = new FeTriFaceMesh3D(name, nodes, tris);
   triMesh->deleteRedundantNodes();
   triMesh->resizeAttributes();
   return triMesh;
}

/*======================================================================*/
UbTuple<string,string> FeTriFaceMesh3D::writeMesh(string filename, WbWriter* writer, bool writeNormals, std::vector< std::string >* datanames, std::vector< std::vector < double > >* nodedata)
{
   if(datanames || nodedata)
   {
      UBLOG(logWARNING,"FeTriFaceMesh3D::writeMesh - no support for extra datanames || nodedata");
   }

   UBLOG2(logDEBUG1,std::cout,"FeTriFaceMesh3D::writeMesh - start");

   UbTuple<string,string> filenames;

   if( dynamic_cast<WbWriterVtkXmlBinary*>(writer) || dynamic_cast<WbWriterVtkXmlASCII*>(writer))
   {
      vector< UbTupleFloat3 > triNodes( nodes->size() );
      vector< UbTupleInt3   > tris( triangles->size() );

      for(size_t i=0; i<nodes->size(); i++)
         triNodes[i] = makeUbTuple( (*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z );

      for(size_t i=0; i<triangles->size(); i++)
         tris[i] = makeUbTuple( (*triangles)[i].v1, (*triangles)[i].v2, (*triangles)[i].v3 ) ;

      vector<string> localDataNames;
      localDataNames.push_back("Fx"      );		
      localDataNames.push_back("Fy"      );		
      localDataNames.push_back("Fz"      );		
      localDataNames.push_back("sumFx"   );		
      localDataNames.push_back("sumFy"   );		
      localDataNames.push_back("sumFz"   );		
      localDataNames.push_back("vx"      );		
      localDataNames.push_back("vy"      );		
      localDataNames.push_back("vz"      );		
      localDataNames.push_back("S1"      );		
      localDataNames.push_back("S2"      );		
      localDataNames.push_back("S3"      );		
      localDataNames.push_back("S4"      );		
      localDataNames.push_back("S5"      );		
      localDataNames.push_back("S6"      );		
      localDataNames.push_back("Pressure");		

      std::vector< std::vector < double > > localNodedata( localDataNames.size(), std::vector<double>( nodes->size() ) );
      for(size_t n=0; n<nodes->size(); n++)
      {
         FeTriFaceMesh3D::VertexAttributes& attribut = (*this->attributes)[n];

         localNodedata[ 0][n] = attribut.getFX();
         localNodedata[ 1][n] = attribut.getFY();
         localNodedata[ 2][n] = attribut.getFZ();
         localNodedata[ 3][n] = attribut.getSumFX();
         localNodedata[ 4][n] = attribut.getSumFY();
         localNodedata[ 5][n] = attribut.getSumFZ();
         localNodedata[ 6][n] = attribut.getVelocityX();
         localNodedata[ 7][n] = attribut.getVelocityY();
         localNodedata[ 8][n] = attribut.getVelocityZ();
         localNodedata[ 9][n] = val<1>(attribut.getStresses());
         localNodedata[10][n] = val<2>(attribut.getStresses());
         localNodedata[11][n] = val<3>(attribut.getStresses());
         localNodedata[12][n] = val<4>(attribut.getStresses());
         localNodedata[13][n] = val<5>(attribut.getStresses());
         localNodedata[14][n] = val<6>(attribut.getStresses());
         localNodedata[15][n] = attribut.getPressure();
      }
      val<1>(filenames) = writer->writeTrianglesWithNodeData(filename,triNodes,tris,localDataNames,localNodedata);
   }
   else
   {
      string avsfilename = filename+writer->getFileExtension();
      val<1>(filenames)=avsfilename;

      ofstream out(avsfilename.c_str(),ios::out);
      if(!out)
      { 
         out.clear(); //flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
         string path = UbSystem::getPathFromString(filename);
         if(path.size()>0){UbSystem::makeDirectory(path);out.open(filename.c_str(),ios::out);}
         if(!out) throw UbException(UB_EXARGS,"file konnte nicht geschrieben werden");
      }

      //cout<<"AvsASCII - writeLines to "<<avsfilename<<" ...";

      int nofNodes = (int)nodes->size(); 
      int nofTrian = (int)triangles->size(); 
      int dataSize = 16;

      out<<"# UCD-File created by WbWriterAvsASCII\n";
      out<<nofNodes<<" "<<nofTrian<<" "<<dataSize<<" 0 0 "<<endl;

      for(int i=0; i<nofNodes; i++)
         out<<i+1<<" "<< (*nodes)[i].x<<" "<< (*nodes)[i].y<<" "<< (*nodes)[i].z<<" \n";

      for(int i=0; i<nofTrian; i++)
         out<<i+1<<" 2 tri "<<(*triangles)[i].v1+1<<" "<<(*triangles)[i].v2+1<<" "<<(*triangles)[i].v3+1<<" "<<endl;

      out<<dataSize;	
      for(int i=0;i<dataSize;i++) out<<" "<<1;
      out<<endl;

      out<<"Fx, no_unit"<<endl;		
      out<<"Fy, no_unit"<<endl;		
      out<<"Fz, no_unit"<<endl;		
      out<<"sumFx, no_unit"<<endl;		
      out<<"sumFy, no_unit"<<endl;		
      out<<"sumFz, no_unit"<<endl;		
      out<<"vx, no_unit"<<endl;		
      out<<"vy, no_unit"<<endl;		
      out<<"vz, no_unit"<<endl;		
      out<<"S1, no_unit"<<endl;		
      out<<"S2, no_unit"<<endl;		
      out<<"S3, no_unit"<<endl;		
      out<<"S4, no_unit"<<endl;		
      out<<"S5, no_unit"<<endl;		
      out<<"S6, no_unit"<<endl;		
      out<<"Pressure, no_unit"<<endl;		

      for(int n=0; n<nofNodes; n++)
      {
         FeTriFaceMesh3D::VertexAttributes& attribut = (*this->attributes)[n];

         double Fx = attribut.getFX();
         double Fy = attribut.getFY();
         double Fz = attribut.getFZ();
         double sumFx = attribut.getSumFX();
         double sumFy = attribut.getSumFY();
         double sumFz = attribut.getSumFZ();
         double vx = attribut.getVelocityX();
         double vy = attribut.getVelocityY();
         double vz = attribut.getVelocityZ();
         double p = attribut.getPressure();
         UbTupleDouble6& stresses = attribut.getStresses();
         out<<n+1<<" "<<Fx<<" "<<Fy<<" "<<Fz;
         out<<" "<<sumFx<<" "<<sumFy<<" "<<sumFz;
         out<<" "<<vx<<" "<<vy<<" "<<vz;
         out<<" "<<val<1>(stresses)<<" "<<val<2>(stresses)<<" "<<val<3>(stresses);
         out<<" "<<val<4>(stresses)<<" "<<val<5>(stresses)<<" "<<val<6>(stresses);
         out<<" "<<p<<endl;
      }
      out.close();
   }

   if(writeNormals)
   {
      vector<UbTupleFloat3 > lineNodes(triangles->size()*2);
      vector<UbTupleInt2 >   lines(triangles->size());
      for(size_t i=0; i<triangles->size(); i++)
      {
         TriFace& triangle = (*triangles)[i];
         lineNodes[i*2  ] = makeUbTuple( triangle.getX1Centroid(*nodes)
                                        ,triangle.getX2Centroid(*nodes)
                                        ,triangle.getX3Centroid(*nodes));
         lineNodes[i*2+1] = makeUbTuple( (float)(triangle.getX1Centroid(*nodes)+triangle.nx)
                                        ,(float)(triangle.getX2Centroid(*nodes)+triangle.ny)
                                        ,(float)(triangle.getX3Centroid(*nodes)+triangle.nz));

         lines[i] = makeUbTuple((int)i*2,(int)i*2+1);
      }
      val<2>(filenames) = writer->writeLines(filename+"_normals",lineNodes,lines);
   }


   if(writeNormals)
   {
      vector<UbTupleFloat3 > lineNodes(nodes->size()*2);
      vector<UbTupleInt2 >   lines(nodes->size());
      for(size_t i=0; i<nodes->size(); i++)
      {
   	    FeTriFaceMesh3D::VertexAttributes& attribut = (*this->attributes)[i];
         lineNodes[i*2  ] = makeUbTuple((*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z );
         lineNodes[i*2+1] = makeUbTuple((*nodes)[i].x+(float)attribut.getNormalX()
                                       ,(*nodes)[i].y+(float)attribut.getNormalY()
                                       ,(*nodes)[i].z+(float)attribut.getNormalZ());

         lines[i] = makeUbTuple((int)i*2,(int)i*2+1);
      }
      writer->writeLines(filename+"_PointNormals",lineNodes,lines);
   }

   UBLOG2(logDEBUG1,std::cout,"FeTriFaceMesh3D::writeMesh - end");

   return filenames;
}
