#include <basics/utilities/UbTiming.h>
#include <basics/utilities/UbRandom.h>

#include <basics/writer/WbWriterAvsASCII.h>
#include <basics/writer/WbWriterAvsBinary.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>

#include <basics/container/CbUniformMatrix4D.h>

#include <numerics/geometry3d/GbTriFaceMesh3D.h>
#include <numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h>

#include <numerics/geometry3d/KdTree/KdTree.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>

using namespace std;

void KdTreeTest         (std::string meshfile, int maxNofPointsPerDir, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM pio, bool writeFiles = true, std::string outpath = "g:/temp");
void KdTreeTestWithLines(std::string meshfile, int maxNofPointsPerDir, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM pio, bool writeFiles = true, std::string outpath = "g:/temp");

int main()
{
   try
   {
      //KdTreeTest("c:/temp/clumps.stl", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/50spheres.inp", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/Sphere5040.inp", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/cooling_2.inp", 400, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/NDR-Konzertsaal.inp", 100, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      
      //KdTreeTest("c:/temp/Campus-Details-W3.inp", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/Boein707.mesh", 100, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      
      //KdTreeTest("c:/temp/dolphin.mesh", 400, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/box.ply", 10, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/bodyRight.stl", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/flamingo.mesh", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      
      //KdTreeTest("c:/temp/torus.inp", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTestWithLines("c:/temp/box.ply", 10, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      
      KdTreeTest("c:/temp/doppelwandbox.ply", 100, GbTriFaceMesh3D::KDTREE_SPATIALSPLIT, true, "g:/temp");
      
      //KdTreeTestWithLines("c:/temp/torus.inp", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT, true, "g:/temp");

      //KdTreeTestWithLines("c:/temp/bodyRight.stl", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");

      //KdTreeTest("c:/temp/jetta.stl", 200, GbTriFaceMesh3D::KDTREE_SPATIALSPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/jetta.stl", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT, true, "g:/temp");
      //KdTreeTest("c:/temp/VW_body.ply", 200, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/kugel.stl", 50, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp");
      //KdTreeTest("c:/temp/ship-2.mesh", 100, GbTriFaceMesh3D::KDTREE_SAHPLIT , true, "g:/temp/schiff2");                                                   
   }
   catch(const std::exception& e)
   {
      UBLOG2(  logERROR, std::cerr, "Caught exception:");
      UBLOG2(  logERROR, std::cerr, "Type: " << typeid(e).name() );
      UBLOG2ML(logERROR, std::cerr, "What: " << e.what() );
   }
   catch(...)
   {
      UBLOG2(logERROR, std::cerr, "unknown exception occurs in "<< UB_FUNCTION)
   }



  
}

//////////////////////////////////////////////////////////////////////
void KdTreeTest(std::string meshfile, int maxNofPointsPerDir, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM pio, bool writeFiles, std::string outpath)
{
   UbLog::setReportingLevel(logDEBUG5);
   std::string filename = UbSystem::getFilenameFromString(meshfile);

   GbTriFaceMesh3D* mesh = GbTriFaceMesh3DCreator::getInstance()->readMeshFromFile(meshfile,"mesh",pio);
   mesh->scale(10000,10000,10000);
   //dummy test, damit der baum erstellt wird
   mesh->isPointInGbObject3D(0,0,0);

   UBLOG(logINFO, "############################################################");
   UBLOG(logINFO, "nodes of TriFaceMesh....... "<<mesh->getNodes()->size()              );
   UBLOG(logINFO, "triFaces of TriFaceMesh.... "<<mesh->getTriangles()->size()          );
   UBLOG(logINFO, "triFace copies in KdTree... "<<mesh->getKdTree()->getNumOfTriFaces() );
   UBLOG(logINFO, "nodes of kdNodes of KdTree. "<<mesh->getKdTree()->getNumOfNodes()   );
   UBLOG(logINFO, "");


   const float percentOverLap = 0.05f; //=5%
   const float maxLength = (1.0f+percentOverLap)*UbMath::max( mesh->getLengthX1(), mesh->getLengthX2(), mesh->getLengthX3() );
   const float dx1 = maxLength/(maxNofPointsPerDir-1);
   const float dx2 = dx1;
   const float dx3 = dx1;

   const int nx1 = 1 + int( std::ceil(mesh->getLengthX1()*(1.0f+percentOverLap)/dx1)+0.5 );
   const int nx2 = 1 + int( std::ceil(mesh->getLengthX2()*(1.0f+percentOverLap)/dx2)+0.5 );
   const int nx3 = 1 + int( std::ceil(mesh->getLengthX3()*(1.0f+percentOverLap)/dx3)+0.5 );

   CbUniformMatrix4D<int> solids(nx1,nx2,nx3,1,0);

   const float orgx1 = -0.5*percentOverLap*mesh->getLengthX1()+mesh->getX1Minimum();
   const float orgx2 = -0.5*percentOverLap*mesh->getLengthX2()+mesh->getX2Minimum();
   const float orgx3 = -0.5*percentOverLap*mesh->getLengthX3()+mesh->getX3Minimum();

   const float outX1 = 2*mesh->getX1Maximum();
   const float outX2 = 0;//2*mesh->getX2Maximum();
   const float outX3 = 0;//2*mesh->getX3Maximum();

   UBLOG( logINFO, "performing " << nx1*nx2*nx3  <<" point-in-object(PIO)-tests");
   UbTimer ff;
   ff.start();
   for(int x3=0; x3<solids.getNX3(); x3++)
      for(int x2=0; x2<solids.getNX2(); x2++)
         for(int x1=0; x1<solids.getNX1(); x1++)
         {
            solids(x1,x2,x3,0) = mesh->isPointInGbObject3D(orgx1+x1*dx1, orgx2+x2*dx2, orgx3+x3*dx3);
         }
   UBLOG( logINFO, nx1*nx2*nx3 <<" point-in-object(PIO)-tests done in "<<ff.stop()<<"s" );
   UBLOG(logINFO, "############################################################");


   /* ======================================================================================= */
   if(writeFiles) 
   {
      UBLOG( logINFO, "writeFiles - start");
      string subfiledir = outpath+"/"+filename+"_solid_node_files";
      UbSystem::makeDirectory( subfiledir );

      std::vector<UbTupleFloat3 > nodes;
      std::vector<std::string > datanames(solids.getNX4(),"data");
      datanames[0] = "solid";
      //datanames[1] = "solid";

      std::vector< std::string > outFilenames;

      std::vector<std::vector<double > > nodedata( datanames.size() );
      for(int x3=0; x3<solids.getNX3(); x3++)
         for(int x2=0; x2<solids.getNX2(); x2++)
            for(int x1=0; x1<solids.getNX1(); x1++)
            {
               if( solids(x1  ,x2  ,x3  , 0)  ) 
               {
                  nodes.push_back( makeUbTuple(  orgx1+x1*dx1, orgx2+x2*dx2, orgx3+x3*dx3 ) );

                  for(int i=0; i<solids.getNX4(); i++)
                  {
                     nodedata[i].push_back( solids(x1  ,x2  ,x3  ,i) );
                  }
               }          

               if(    nodes.size() > 2000000  
                   || ( x1==(solids.getNX1()-1) && x2==(solids.getNX2()-1) && x3==(solids.getNX3()-1) ) ) 
               {
                  outFilenames.push_back( WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(subfiledir+"/"+filename+"_solid_nodes_"+"_part"+UbSystem::toString(outFilenames.size()+1),nodes,datanames,nodedata) );
                  nodes.clear();
                  nodedata.clear();
                  nodedata.resize( datanames.size() );
               }
            }
            
      WbWriterVtkXmlBinary::getInstance()->writeCollection(outpath+"/"+filename+"_solids_nodes",outFilenames,0,false);
      
      
      mesh->writeMesh(outpath+"/"+filename+"_mesh",WbWriterVtkXmlBinary::getInstance(),true);
      mesh->writeMesh(outpath+"/"+filename+"_mesh",WbWriterAvsASCII::getInstance(),true);
      mesh->getKdTree()->writeTree(outpath+"/"+filename+"_kdTree",WbWriterVtkXmlBinary::getInstance());

      UBLOG( logINFO, "writeFiles - end")
   }
}

//////////////////////////////////////////////////////////////////////
void KdTreeTestWithLines(std::string meshfile, int maxNofPointsPerDir, GbTriFaceMesh3D::KDTREE_SPLITAGORITHM pio, bool writeFiles, std::string outpath)
{
   UbLog::setReportingLevel(logDEBUG5);
   std::string filename = UbSystem::getFilenameFromString(meshfile);

   GbTriFaceMesh3D* mesh = GbTriFaceMesh3DCreator::getInstance()->readMeshFromFile(meshfile,"mesh",pio);

   //dummy test, damit der baum erstellt wird
   mesh->isPointInGbObject3D(0,0,0);

   UBLOG(logINFO, "############################################################");
   UBLOG(logINFO, "nodes of TriFaceMesh....... "<<mesh->getNodes()->size()              );
   UBLOG(logINFO, "triFaces of TriFaceMesh.... "<<mesh->getTriangles()->size()          );
   UBLOG(logINFO, "triFace copies in KdTree... "<<mesh->getKdTree()->getNumOfTriFaces() );
   UBLOG(logINFO, "nodes of kdNodes of KdTree. "<<mesh->getKdTree()->getNumOfNodes()   );
   UBLOG(logINFO, "");


   const float percentOverLap = 0.05f; //=5%
   const float maxLength = (1.0f+percentOverLap)*UbMath::max( mesh->getLengthX1(), mesh->getLengthX2(), mesh->getLengthX3() );
   const float dx1 = maxLength/(maxNofPointsPerDir-1);
   const float dx2 = dx1;
   const float dx3 = dx1;

   const int nx1 = 1 + /*UbMath::integerRounding*/( std::ceil(mesh->getLengthX1()*(1.0f+percentOverLap)/dx1) );
   const int nx2 = 1 + /*UbMath::integerRounding*/( std::ceil(mesh->getLengthX2()*(1.0f+percentOverLap)/dx2) );
   const int nx3 = 1 + /*UbMath::integerRounding*/( std::ceil(mesh->getLengthX3()*(1.0f+percentOverLap)/dx3) );

   CbUniformMatrix4D<int> solids(nx1,nx2,nx3,2,0);

   const float orgx1 = -0.5*percentOverLap*mesh->getLengthX1()+mesh->getX1Minimum();
   const float orgx2 = -0.5*percentOverLap*mesh->getLengthX2()+mesh->getX2Minimum();
   const float orgx3 = -0.5*percentOverLap*mesh->getLengthX3()+mesh->getX3Minimum();

//    const float outX1 = 2*mesh->getX1Maximum();
//    const float outX2 = 2*mesh->getX2Maximum();
//    const float outX3 = 2*mesh->getX3Maximum();

   
   Kd::Tree<double>* kdTree = mesh->getKdTree();
   
   std::vector<GbTriFaceMesh3D::TriFace>& triFaces = *mesh->getTriangles();
   std::vector<GbTriFaceMesh3D::Vertex>&  vertices = *mesh->getNodes();

   float x1center = float( mesh->getX1Centroid() );
   float x2center = float( mesh->getX2Centroid() );
   float x3center = float( mesh->getX3Centroid() );

   UBLOG( logINFO, "performing point-in-object(PIO)-tests");
   UbTimer ff;
   ff.start();
   long counter1=0, counter2 = 0;
   for(size_t t=0; t<triFaces.size(); t++)   
   {
      int x1Min = /*UbMath::integerRounding*/( std::floor( (triFaces[t].getMinX(vertices)-orgx1) / dx1 ) );
      int x2Min = /*UbMath::integerRounding*/( std::floor( (triFaces[t].getMinY(vertices)-orgx2) / dx2 ) );
      int x3Min = /*UbMath::integerRounding*/( std::floor( (triFaces[t].getMinZ(vertices)-orgx3) / dx3 ) );

      int x1Max = /*UbMath::integerRounding*/( std::ceil(  (triFaces[t].getMaxX(vertices)-orgx1) / dx1 ) );
      int x2Max = /*UbMath::integerRounding*/( std::ceil(  (triFaces[t].getMaxY(vertices)-orgx2) / dx2 ) );
      int x3Max = /*UbMath::integerRounding*/( std::ceil(  (triFaces[t].getMaxZ(vertices)-orgx3) / dx3 ) );

      for(int x3=x3Min; x3<=x3Max; x3++)
         for(int x2=x2Min; x2<=x2Max; x2++)
            for(int x1=x1Min; x1<=x1Max; x1++)
            {
               counter1++;
               
               if( !solids.indicesInRange(x1,x2,x3,0) 
                  || solids(x1,x2,x3,1) == 1 )  //doppeltes Testeb vermeiden
               {
                  continue;
               }

               counter2++;
               
               double x1w = orgx1+x1*dx1;
               double x2w = orgx2+x2*dx2;
               double x3w = orgx3+x3*dx3;

               //eigentlicher PIO-Test
               bool testFailed = true;
               for(int i=0; i<100; i++)
               {
                  UbTupleDouble3 n1(x1w, x2w, x3w);
                  UbTupleDouble3 n2(  double( x1w < x1center ? mesh->getX1Minimum()-UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX1() : mesh->getX1Maximum()+UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX1() )
                                    , double( x2w < x2center ? mesh->getX2Minimum()-UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX2() : mesh->getX2Maximum()+UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX2() )
                                    , double( x3w < x3center ? mesh->getX3Minimum()-UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX3() : mesh->getX3Maximum()+UbRandom::rand(0.5, 1.0, 10)*mesh->getLengthX3() ) );

                  int iSec = kdTree->intersectLine( n1, n2, Kd::CountLineIntersectionHandler<double>() );
                  
                  if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
                  {
                     if(iSec == Kd::Intersection::ON_BOUNDARY )
                     {
                        solids(x1,x2,x3,0) = true;
                     }
                     else
                     {
                        solids(x1,x2,x3,0) = (iSec&1);  //ungerade anzahl an schnitten --> drinnen
                     }
                     testFailed = false;
                     break;
                  }
                  else
                  {
                     UBLOG(logDEBUG3, "GbTriFaceMesh3D.isPointInGbObject3D.if  - an edge was hit ");
                  }
               }
               if( testFailed ) throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
               solids(x1,x2,x3,1) = 1;
             }
   }
   UBLOG( logINFO,counter2 <<" point-in-object(PIO)-tests done in "<<ff.stop()<<"s" );
   UBLOG( logINFO,counter1-counter2 <<" point-in-object(PIO)-tests uebersprungen" );
   UBLOG(logINFO, "############################################################");

   /* ======================================================================================= */
   if(writeFiles) 
   {
      UBLOG( logINFO, "writeFiles - start");
      string subfiledir = outpath+"/"+filename+"_solid_node_files";
      UbSystem::makeDirectory( subfiledir );

      std::vector<UbTupleFloat3 > nodes;
      std::vector<std::string > datanames(solids.getNX4(),"data");
      datanames[0] = "solid";
      //datanames[1] = "solid";

      std::vector< std::string > outFilenames;

      std::vector<std::vector<double > > nodedata( datanames.size() );
      for(int x3=0; x3<solids.getNX3(); x3++)
         for(int x2=0; x2<solids.getNX2(); x2++)
            for(int x1=0; x1<solids.getNX1(); x1++)
            {
               if( solids(x1  ,x2  ,x3  , 0)  ) 
               {
                  nodes.push_back( makeUbTuple(  orgx1+x1*dx1, orgx2+x2*dx2, orgx3+x3*dx3 ) );

                  for(int i=0; i<solids.getNX4(); i++)
                  {
                     nodedata[i].push_back( solids(x1  ,x2  ,x3  ,i) );
                  }
               }          

               if(    nodes.size() > 2000000  
                   || ( x1==(solids.getNX1()-1) && x2==(solids.getNX2()-1) && x3==(solids.getNX3()-1) ) ) 
               {
                  outFilenames.push_back( WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(subfiledir+"/"+filename+"_solid_nodes_"+"_part"+UbSystem::toString(outFilenames.size()+1),nodes,datanames,nodedata) );
                  nodes.clear();
                  nodedata.clear();
                  nodedata.resize( datanames.size() );
               }
            }
   
      WbWriterVtkXmlBinary::getInstance()->writeCollection(outpath+"/"+filename+"_solids_nodes",outFilenames,0,false);
      
      
      mesh->writeMesh(outpath+"/"+filename+"_mesh",WbWriterVtkXmlBinary::getInstance());
      mesh->writeMesh(outpath+"/"+filename+"_mesh",WbWriterAvsASCII::getInstance());
      mesh->getKdTree()->writeTree(outpath+"/"+filename+"_kdTree",WbWriterVtkXmlBinary::getInstance());

      UBLOG( logINFO, "writeFiles - end")
   }
}
