#include <basics/utilities/UbTiming.h>
#include <basics/utilities/UbRandom.h>
#include <basics/utilities/UbTuple.h>

#include <basics/writer/WbWriterAvsASCII.h>
#include <basics/writer/WbWriterAvsBinary.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/writer/WbWriterVtkXmlASCII.h>

#include <basics/container/CbUniformMatrix3D.h>

#include <numerics/geometry3d/GbTriFaceMesh3D.h>
#include <numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h>
#include <numerics/geometry3d/GbHalfSpace3D.h>

#include <numerics/geometry3d/KdTree/KdTree.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>

#include <stack>
#include <list>

using namespace std;

void createGrid(std::string meshfile, int maxNofPointsPerDir, Kd::SplitAlgorithm<float>& splitAlg, bool writeFiles = true, std::string outpath = "g:/temp");
void recursiveGridFill(CbUniformMatrix3D<short>& grid, const short& xs, const short& ys, const short& zs, const short& type);
void iterativeGridFill(CbUniformMatrix3D<short>& grid, const short& xs, const short& ys, const short& zs, const short& type);

#include <3rdParty/dirstream/dirstream.h>
#include <3rdParty/dirstream/filter_utils.h>  // enthält die Definiton der bool'schen Ops für Filter


using namespace std;

int main()
{
   try
   {
//       //////////////////////////////////////////////////////////////////////////
//       //STL File auswaehlen
//       //////////////////////////////////////////////////////////////////////////
//       string pathname = "c:/temp";
//       string stlPath = "c:/temp";
//       dirstr::dirstream str(stlPath.c_str(), dirstr::op(dirstr::pattern_f("*.stl")) 
//                            || dirstr::op(dirstr::pattern_f("*.ply"))  || dirstr::op(dirstr::pattern_f("*.inp"))
//                            || dirstr::op(dirstr::pattern_f("*.mesh")));
// 
// //      UbLog::setReportingLevel(logDEBUG5);
//       UbLog::setReportingLevel(logINFO);
// 
//       vector<string> filenames;
//       for(string entry; str >> entry;)
//       {
//          GbTriFaceMesh3D* mesh = GbTriFaceMesh3DCreator::getInstance()->readMeshFromFile(entry,"mesh");
// 
//          string fn = UbSystem::getFilenameFromString(entry);
//          mesh->writeMeshPly(pathname+"/"+fn+".ply");
// 
//          delete mesh;
//       }
// 
//       exit(0);

      //createGrid("c:/temp/clumps.stl", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/50spheres.inp", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/Sphere5040.inp", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/cooling_2.inp", 400, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/NDR-Konzertsaal.inp", 100, Kd::SAHSplit<float>() , true, "g:/temp");
      
      //createGrid("c:/temp/Campus-Details-W3.inp", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/Boein707.mesh", 100, Kd::SAHSplit<float>() , true, "g:/temp");
      
      //createGrid("c:/temp/dolphin.mesh", 400, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/box.ply", 10, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/bodyRight.stl", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/flamingo.mesh", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      
      //createGrid("c:/temp/torus.inp", 256, Kd::SAHSplit<float>() , true, "g:/temp");
      createGrid("c:/temp/xzx_dragon.stl", 512, Kd::SAHSplit<float>() , true, "g:/temp");
//       createGrid("c:/temp/bunny_ascii.ply", 256, Kd::SAHSplit<float>() , true, "g:/temp");
//       createGrid("c:/temp/dragon_ascii.ply", 256, Kd::SAHSplit<float>() , true, "g:/temp");
//       createGrid("c:/temp/budda_ascii.ply", 256, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGridWithLines("c:/temp/box.ply", 10, Kd::SAHSplit<float>() , true, "g:/temp");
      
//       createGrid("c:/temp/beatle.mesh",200, Kd::SAHSplit<float>(), true, "g:/temp");
//       createGrid("c:/temp/atrium-30000tri.inp",200, Kd::SAHSplit<float>(), true, "g:/temp");
//       createGrid("c:/temp/Buero.inp",200, Kd::SAHSplit<float>(), true, "g:/temp");
//       createGrid("c:/temp/office_space.inp",200, Kd::SAHSplit<float>(), true, "g:/temp");

      //createGrid("d:/meshes/50spheres.inp",200, Kd::SAHSplit<float>(), true, "d:/temp");

      //createGrid("c:/temp/torus.inp",200, Kd::SAHSplit<float>(), true, "g:/temp");
      //createGrid("c:/temp/bodyRight.stl", 200, Kd::SAHSplit<float>() , true, "g:/temp");

      //createGrid("c:/temp/jetta.stl", 200, GbTriFaceMesh3D::KDTREE_SPATIALSPLIT , true, "g:/temp");
      //createGrid("c:/temp/jetta.stl", 200, Kd::SAHSplit<float>(), true, "g:/temp");
      //createGrid("c:/temp/VW_body.ply", 200, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/kugel.stl", 50, Kd::SAHSplit<float>() , true, "g:/temp");
      //createGrid("c:/temp/ship-2.mesh", 100, Kd::SAHSplit<float>() , true, "g:/temp/schiff2");                                                   
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

namespace Flag
{
   const short UNDEF = 2;
   const short SOLID = 1;
   const short FLUID = 0;
}

//////////////////////////////////////////////////////////////////////
void createGrid(std::string meshfile, int maxNofPointsPerDir, Kd::SplitAlgorithm<float>& splitAlg, bool writeFiles, std::string outpath)
{
   UbLog::setReportingLevel(logDEBUG5);
   std::string filename = UbSystem::getFilenameFromString(meshfile);

   std::list< UbTuple<string, double> > timerVals;
   UbTimer timer;
   timer.start();
   GbTriFaceMesh3D* mesh = GbTriFaceMesh3DCreator::getInstance()->readMeshFromFile(meshfile,"mesh",GbTriFaceMesh3D::KDTREE_SAHPLIT);
   timerVals.push_back( UbTuple<string, double>("mesh", timer.stop() ) );
   UBLOG( logINFO, "read mesh in "<<val<2>(timerVals.back())<<"s" );
   
   timer.start();
   Kd::Tree<float> kdTree( *mesh, splitAlg  );
   timerVals.push_back( UbTuple<string, double>("kdTree", timer.stop() ) );
   UBLOG( logINFO, "build tree in "<<val<2>(timerVals.back())<<"s" );

   UBLOG(logINFO, "############################################################");
   UBLOG(logINFO, "nodes of TriFaceMesh....... "<<mesh->getNodes()->size()      );
   UBLOG(logINFO, "triFaces of TriFaceMesh.... "<<mesh->getTriangles()->size()  );
   UBLOG(logINFO, "triFace copies in KdTree... "<<kdTree.getNumOfTriFaces()     );
   UBLOG(logINFO, "nodes of kdNodes of KdTree. "<<kdTree.getNumOfNodes()        );
   UBLOG(logINFO, "");

   //////////////////////////////////////////////////////////////////////////
   // Ausgangs 3-D_Feld erstellen
   //////////////////////////////////////////////////////////////////////////
   const float percentOverLap = 0.05f; //=5%
   const float maxLength = (1.0f+percentOverLap)*(float)UbMath::max( mesh->getLengthX1(), mesh->getLengthX2(), mesh->getLengthX3() );
   const float dx1 = maxLength/(maxNofPointsPerDir-1);
   const float dx2 = dx1;
   const float dx3 = dx1;

   const int nx1 = 1 + int( std::ceil(mesh->getLengthX1()*(1.0f+percentOverLap)/dx1)+0.5 );
   const int nx2 = 1 + int( std::ceil(mesh->getLengthX2()*(1.0f+percentOverLap)/dx2)+0.5 );
   const int nx3 = 1 + int( std::ceil(mesh->getLengthX3()*(1.0f+percentOverLap)/dx3)+0.5 );

   CbUniformMatrix3D<short> solids(nx1,nx2,nx3,Flag::UNDEF);

   //////////////////////////////////////////////////////////////////////////
   // Knoten typisieren
   //////////////////////////////////////////////////////////////////////////
   const float orgx1 = (float)( -0.5*percentOverLap*mesh->getLengthX1()+mesh->getX1Minimum() );
   const float orgx2 = (float)( -0.5*percentOverLap*mesh->getLengthX2()+mesh->getX2Minimum() );
   const float orgx3 = (float)( -0.5*percentOverLap*mesh->getLengthX3()+mesh->getX3Minimum() );

   std::vector<GbTriFaceMesh3D::TriFace>& triFaces = *mesh->getTriangles();
   std::vector<GbTriFaceMesh3D::Vertex>&  vertices = *mesh->getNodes();

   float x1center = float( mesh->getX1Centroid() );
   float x2center = float( mesh->getX2Centroid() );
   float x3center = float( mesh->getX3Centroid() );

   UBLOG( logINFO, "performing point-in-object(PIO)-tests");
   long  counter1=0, counter2 = 0;
   float x1w, x2w, x3w;
   int   x1Min, x2Min, x3Min, x1Max, x2Max, x3Max;
   float einflussBereichKnoten_sq = dx1*dx1+dx2*dx2+dx3*dx3;

   timer.start();
   for(size_t t=0; t<triFaces.size(); t++)   
   {
      GbTriFaceMesh3D::TriFace& triangle = triFaces[t];
      GbTriFaceMesh3D::Vertex& v1 = vertices[triangle.v1];
      GbTriFaceMesh3D::Vertex& v2 = vertices[triangle.v2];
      GbTriFaceMesh3D::Vertex& v3 = vertices[triangle.v3];

      //////////////////////////////////////////////////////////////////////////
      // AABB riangle
      //////////////////////////////////////////////////////////////////////////
      x1Min = /*UbMath::integerRounding*/( std::floor( ( UbMath::min( v1.x, v2.x, v3.x ) - orgx1) / dx1 ) );
      x2Min = /*UbMath::integerRounding*/( std::floor( ( UbMath::min( v1.y, v2.y, v3.y ) - orgx2) / dx2 ) );
      x3Min = /*UbMath::integerRounding*/( std::floor( ( UbMath::min( v1.z, v2.z, v3.z ) - orgx3) / dx3 ) );

      x1Max = /*UbMath::integerRounding*/( std::ceil(  ( UbMath::max( v1.x, v2.x, v3.x ) - orgx1) / dx1 ) );
      x2Max = /*UbMath::integerRounding*/( std::ceil(  ( UbMath::max( v1.y, v2.y, v3.y ) - orgx2) / dx2 ) );
      x3Max = /*UbMath::integerRounding*/( std::ceil(  ( UbMath::max( v1.z, v2.z, v3.z ) - orgx3) / dx3 ) );

      GbHalfSpace3D halfSpace(  v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z );

      for(int x3=x3Min; x3<=x3Max; x3++)
      {
         for(int x2=x2Min; x2<=x2Max; x2++)
         {
            for(int x1=x1Min; x1<=x1Max; x1++)
            {
               counter1++;
               
               short& solidVal = solids(x1,x2,x3);

               if( solidVal != Flag::UNDEF )  //doppeltes Testen vermeiden
               {
                  continue;
               }
 
               counter2++;
               
               //Weltkoords
               x1w = orgx1+x1*dx1;
               x2w = orgx2+x2*dx2;
               x3w = orgx3+x3*dx3;

               float dist = (float)halfSpace.getDistance( x1w, x2w, x3w );
               if( UbMath::greater( dist, 0.0f) )
               {
                  continue;
               }
               if( UbMath::greater(dist*dist, einflussBereichKnoten_sq))
               {
                  continue;
               }
               
               //eigentlicher PIO-Test
               bool testFailed = true;
               for(int i=0; i<100; i++ )
               {
                  Kd::Ray<float> ray(  x1w, x2w, x3w  //, 1, 0 ,0 );
                                     , ( x1w < x1center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                     , ( x2w < x2center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                     , ( x3w < x3center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) ) );

                  int iSec = kdTree.intersectRay( ray, Kd::CountRayIntersectionHandler<float>() );
                  
                  if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
                  {
                     if(iSec == Kd::Intersection::ON_BOUNDARY )
                     {
                        solidVal = Flag::SOLID;
                     }
                     else
                     {
                        solidVal = (iSec&1);  //ungerade anzahl an schnitten --> drinnen
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
             }
         }
      }
   }
   timerVals.push_back( UbTuple<string, double>("PiO-Test", timer.stop() ) );
   UBLOG( logINFO,counter2 <<" point-in-object(PIO)-tests done in "<<val<2>(timerVals.back())<<"s" );
   UBLOG( logINFO,counter1-counter2 <<" point-in-object(PIO)-tests uebersprungen" );

   //////////////////////////////////////////////////////////////////////////
   // FLOOD FILL
   //////////////////////////////////////////////////////////////////////////

   if( false) //using just one seed point
   {
      //FUELL
      bool foundSeedPoint         = false;
      int  seedPointSearchCounter = 0;
      int  seedX1 = Ub::inf;
      int  seedX2 = Ub::inf;
      int  seedX3 = Ub::inf;

      timer.start();
      for(size_t t=0; t<triFaces.size(); t++)   
      {
          GbTriFaceMesh3D::TriFace& triangle = triFaces[t];
          
          float& nx = triangle.nx;
          float& ny = triangle.ny;
          float& nz = triangle.nz;

          float cx1 = triangle.getX1Centroid(vertices);
          float cx2 = triangle.getX2Centroid(vertices);
          float cx3 = triangle.getX3Centroid(vertices);

          for(int k=0; k<5; k++) 
          {
             seedPointSearchCounter++;

             cx1 -= nx * dx1;
             cx2 -= ny * dx2;
             cx3 -= nz * dx3;

             int ix1 = UbMath::integerRounding( (cx1-orgx1)/dx1 );
             int ix2 = UbMath::integerRounding( (cx2-orgx2)/dx2 );
             int ix3 = UbMath::integerRounding( (cx3-orgx3)/dx3 );

             if(   solids.indicesInRange(ix1,ix2,ix3)
                && solids(ix1, ix2, ix3 ) == Flag::UNDEF )
             {
                x1w = orgx1+ix1*dx1;
                x2w = orgx2+ix2*dx2;
                x3w = orgx3+ix3*dx3;

                Kd::Ray<float> ray(  x1w, x2w, x3w  //, 1, 0 ,0 );
                                  , ( x1w < x1center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                  , ( x2w < x2center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                  , ( x3w < x3center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) ) );

                int iSec = kdTree.intersectRay( ray, Kd::CountRayIntersectionHandler<float>() );

                if( iSec>0 && (iSec&1) )
                {
                   seedX1 = ix1;
                   seedX2 = ix2;
                   seedX3 = ix3;
                   foundSeedPoint = true;
                   break;
                }
              }
          }
          if(foundSeedPoint) break;
      }
      if(!foundSeedPoint)
         throw UbException(UB_EXARGS,"fuck no seed point found");
      timerVals.push_back( UbTuple<string, double>("Seed found in", timer.stop() ) );
      UBLOG( logINFO,"found seed Point in "<<val<2>(timerVals.back())<<"s with "<<seedPointSearchCounter << " tested points" );

      cout<<nx1<<","<<nx2<<","<<nx3<<endl;
      bool recursiveFloodFill = ( nx1*nx2*nx3 < 100*100*30 );
      if(recursiveFloodFill)
      {
         timer.start();
         recursiveGridFill(solids, seedX1, seedX2, seedX3, Flag::SOLID);
         timerVals.push_back( UbTuple<string, double>("flood fill (r)", timer.stop() ) );
         UBLOG( logINFO,"recursive flood fill in "<<val<2>(timerVals.back())<<"s with "<<seedPointSearchCounter << " tested points" );

         CbUniformMatrix3D<short> solidsCpy(solids);
         timer.start();
         iterativeGridFill(solidsCpy, seedX1, seedX2, seedX3, Flag::SOLID);
         timerVals.push_back( UbTuple<string, double>("flood fill (i)", timer.stop() ) );
         UBLOG( logINFO,"iterative flood fill in "<<val<2>(timerVals.back())<<"s with "<<seedPointSearchCounter << " tested points" );
      }
      else
      {
         timer.start();
         iterativeGridFill(solids, seedX1, seedX2, seedX3, Flag::SOLID);
         timerVals.push_back( UbTuple<string, double>("flood fill (r)", timer.stop() ) );
         UBLOG( logINFO,"recursive flood fill in "<<val<2>(timerVals.back())<<"s with "<<seedPointSearchCounter << " tested points" );
      }
      
      UBLOG(logINFO, "############################################################");

   }
   else //verifying complete arry
   {
      bool recursiveFloodFill = ( nx1*nx2*nx3 < 100*100*30 );
      int solidCounter  = 0;

      timer.start();
      for(int x3=0; x3<solids.getNX3(); x3++)
         for(int x2=0; x2<solids.getNX2(); x2++)
            for(int x1=0; x1<solids.getNX1(); x1++)
            {
               if( solids(x1  ,x2  ,x3 ) == Flag::UNDEF ) 
               {
                  x1w = orgx1+x1*dx1;
                  x2w = orgx2+x2*dx2;
                  x3w = orgx3+x3*dx3;

                  int iSec = -1;
                  do{
                     Kd::Ray<float> ray(  x1w, x2w, x3w  //, 1, 0 ,0 );
                                       , ( x1w < x1center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                       , ( x2w < x2center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) )
                                       , ( x3w < x3center ? (float)UbRandom::rand(-1.0,-0.001, 10) : (float)UbRandom::rand(0.001, 1.0, 10) ) );

                     iSec = kdTree.intersectRay( ray, Kd::CountRayIntersectionHandler<float>() );
                  }while(iSec<0);

                  if( iSec&1 )
                  {
                     if(recursiveFloodFill) recursiveGridFill(solids,x1,x2,x3,Flag::SOLID);
                     else                   iterativeGridFill(solids,x1,x2,x3,Flag::SOLID);
                  }
                  else
                  {
                     if(recursiveFloodFill) recursiveGridFill(solids,x1,x2,x3,Flag::FLUID);
                     else                   iterativeGridFill(solids,x1,x2,x3,Flag::FLUID);
                  }
               }
            }
      if(recursiveFloodFill) timerVals.push_back( UbTuple<string, double>("flood fill (r)", timer.stop() ) );
      else                   timerVals.push_back( UbTuple<string, double>("flood fill (i)", timer.stop() ) );
      UBLOG( logINFO,"recursive flood fill in "<<val<2>(timerVals.back())<<"s " );
   }

   list< UbTuple< string, double > >::iterator iter;
   for(iter = timerVals.begin(); iter!=timerVals.end(); ++iter)
   {
      UBLOG( logINFO, setw(16) << val<1>(*iter) << " in " << setw(8) << setprecision(8) << val<2>(*iter) << "s" );
   }

   int solidCounter  = 0;
   for(int x3=0; x3<solids.getNX3(); x3++)
      for(int x2=0; x2<solids.getNX2(); x2++)
         for(int x1=0; x1<solids.getNX1(); x1++)
         {
            if( solids(x1  ,x2  ,x3 ) == Flag::SOLID ) 
            {
               solidCounter++;
            }
         }

   UBLOG( logINFO, "SOLIDS = " <<solidCounter);
   UBLOG( logINFO, "SOLIDS = " <<solidCounter);
   UBLOG( logINFO, "SOLIDS = " <<solidCounter);
   UBLOG( logINFO, "SOLIDS = " <<solidCounter);

   /* ======================================================================================= */
   if(writeFiles) 
   {
      UBLOG( logINFO, "writeFiles - start");
      string subfiledir = outpath+"/"+filename+"_solid_node_files";
      UbSystem::makeDirectory( subfiledir );

      std::vector<UbTupleFloat3 > nodes;
      std::vector<std::string > datanames(1,"data");
      datanames[0] = "solid";

      std::vector< std::string > outFilenames;

      std::vector<std::vector<double > > nodedata( datanames.size() );

      for(int x3=0; x3<solids.getNX3(); x3++)
         for(int x2=0; x2<solids.getNX2(); x2++)
            for(int x1=0; x1<solids.getNX1(); x1++)
            {
               if( solids(x1  ,x2  ,x3 ) == Flag::SOLID ) 
               {
                  nodes.push_back( makeUbTuple(  orgx1+x1*dx1, orgx2+x2*dx2, orgx3+x3*dx3 ) );
                  nodedata[0].push_back( solids(x1  ,x2  ,x3 ) );
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
      kdTree.writeTree(outpath+"/"+filename+"_kdTree",WbWriterVtkXmlBinary::getInstance());

      UBLOG( logINFO, "writeFiles - end")
   }

   delete mesh;
}

namespace Dirs
{
   const int X1[] = { 1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  0 };
   const int X2[] = { 0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,  0 };
   const int X3[] = { 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  0 };

   const int START = 0;
   const int END6  = 5;
   const int END18 = 17;
}
/*==================================================================*/
bool floodFillCheck(CbUniformMatrix3D<short>& grid, const short& x, const short& y, const short& z)
{
   return grid.indicesInRange( x, y, z ) && grid(x,y,z)==Flag::UNDEF;
}
int g_counter = 0;
void recursiveGridFill(CbUniformMatrix3D<short>& grid, const short& xs, const short& ys, const short& zs, const short& type)
{
   // Algorithmus zum Füllen eines Polyeders, ausgehend vom Saatpunkt xs,ys,zs

   //Saatknoten einfärben
   short& val = grid(xs,ys,zs);
   if( val==Flag::UNDEF )
   {
      val = type;
   }
   if(   floodFillCheck( grid, xs+1, ys  , zs   ) ) recursiveGridFill( grid, xs+1, ys  , zs  , type );
   if(   floodFillCheck( grid, xs  , ys+1, zs   ) ) recursiveGridFill( grid, xs  , ys+1, zs  , type );
   if(   floodFillCheck( grid, xs  , ys  , zs+1 ) ) recursiveGridFill( grid, xs  , ys  , zs+1, type );
   if(   floodFillCheck( grid, xs-1, ys  , zs   ) ) recursiveGridFill( grid, xs-1, ys  , zs  , type );
   if(   floodFillCheck( grid, xs  , ys-1, zs   ) ) recursiveGridFill( grid, xs  , ys-1, zs  , type );
   if(   floodFillCheck( grid, xs  , ys  , zs-1 ) ) recursiveGridFill( grid, xs  , ys  , zs-1, type );
}
/*==================================================================*/
void iterativeGridFill(CbUniformMatrix3D<short>& grid, const short& xs, const short& ys, const short& zs, const short& type)
{
   std::stack< UbTupleInt3 > stck;
   stck.push( UbTupleInt3(xs,ys,zs) );

   int x,y,z;

   while( !stck.empty() )  
   {
      x = val<1>( stck.top() );
      y = val<2>( stck.top() );
      z = val<3>( stck.top() );
      stck.pop();

      short& flagType = grid( x, y, z );

      if( flagType == Flag::UNDEF ) 
      {     
         flagType = type;

         if ( grid.indicesInRange( x+1, y  , z   ) ) stck.push( UbTupleInt3( x+1, y  , z   ) );
         if ( grid.indicesInRange( x  , y+1, z   ) ) stck.push( UbTupleInt3( x  , y+1, z   ) );
         if ( grid.indicesInRange( x  , y  , z+1 ) ) stck.push( UbTupleInt3( x  , y  , z+1 ) );
         if ( grid.indicesInRange( x-1, y  , z   ) ) stck.push( UbTupleInt3( x-1, y  , z   ) );
         if ( grid.indicesInRange( x  , y-1, z   ) ) stck.push( UbTupleInt3( x  , y-1, z   ) );
         if ( grid.indicesInRange( x  , y  , z-1 ) ) stck.push( UbTupleInt3( x  , y  , z-1 ) );
      }
   }
}
