#include <iostream>
#include <string>

#include "geometry3d/CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "GenBlocksGridVisitor.h"
#include "geometry3d/GbSystem3D.h"
#include "geometry3d/GbCuboid3D.h"
#include "geometry3d/GbCylinder3D.h"
#include <geometry3d/GbSphere3D.h>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "RatioBlockVisitor.h"
#include "RatioSmoothBlockVisitor.h"
#include "OverlapBlockVisitor.h"
#include "RefineInterGbObjectsVisitor.h"
#include "RefineCrossAndInsideGbObjectBlockVisitor.h"
#include "SetKernelBlockVisitor.h"
#include "LBMKernelETD3Q27Cascaded.h"
#include "D3Q27MacroscopicQuantitiesPostprocessor.h"
#include "MPICommunicator.h"
#include "D3Q27ETBCProcessor.h"
#include "SimulationParameters.h"
#include "D3Q27SetUndefinedNodesBlockVisitor.h"
#include "SetInterpolationDirsBlockVisitor.h"
#include "D3Q27SetConnectorsBlockVisitor.h"
#include "NullCommunicator.h"
#include "D3Q27ETInitDistributionsBlockVisitor.h"
#include "CalculationManager.h"
#include "PQueuePartitioningGridVisitor.h"
#include "MetisPartitioningGridVisitor.h"
#include "D3Q27Interactor.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27BoundaryConditionAdapter.h"
#include "StringUtil.hpp"
#include "D3Q27OffsetInterpolationProcessor.h"
#include "D3Q27CompactInterpolationProcessor.h"
#include "SyncBcBlockVisitor.h"
#include "geometry3d/creator/GbTriFaceMesh3DCreator.h"
#include "geometry3d/GbTriFaceMesh3D.h"
#include "D3Q27TriFaceMeshInteractor.h"
#include "basics/utilities/UbFileOutputASCII.h"
#include "basics/utilities/UbFileInputASCII.h"
#include "basics/utilities/UbFileInputBinary.h"
#include "basics/container/CbArray3D.h"
#include "geometry3d/GbVoxelMatrix3D.h"


/* 
 *  The first 3 bits contain node type (fluid, inlet, etc.).
 *  The remaining 5 bits contain the unique geometry number 
 *  in case of solid nodes.
 *
 *  0 0 0 0 0 | 0 0 0
 */
#define getGeoType(geo) (geo & 0x7) /*= 00000111*/
#define getNormal(geo) (geo >> 3)
#define setGeoType(dest, geo_type) dest = (dest & 0xF8) + geo_type
#define setGeoNormal(dest, geo_id) dest = (geo_id << 3) + getGeoType(dest)

#define GEO_INVALID 0
#define GEO_FLUID 1
#define GEO_INLET 2
#define GEO_HULL 3           //hull
#define GEO_FLUID_IN_HULL 4  //fluid inside hull
#define GEO_BANANAS 5        //bananas    
#define GEO_BOX 6            //box

#define NORMAL_POS_X1 1
#define NORMAL_NEG_X1 2

#define NORMAL_POS_X2 3
#define NORMAL_NEG_X2 4

#define NORMAL_POS_X3 5
#define NORMAL_NEG_X3 6

#define CONVEXHULL

using namespace std;

typedef CbArray3D<int> VoxelMatrix;

//index             0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18
//f:                E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, ZERO
const int EX1[] = { 1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  0 };
const int EX2[] = { 0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,  0 };
const int EX3[] = { 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  0 };

//////////////////////////////////////////////////////////////////////////
void writeMatrixToVtkImageFile(const std::string& fileName, const VoxelMatrix& geoMatrix,
                               double itsDeltaXWorld, double orgX1, double orgX2, double orgX3)
{
   UbFileOutputASCII out(fileName);

   int NX1 = (int)geoMatrix.getNX1();	
   int NX2 = (int)geoMatrix.getNX2();	
   int NX3 = (int)geoMatrix.getNX3();
   int nn = NX1*NX2*NX3;
   out.writeLine("# vtk DataFile Version 3.0");
   out.writeLine(fileName);
   out.writeLine("ASCII");
   out.writeLine("DATASET STRUCTURED_POINTS");
   out.writeString("DIMENSIONS");
   out.writeInteger(NX1);
   out.writeInteger(NX2);
   out.writeInteger(NX3);
   out.writeLine();
   out.writeString("ORIGIN");
   out.writeDouble(orgX1);
   out.writeDouble(orgX2);
   out.writeDouble(orgX3);
   out.writeLine();
   out.writeString("SPACING");
   out.writeDouble(itsDeltaXWorld);
   out.writeDouble(itsDeltaXWorld);
   out.writeDouble(itsDeltaXWorld);
   out.writeLine();
   out.writeString("POINT_DATA");
   out.writeInteger(nn);
   out.writeLine();
   out.writeLine("SCALARS Geo integer");
   out.writeLine("LOOKUP_TABLE default");

   for(int k=0 ; k<NX3 ; k++){
      for(int j=0 ; j<NX2 ; j++){
         for(int i=0 ; i<NX1 ; i++){
            out.writeInteger( geoMatrix(i,j,k) );
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void readDimensionsFromFldFile(const std::string& fileName, int& d1, int& d2, int& d3)
{
   UbFileInputASCII in(fileName);
   // read grid nx3
   int dim   = in.readIntegerAfterString("ndim=");

   if (dim != 3) throw UbException(UB_EXARGS,"readGeoMatrixFromFldFile() - Wrong number of dimensions.");

   d1 = in.readIntegerAfterString("dim1=");
   d2 = in.readIntegerAfterString("dim2=");
   d3 = in.readIntegerAfterString("dim3=");
}
//////////////////////////////////////////////////////////////////////////
void readGeoMatrixFromFldFile(const std::string& fileName, GbVoxelMatrix3DPtr geoMatrix)
{
   UbFileInputASCII in(fileName);
   // read grid nx3
   int dim   = in.readIntegerAfterString("ndim=");

   if (dim != 3) throw UbException(UB_EXARGS,"readGeoMatrixFromFldFile() - Wrong number of dimensions.");

   int sizeX = in.readIntegerAfterString("dim1=");
   int sizeY = in.readIntegerAfterString("dim2=");
   int sizeZ = in.readIntegerAfterString("dim3=");

   std::string binFileName = in.readStringAfterString("variable 1 file=");

   //separate name from path
   std::string path = fileName.substr( 0, fileName.find_last_of('//')+1 );

   binFileName = path.append(binFileName);

   UbFileInputBinary binIn(binFileName);

   for (int i=0; i<2048; i++) 
   {
      binIn.readChar();
   }

   int x, y, z, val;

   for(z=0; z<sizeZ; z++)
   {
      for(y=0; y<sizeY; y++)
      {
         for(x=0; x<sizeX; x++)
         {
            val = binIn.readChar();

            if(x!=0 && x!=sizeX-1 && 
               y!=0 && y!=sizeY-1 &&
               z!=0 && z!=sizeZ-1   )
            {
               if(val == 0)
               {
                   (*geoMatrix)(x,y,z) = GbVoxelMatrix3D::SOLID;
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void discretizeGeoObject(GbObject3DPtr geoObject, VoxelMatrix& geoMatrix, double delta, 
                         double orgX1, double orgX2, double orgX3, 
                         int inValue, int outValue, bool bothSides,
                         bool overwriteIn, bool overwriteOut,
                         int noInValue, int noOutValue)
{
   int nx1 = (int)geoMatrix.getNX1();
   int nx2 = (int)geoMatrix.getNX2();
   int nx3 = (int)geoMatrix.getNX3();

   for(int k=0 ; k<nx3 ; k++)
   {
      for(int j=0 ; j<nx2 ; j++)
      {
         for(int i=0 ; i<nx1 ; i++)
         {
            double x = orgX1 + i*delta;
            double y = orgX2 + j*delta;
            double z = orgX3 + k*delta;

            int temp = 0;
            int gm = geoMatrix(i,j,k);

            if(geoObject->isPointInGbObject3D(x, y, z))  
            {
               setGeoType(temp, inValue);
               if (overwriteIn)
               {
                  geoMatrix(i,j,k) = temp;
               }
               else
               {
                  if(gm != noInValue) 
                  {
                     geoMatrix(i,j,k) = temp;
                  }
               }
            }
            else if(bothSides)
            {
               setGeoType(temp, outValue);
               if (overwriteOut)
               {
                  geoMatrix(i,j,k) = temp;
               }
               else
               {
                  if(gm != noOutValue) 
                  {
                     geoMatrix(i,j,k) = temp;
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
bool hasNeighbor(VoxelMatrix& geoMatrix, int x1, int x2, int x3)
{
   bool result = false;
   for( int dir = 0; dir < 18; dir++)
   {
      int temp = geoMatrix(x1+EX1[dir], x2+EX2[dir], x3+EX3[dir]);
      if(temp == GEO_BANANAS || temp == GEO_FLUID_IN_HULL)
      {
         result = true;
         break;
      }
   }
   return result;
}
//////////////////////////////////////////////////////////////////////////
void createHull(VoxelMatrix& geoMatrix)
{
   int nx1 = (int)geoMatrix.getNX1();
   int nx2 = (int)geoMatrix.getNX2();
   int nx3 = (int)geoMatrix.getNX3();

   for(int k=1 ; k<nx3-1 ; k++)
   {
      for(int j=1 ; j<nx2-1 ; j++)
      {
         for(int i=1 ; i<nx1-1 ; i++)
         {
            int val = geoMatrix(i,j,k);
            if(val == GEO_FLUID)
            {
               if(hasNeighbor(geoMatrix, i, j, k))
               {
                  int temp = 0;
                  setGeoType(temp, GEO_HULL);
                  geoMatrix(i,j,k) = temp;
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void writeGbVoxelMatrix3DtoVtuXmlASCII(const std::string& fileName, GbVoxelMatrix3DPtr voxelMatrix, 
                                       double worldDeltaX1, double worldDeltaX2, double worldDeltaX3,
                                       int nx1, int nx2, int nx3)
{
   std::vector< UbTupleFloat3 > nodes;
   std::vector<std::string > datanames;
   std::vector<std::vector<double > > nodedata;
   
   datanames.resize(0);
   datanames.push_back("Solid");
   nodes.resize(0);
   nodedata.resize(datanames.size());

   double orgX1 = voxelMatrix->getX1Minimum();
   double orgX2 = voxelMatrix->getX2Minimum();
   double orgX3 = voxelMatrix->getX3Minimum();

   int index = 0;
   double x1KO,x2KO,x3KO;
   
      for (int x3=0; x3<nx3;x3++){
         for (int x2=0; x2<nx2;x2++){
            for(int x1=0; x1<nx1;x1++)
            {
                  x1KO = orgX1 + worldDeltaX1*(double)x1;
                  x2KO = orgX2 + worldDeltaX2*(double)x2;
                  x3KO = orgX3 + worldDeltaX3*(double)x3;
                  nodes.push_back( makeUbTuple(float(x1KO), float(x2KO), float(x3KO)) );
                  nodedata[0].push_back((*voxelMatrix)(x1,x2,x3));
            }
         }
      }
   WbWriterVtkXmlASCII::getInstance()->writeNodesWithNodeData(fileName, nodes,  datanames, nodedata); 
}
//////////////////////////////////////////////////////////////////////////
void writeGbVoxelMatrix3DtoLegacyVTK(const std::string& fileName, GbVoxelMatrix3DPtr voxelMatrix,
                                       double worldDeltaX1, double worldDeltaX2, double worldDeltaX3,
                                       int nx1, int nx2, int nx3)
{
   UbFileOutputASCII out(fileName);

   int nn = nx1*nx2*nx3;
   out.writeLine("# vtk DataFile Version 3.0");
   out.writeLine(fileName);
   out.writeLine("ASCII");
   out.writeLine("DATASET STRUCTURED_POINTS");
   out.writeString("DIMENSIONS");
   out.writeInteger(nx1);
   out.writeInteger(nx2);
   out.writeInteger(nx3);
   out.writeLine();
   out.writeString("ORIGIN");
   out.writeDouble(voxelMatrix->getX1Minimum());
   out.writeDouble(voxelMatrix->getX2Minimum());
   out.writeDouble(voxelMatrix->getX3Minimum());
   out.writeLine();
   out.writeString("SPACING");
   out.writeDouble(worldDeltaX1);
   out.writeDouble(worldDeltaX2);
   out.writeDouble(worldDeltaX3);
   out.writeLine();
   out.writeString("POINT_DATA");
   out.writeInteger(nn);
   out.writeLine();
   out.writeLine("SCALARS Geo integer");
   out.writeLine("LOOKUP_TABLE default");

   for(int k=0 ; k<nx3 ; k++){
      for(int j=0 ; j<nx2 ; j++){
         for(int i=0 ; i<nx1 ; i++){
            out.writeInteger( (int)(*voxelMatrix)(i,j,k) );
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void setNormalsOnBoundary(int minX1, int minX2, int minX3, int maxX1, int maxX2, int maxX3, VoxelMatrix& matrix, int dir)
{
   for (int ix3 = minX3; ix3 <= maxX3; ix3++)
      for (int ix2 = minX2; ix2 <= maxX2; ix2++)
         for (int ix1 = minX1; ix1 <= maxX1; ix1++)
         {

            int temp = 0;
            temp = getGeoType(matrix(ix1, ix2, ix3));
            setGeoNormal(temp, dir);
            matrix(ix1, ix2, ix3) = temp;
         }
}
//////////////////////////////////////////////////////////////////////////
void run(const char *cstr)
{
   try
   {
      std::string opt;
      if(cstr!= NULL)
         opt = std::string(cstr);
      else
      {
         UBLOG(logINFO,"no option: x, y or z");
         return;
      }

      string pathnameGeo = "/home/koskuche/data/bananas";
      string pathname;

      if(opt == "z") pathname = "/work/koskuche/scratch/bananas2/setupZ";

      if(opt == "x") pathname = "/work/koskuche/scratch/bananas2/setupX";

      if(opt == "y") pathname = "/work/koskuche/scratch/bananas2/setupY";

      CommunicatorPtr comm(new MPICommunicator());
     
      //////////////////////////////////////////////////////////////////////////
      // Geometries
      //////////////////////////////////////////////////////////////////////////
      //bananas box geometry
      UBLOG(logINFO,"Start read bananas box geometry");
      GbTriFaceMesh3DPtr bananaBox (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathnameGeo+"/Banana_boxD.stl","banana_box"));
      UBLOG(logINFO,"Stop read bananas box geometry");
      bananaBox->rotate(90.0, 0.0, 0.0); //around Z

      double b_minX1 = bananaBox->getX1Minimum();
      double b_minX2 = bananaBox->getX2Minimum();
      double b_minX3 = bananaBox->getX3Minimum();

      double b_maxX1 = bananaBox->getX1Maximum();
      double b_maxX2 = bananaBox->getX2Maximum();
      double b_maxX3 = bananaBox->getX3Maximum();

      if(opt == "x") bananaBox->rotate(0.0, 0.0, -90.0); //around X

      if(opt == "y") bananaBox->rotate(0.0, -90.0, 0.0); //around Y

      //after rotation for setup 2-3
      double bb_minX1 = bananaBox->getX1Minimum();
      double bb_minX2 = bananaBox->getX2Minimum();
      double bb_minX3 = bananaBox->getX3Minimum();

      double bb_maxX1 = bananaBox->getX1Maximum();
      double bb_maxX2 = bananaBox->getX2Maximum();
      double bb_maxX3 = bananaBox->getX3Maximum();

      UBLOG(logINFO,"Start write bananas box geometry");
      GbSystem3D::writeGeoObject(bananaBox.get(), pathname+"/banana_box", WbWriterVtkXmlASCII::getInstance());
      UBLOG(logINFO,"Stop write bananas box geometry");

      //distances for bounding box
      double dist_z = 0.022;
      double site   = 0.011;

      //bounding box of simulation
      double g_minX1 = bananaBox->getX1Minimum()-site;
      double g_minX2 = bananaBox->getX2Minimum()-site;
      double g_minX3 = bananaBox->getX3Minimum()-dist_z*2.0;

      double g_maxX1 = bananaBox->getX1Maximum()+site;
      double g_maxX2 = bananaBox->getX2Maximum()+site;
      double g_maxX3 = bananaBox->getX3Maximum()+dist_z*2.0;

      const double gridOriginX1 = g_minX1;
      const double gridOriginX2 = g_minX2;
      const double gridOriginX3 = g_minX3;

      const double dx = 2.20183486239e-3;
      UBLOG(logINFO,"DeltaX = " << dx);

      GbCuboid3DPtr addWall1 (new GbCuboid3D(g_minX1, g_minX2, bb_minX3, bb_minX1, g_maxX2, bb_minX3+2*dx));
      GbSystem3D::writeGeoObject(addWall1.get(), pathname+"/addWall1", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWall2 (new GbCuboid3D(bb_maxX1, g_minX2, bb_minX3, g_maxX1, g_maxX2, bb_minX3+2*dx));
      GbSystem3D::writeGeoObject(addWall2.get(), pathname+"/addWall2", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWall3 (new GbCuboid3D(g_minX1, g_minX2, bb_minX3, g_maxX1, bb_minX2, bb_minX3+2*dx));
      GbSystem3D::writeGeoObject(addWall3.get(), pathname+"/addWall3", WbWriterVtkXmlASCII::getInstance());

      GbCuboid3DPtr addWall4 (new GbCuboid3D(g_minX1, bb_maxX2, bb_minX3, g_maxX1, g_maxX2, bb_minX3+2*dx));
      GbSystem3D::writeGeoObject(addWall4.get(), pathname+"/addWall4", WbWriterVtkXmlASCII::getInstance());

      VoxelMatrix grid(int((g_maxX1-g_minX1)/dx)+1, int((g_maxX2-g_minX2)/dx)+1, int((g_maxX3-g_minX3)/dx)+1, GEO_FLUID);

      UBLOG(logINFO,"Start write geo matrix empty");
      writeMatrixToVtkImageFile(pathname + "/geo_matrix_empty.vtk", grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop write geo matrix empty");

#ifdef CONVEXHULL
      UBLOG(logINFO,"Start read bananas box geometry");
      GbTriFaceMesh3DPtr bananasHull (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathnameGeo+"/convexhullASCII.stl","banana_hull"));
      UBLOG(logINFO,"Stop read bananas box geometry");
      double tr1 = bananasHull->getX3Minimum() - bananaBox->getX3Minimum();
      bananasHull->translate(0.0, 0.0, 5.0*dx);
#endif

      //reed bananas
      UBLOG(logINFO,"Start read bananas geometry");
      int d1, d2, d3;
      readDimensionsFromFldFile(pathnameGeo + "/BANANA_8binn_Binear_A.fld", d1, d2, d3);
      UBLOG(logINFO,"Dimensions of bananas geometry: " << d1 << ", " << d2 << ", " << d3);
      GbVoxelMatrix3DPtr bananas(new GbVoxelMatrix3D(d1, d2, d3, float(GbVoxelMatrix3D::FLUID))); 
      readGeoMatrixFromFldFile(pathnameGeo + "/BANANA_8binn_Binear_A.fld", bananas);
      UBLOG(logINFO,"Stop read bananas geometry");
      double bananasDx1 = (b_maxX1 - b_minX1) / float(d1);
      double bananasDx2 = (b_maxX2 - b_minX2) / float(d2);
      double bananasDx3 = (b_maxX3 - b_minX3) / float(d3);
      bananas->setVoxelMatrixDelta(float(bananasDx1), float(bananasDx2), float(bananasDx3));
      bananas->setCenterCoordinates(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
      bananas->setVoxelMatrixMininum(float(b_minX1), float(b_minX2), float(b_minX3));
      
      bananas->rotate90aroundY();
      bananas->rotate90aroundY();
      bananas->rotate90aroundZ();
      bananas->rotate90aroundZ();
      
#ifdef CONVEXHULL
      std::cout << "translate bananas: " <<bananasHull->getX3Minimum() - bananas->getX3Minimum() - tr1<<"\n";
      bananas->translate(0.0, 0.0, bananasHull->getX3Minimum() - bananas->getX3Minimum() - tr1);
      if(opt == "x") bananasHull->rotateAroundPoint(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid(), 0.0, 0.0, -90.0); //around X
      if(opt == "y") bananasHull->rotateAroundPoint(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid(), 0.0, -90.0, 0.0); //around Y
      UBLOG(logINFO,"Start write banana hull geometry");
      GbSystem3D::writeGeoObject(bananasHull.get(), pathname+"/banana_hull", WbWriterVtkXmlASCII::getInstance());
      UBLOG(logINFO,"Stop write banana hull geometry");
#endif
      
      if(opt == "x")
      {
         bananas->rotate90aroundX(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
         bananas->rotate90aroundX(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
         bananas->rotate90aroundX(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
      }
      else if(opt == "y")
      {
         bananas->rotate90aroundY(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
         bananas->rotate90aroundY(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
         bananas->rotate90aroundY(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid());
      }

      UBLOG(logINFO,"Start write bananas geometry");
      bananas->writeToLegacyVTK(pathname + "/bananas.vtk");
      UBLOG(logINFO,"Stop write bananas geometry");

#ifdef BANANAS
      UBLOG(logINFO,"Start discretization of banana box");
      discretizeGeoObject(bananaBox, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX);
      UBLOG(logINFO,"Stop discretization of banana box");
      UBLOG(logINFO,"Start discretization of bananas");
      discretizeGeoObject(bananas, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BANANAS);
      UBLOG(logINFO,"Stop discretization of bananas");
#endif

#ifdef CONVEXHULL
      UBLOG(logINFO,"Start discretization of bananas");
      discretizeGeoObject(bananas, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BANANAS, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      UBLOG(logINFO,"Stop discretization of bananas");
      UBLOG(logINFO,"Start discretization of hull");
      discretizeGeoObject(bananasHull, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_FLUID_IN_HULL, GEO_FLUID, true, false, true, GEO_BANANAS, GEO_INVALID);
      UBLOG(logINFO,"Stop discretization of hull");
      UBLOG(logINFO,"Start creation of hull film");
      createHull(grid);
      UBLOG(logINFO,"Stop creation of hull film");
      UBLOG(logINFO,"Start discretization of banana box");
      discretizeGeoObject(bananaBox, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      UBLOG(logINFO,"Stop discretization of banana box");
#endif

      UBLOG(logINFO,"Start discretization of add walls");
      discretizeGeoObject(addWall1, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      discretizeGeoObject(addWall2, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      discretizeGeoObject(addWall3, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      discretizeGeoObject(addWall4, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3, GEO_BOX, GEO_FLUID, false, true, true, GEO_INVALID, GEO_INVALID);
      UBLOG(logINFO,"Stop discretization of add walls");

      UBLOG(logINFO,"Start set normals");
      int boxNX1 = int(bananaBox->getLengthX1() / dx);
      int boxNX2 = int(bananaBox->getLengthX2() / dx);
      int boxNX3 = int(bananaBox->getLengthX3() / dx);

      int minX1 = int((bb_minX1 - gridOriginX1) / dx)+1;
      int minX2 = int((bb_minX2 - gridOriginX2) / dx)+1;
      int minX3 = int((bb_minX3 - gridOriginX3) / dx)+1;

      int maxX1 = minX1 + boxNX1;
      int maxX2 = minX2 + boxNX2;
      int maxX3 = minX3 + boxNX3;

      UBLOG(logINFO,"minX1="<<minX1<<",minX2= "<<minX2<<",minX3="<<minX3);
      UBLOG(logINFO,"maxX1="<<maxX1<<",maxX2= "<<maxX2<<",maxX3="<<maxX3);


      for (int ix3 = 0; ix3 < grid.getNX3(); ix3++)
         for (int ix2 = 0; ix2 < grid.getNX2(); ix2++)
            for (int ix1 = 0; ix1 < grid.getNX1(); ix1++)
            {
               int temp = grid(ix1, ix2, ix3);
               setGeoNormal(temp, 0);
               grid(ix1, ix2, ix3) = temp;
            }


      setNormalsOnBoundary(minX1, minX2, minX3, minX1, maxX2, maxX3, grid, NORMAL_NEG_X1);
      setNormalsOnBoundary(maxX1, minX2, minX3, maxX1, maxX2, maxX3, grid, NORMAL_POS_X1);
      setNormalsOnBoundary(minX1, minX2, minX3, maxX1, minX2, maxX3, grid, NORMAL_NEG_X2);
      setNormalsOnBoundary(minX1, maxX2, minX3, maxX1, maxX2, maxX3, grid, NORMAL_POS_X2);
      setNormalsOnBoundary(minX1, minX2, minX3, maxX1, maxX2, minX3, grid, NORMAL_NEG_X3);
      setNormalsOnBoundary(minX1, minX2, maxX3, maxX1, maxX2, maxX3, grid, NORMAL_POS_X3);
      UBLOG(logINFO,"Stop set normals");


      UBLOG(logINFO,"Start write geo matrix");
      writeMatrixToVtkImageFile(pathname + "/geo_matrix.vtk", grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop write geo matrix");
   }
   catch(std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch(std::string& s)
   {
      cerr << s << endl;
   }
   catch(...)
   {
      cerr << "unknown exception" << endl;
   }

}
int main(int argc, char* argv[])
{

   run(argv[1]);

   return 0;
}

