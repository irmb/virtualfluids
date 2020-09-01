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

#define CONVEXHULL

using namespace std;

const int FLUID = 1;
const int SOLID = 15;

//////////////////////////////////////////////////////////////////////////
void writeMatrixToVtkImageFile(const std::string& fileName, const CbArray3D <int>& geoMatrix,
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
void discretizeGeoObject(GbObject3DPtr geoObject, CbArray3D<int>& geoMatrix, double delta, double orgX1, double orgX2, double orgX3)
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
            if(geoObject->isPointInGbObject3D(x, y, z)) geoMatrix(i,j,k) = SOLID;
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
void run(const char *cstr)
{
   try
   {
      //string pathname = "c:/temp/bananas/out";
      //string pathnameGeo = "c:/temp/bananas/geo";

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

      if(opt == "z") pathname = "/work/koskuche/scratch/bananas/setupZ/out";

      if(opt == "x") pathname = "/work/koskuche/scratch/bananas/setupX/out";

      if(opt == "y") pathname = "/work/koskuche/scratch/bananas/setupY/out";

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

      UBLOG(logINFO,"Start write bananas box geometry");
      GbSystem3D::writeGeoObject(bananaBox.get(), pathname+"/banana_box", WbWriterVtkXmlASCII::getInstance());
      UBLOG(logINFO,"Stop write bananas box geometry");

      //distances for bounding box
      double dist_z = 0.12;
      double dist_x = 0.26;
      double dist_y = 0.195;

      double g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3;

      //bounding box of simulation
      //setup1 - z
      if(opt == "z")
      {
         g_minX1 = bananaBox->getX1Minimum();
         g_minX2 = bananaBox->getX2Minimum();
         g_minX3 = bananaBox->getX3Minimum()-dist_z;

         g_maxX1 = bananaBox->getX1Maximum();
         g_maxX2 = bananaBox->getX2Maximum();
         g_maxX3 = bananaBox->getX3Maximum()+dist_z*2.0;
      }

      //setup2 - x
      if(opt == "x")
      {
         g_minX1 = bananaBox->getX1Minimum();
         g_minX2 = bananaBox->getX2Minimum();
         g_minX3 = bananaBox->getX3Minimum()-dist_x;

         g_maxX1 = bananaBox->getX1Maximum();
         g_maxX2 = bananaBox->getX2Maximum();
         g_maxX3 = bananaBox->getX3Maximum()+dist_x*2.0;
      }

      //setup3 - y
      if(opt == "y")
      {
         g_minX1 = bananaBox->getX1Minimum();
         g_minX2 = bananaBox->getX2Minimum();
         g_minX3 = bananaBox->getX3Minimum()-dist_y;

         g_maxX1 = bananaBox->getX1Maximum();
         g_maxX2 = bananaBox->getX2Maximum();
         g_maxX3 = bananaBox->getX3Maximum()+dist_y*2.0;
      }

      const double gridOriginX1 = g_minX1;
      const double gridOriginX2 = g_minX2;
      const double gridOriginX3 = g_minX3;

      //int gridNX1 = 170;
      //int gridNX2 = 226;
      //int gridNX3 = 104;

      const double dx = 2.20183486239e-3; //blockLentghX1/static_cast<double>(blocknx1);

      UBLOG(logINFO,"DeltaX = " << dx);

      CbArray3D<int> grid(int((g_maxX1-g_minX1)/dx)+1, int((g_maxX2-g_minX2)/dx)+1, int((g_maxX3-g_minX3)/dx)+1, FLUID);

      UBLOG(logINFO,"Start write geo matrix empty");
      writeMatrixToVtkImageFile(pathname + "/geo_matrix_empty.vtk", grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop write geo matrix empty");

#ifdef BANANAS
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

      bananas->rotate90aroundX();
      bananas->rotate90aroundY();
      //bananas->rotate90aroundX();

      UBLOG(logINFO,"Start write bananas geometry");
      bananas->writeToLegacyVTK(pathname + "/bananas.vtk");
      UBLOG(logINFO,"Stop write bananas geometry");
#endif

#ifdef CONVEXHULL
      UBLOG(logINFO,"Start read bananas box geometry");
      GbTriFaceMesh3DPtr bananaHull (GbTriFaceMesh3DCreator::getInstance()->readMeshFromSTLFile(pathnameGeo+"/convexhullASCII.stl","banana_hull"));
      UBLOG(logINFO,"Stop read bananas box geometry");
      bananaHull->translate(0.0, 0.0, 5.0*dx);
      if(opt == "x") bananaHull->rotateAroundPoint(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid(), 0.0, 0.0, -90.0); //around X
      if(opt == "y") bananaHull->rotateAroundPoint(bananaBox->getX1Centroid(), bananaBox->getX2Centroid(), bananaBox->getX3Centroid(), 0.0, -90.0, 0.0); //around Y
      UBLOG(logINFO,"Start write banana hull geometry");
      GbSystem3D::writeGeoObject(bananaHull.get(), pathname+"/banana_hull", WbWriterVtkXmlASCII::getInstance());
      UBLOG(logINFO,"Stop write banana hull geometry");
#endif
      ////////////////////////////////////////
      //return;
      /////////////////////////////////////////

      UBLOG(logINFO,"Start discretization of banana box");
      discretizeGeoObject(bananaBox, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop discretization of banana box");

#ifdef BANANAS
      UBLOG(logINFO,"Start discretization of bananas");
      discretizeGeoObject(bananas, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop discretization of bananas");
#endif

#ifdef CONVEXHULL
      UBLOG(logINFO,"Start discretization of banana hull");
      discretizeGeoObject(bananaHull, grid, dx, gridOriginX1, gridOriginX2, gridOriginX3);
      UBLOG(logINFO,"Stop discretization of banana hull");
#endif

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

