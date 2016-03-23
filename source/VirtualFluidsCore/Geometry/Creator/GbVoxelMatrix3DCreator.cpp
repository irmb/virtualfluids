#include <numerics/geometry3d/creator/GbVoxelMatrix3DCreator.h>
#include <numerics/geometry3d/GbVoxelMatrix3D.h>
#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbLogger.h>

using namespace std;

/***************************************************************************/
GbVoxelMatrix3D* GbVoxelMatrix3DCreator::createFromRawFloatFile(  string filename, int nodesX1, int nodesX2, int nodesX3, float threshold)
{
   UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromRawFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");
   ifstream in(filename.c_str(), ios::binary);
   if(!in) throw UbException(UB_EXARGS,"could not open file "+filename);
   
   in.seekg( 0, ios::end );     //Ende springen
   fstream::off_type length = in.tellg(); //Position abfragen
   in.seekg( 0, ios::beg );    //An den Anfang springen 
   if( (nodesX1*nodesX2*nodesX3)*sizeof(float) != (long)length )
   {
      throw UbException(UB_EXARGS,"number of nodes doesn't match filesize");
   }

   UBLOG(logINFO,"  - create GbVoxelMatrix3D");
   GbVoxelMatrix3D* voxelGeo = new GbVoxelMatrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID, threshold);
   
   UBLOG(logINFO,"  - init values");
   float val;
   for(int x3=0; x3<nodesX3; x3++)
      for(int x2=0; x2<nodesX2; x2++)
         for(int x1=0; x1<nodesX1; x1++)
         {
            in.read((char*)&val,sizeof(float));
            //if( !UbMath::equal(val, 0.0f) ) 
            if( UbMath::greater(val, threshold) ) 
            {
               (*voxelGeo)(x1,x2,x3) = GbVoxelMatrix3D::SOLID;
            }
         }
   
   UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromRawFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");

   return voxelGeo;
}
/***************************************************************************/
GbVoxelMatrix3D* GbVoxelMatrix3DCreator::createFromVtiASCIIFloatFile(  string filename, int nodesX1, int nodesX2, int nodesX3, float threshold)
{
   UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromVtiASCIIFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");
   UbFileInputASCII in(filename);
   //ifstream in(filename.c_str(), ios::binary);
   if(!in) throw UbException(UB_EXARGS,"could not open file "+filename);
   in.readLine();
   in.readLine();
   in.readLine();
   in.readLine();
   in.readLine();
   //in.readLine(); !!!manchmal hat das vti file noch die xml version dabei ...

   UBLOG(logINFO,"  - create GbVoxelMatrix3D");
   GbVoxelMatrix3D* voxelGeo = new GbVoxelMatrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID, threshold);

   UBLOG(logINFO,"  - init values");
   int val;
   int u=0;
   for(int x3=0; x3<nodesX3; x3++)
      for(int x2=0; x2<nodesX2; x2++)
         for(int x1=0; x1<nodesX1; x1++)
         {
            val = in.readInteger();
            
            //u++; if(u>125000) UBLOG(logINFO,"val:"<<u<<" "<<val);

            //if( !UbMath::equal(val, 0.0f) ) 
            if( UbMath::greater(val, threshold) ) 
            {
               (*voxelGeo)(x1,x2,x3) = GbVoxelMatrix3D::SOLID;
            }
         }

         UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromVtiASCIIFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");

         return voxelGeo;
}

