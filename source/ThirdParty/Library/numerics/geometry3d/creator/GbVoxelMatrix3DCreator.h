#ifndef GBVOXELMATRIX3DCREATOR_H
#define GBVOXELMATRIX3DCREATOR_H

#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbVoxelMatrix3D.h>     
#include <iostream>
#include <fstream>

class GbVoxelMatrix3DCreator : public GbObject3DCreator              
{               
public:
   enum DataType {t8bit, t16bit};
public:
   static GbVoxelMatrix3DCreator* getInstance()
   {
      static GbVoxelMatrix3DCreator instance;
      return &instance;
   }

   GbVoxelMatrix3D* createGbObject3D() { return new GbVoxelMatrix3D(); }          
   GbVoxelMatrix3D* createFromRawFloatFile(  std::string filename, int nodesX1, int nodesX2, int nodesX3, float threshold=0.0);
   GbVoxelMatrix3D* createFromVtiASCIIFloatFile(  std::string filename, int nodesX1, int nodesX2, int nodesX3, float threshold=0.0);

   std::string getGbObject3DTypeID() { return "GbVoxelMatrix3D"; };
   std::string toString()            { return "GbVoxelMatrix3DCreator"; }

private:
   GbVoxelMatrix3DCreator() : GbObject3DCreator() {}

   GbVoxelMatrix3DCreator( const GbVoxelMatrix3DCreator& );                  //no copy allowed 
   const GbVoxelMatrix3DCreator& operator=( const GbVoxelMatrix3DCreator& ); //no copy allowed

public:
   template< typename T >
   GbVoxelMatrix3D* createFromRawFile(std::string filename, int nodesX1, int nodesX2, int nodesX3, float threshold)
   {
      UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromRawFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");

      std::ifstream in(filename.c_str(), std::ios::binary);
      if(!in) throw UbException(UB_EXARGS,"could not open file "+filename);

      in.seekg( 0, std::ios::end );     //Ende springen
      std::fstream::off_type length = in.tellg(); //Position abfragen
      in.seekg( 0, std::ios::beg );    //An den Anfang springen 
      long m_size = (nodesX1*nodesX2*nodesX3)*sizeof(T);
      if( m_size != (long)length )
      {
         throw UbException(UB_EXARGS,"number of nodes doesn't match filesize: " + UbSystem::toString(length));
      }

      UBLOG(logINFO,"  - create GbVoxelMatrix3D");
      GbVoxelMatrix3D* voxelGeo = new GbVoxelMatrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID, threshold);

      UBLOG(logINFO,"  - init values");
      T val;
      for(int x3=0; x3<nodesX3; x3++)
         for(int x2=0; x2<nodesX2; x2++)
            for(int x1=0; x1<nodesX1; x1++)
            {
               in.read((char*)&val,sizeof(T));
               //if( !UbMath::equal(val, 0.0f) ) 
               //if( UbMath::greater(val, (T)threshold) ) 
               if(val > (T)threshold)
               {
                  (*voxelGeo)(x1,x2,x3) = GbVoxelMatrix3D::SOLID;
               }
            }

      UBLOG(logINFO,"GbVoxelMatrix3DCreator::createFromRawFloatFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");

      return voxelGeo;
   }
};



#ifndef SWIG
UB_AUTO_RUN_NAMED( GbObject3DFactory::getInstance()->addObObjectCreator(GbVoxelMatrix3DCreator::getInstance()), CAB_GbVoxelMatrix3DCreator);
#endif

#endif  //GBVOXELMATRIX3DCREATOR_H 
