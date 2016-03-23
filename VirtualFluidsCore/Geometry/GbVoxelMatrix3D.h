//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBVOXELMATRIX3D_H
#define GBVOXELMATRIX3D_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <vector>
#include <cmath>

#include <numerics/geometry3d/GbObject3D.h>
#include <basics/utilities/UbObserver.h>
#include <basics/container/CbArray3D.h>

#include <basics/memory/MbSharedPointerDefines.h>
class GbVoxelMatrix3D;
typedef VFSharedPtr<GbVoxelMatrix3D> GbVoxelMatrix3DPtr;


class GbLine3D;                    
class GbTriangle3D;                    
class GbObject3DCreator;

class GbVoxelMatrix3D : public GbObject3D, public UbObserver
{
public:              
   typedef CbArray3D<float> Matrix3D;
   static const float SOLID; 
   static const float FLUID; 
   enum  Endian {BigEndian, LittleEndian};

   GbVoxelMatrix3D();
   GbVoxelMatrix3D(int nx1, int nx2, int nx3, float initVal, double lowerThreshold=0, double upperThreshold=0);
   GbVoxelMatrix3D(const Matrix3D& voxelMatrix, double lowerThreshold=0, double upperThreshold=0);
   ~GbVoxelMatrix3D() {}   

   void finalize() {};
   GbVoxelMatrix3D* clone(); 

   /*=======================================================================*/
   Matrix3D::reference operator() (const Matrix3D::size_type& x1, const Matrix3D::size_type& x2, const Matrix3D::size_type& x3)
   {
      return voxelMatrix(x1,x2,x3);
   }
   /*=======================================================================*/
   Matrix3D::const_reference operator() (const Matrix3D::size_type& x1, const Matrix3D::size_type& x2, const Matrix3D::size_type& x3)	const
   {
      return voxelMatrix(x1,x2,x3);
   }
   /*=======================================================================*/
   void setTransferViaFilename(bool transferViaFilename, std::string filename)
   {
      this->filename = filename;
      this->transferViaFilename = transferViaFilename;
   }
   void setThreshold(double lowerThreshold, double upperThreshold) { this->lowerThreshold = lowerThreshold; this->upperThreshold = upperThreshold; }
   void setAddSurfaceTriangleSetFlag(bool flag) { this->addSurfaceTriangleSetFlag = flag; }

   /*=======================================================================*/
   void setVoxelMatrixMininum(double minX1, double minX2, double minX3) { this->minX1 = minX1; this->minX2 = minX2; this->minX3 = minX3; }
   void setVoxelMatrixMinX1(double minX1) { this->minX1 = minX1; }
   void setVoxelMatrixMinX2(double minX2) { this->minX2 = minX2; }
   void setVoxelMatrixMinX3(double minX3) { this->minX3 = minX3; }
   
   /*=======================================================================*/
   void setVoxelMatrixDelta(double deltaX1, double deltaX2, double deltaX3) { this->deltaX1 = deltaX1; this->deltaX2 = deltaX2; this->deltaX3 = deltaX3; }
   void setVoxelMatrixDeltaX1(double deltaX1) { this->deltaX1 = deltaX1; }
   void setVoxelMatrixDeltaX2(double deltaX2) { this->deltaX2 = deltaX2; }
   void setVoxelMatrixDeltaX3(double deltaX3) { this->deltaX3 = deltaX3; }

   /*=======================================================================*/
   double getX1Centroid() { return 0.5 * ( minX1 + this->getX1Maximum() ); } 
   double getX1Minimum()  { return minX1; }
   double getX1Maximum()  { return minX1 + deltaX1*voxelMatrix.getNX1(); }  

   double getX2Centroid() { return 0.5 * ( minX2 + this->getX2Maximum() ); } 
   double getX2Minimum()  { return minX2; }
   double getX2Maximum()  { return minX2 + deltaX2*voxelMatrix.getNX2(); }  

   double getX3Centroid() { return 0.5 * ( this->getX3Minimum() + this->getX3Maximum() ); } 
   double getX3Minimum()  { return minX3; }
   double getX3Maximum()  { return minX3 + deltaX3*voxelMatrix.getNX3(); }  

   double getLengthX1() { return this->getX1Maximum() - minX1; }
   double getLengthX2() { return this->getX2Maximum() - minX2; }
   double getLengthX3() { return this->getX3Maximum() - minX3; }

   void setCenterCoordinates(const double& x1, const double& x2, const double& x3);
   void translate(const double& tx1, const double& tx2, const double& tx3);

   bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary);
   bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p); 
   bool isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   //double getCellVolumeInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);

   GbPoint3D*  calculateInterSectionPoint3D(GbPoint3D& point1, GbPoint3D &point2) { throw UbException(__FILE__,__LINE__, UB_FUNCTION,"not implemented");}
   GbLine3D*   createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2){ throw UbException(__FILE__,__LINE__, UB_FUNCTION,"not implemented");}

   std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   void addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles);

   bool hasRaytracing() { return true; }
   /*|r| must be 1! einheitsvector!!*/
   double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3);

   std::string toString();

   ObObjectCreator* getCreator();
   void write(UbFileOutput* out);
   void read(UbFileInput* in);

   //virtuelle Methoden von UbObserver
   void objectChanged(UbObservable* changedObject) {}
   void objectWillBeDeleted(UbObservable* objectForDeletion) {}

   template <class T>
   void readMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian);
   template <class T>
   void readBufferedMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian);
   void readMatrixFromVtiASCIIFile(std::string filename);

   void rotate90aroundX();
   void rotate90aroundY();
   void rotate90aroundZ();
   void rotate90aroundX(double cX1, double cX2, double cX3);
   void rotate90aroundY(double cX1, double cX2, double cX3);
   void rotate90aroundZ(double cX1, double cX2, double cX3);
   void mirrorX();
   void mirrorY();
   void mirrorZ();

   void writeToLegacyVTKASCII(const std::string& fileName);
   void writeToLegacyVTKBinary( const std::string& fileName );
   void writeToVTKImageDataASCII(const std::string& fileName);
   void writeToVTKImageDataAppended(const std::string& fileName);

   void setClosedVoidSpaceToSolid();

   void calculateNumberOfSolidAndFluid();
   long getNumberOfSolid();
   long getNumberOfFluid();

protected:
   void findFluidNeighbor(int cx1, int cx2, int cx3);
   Matrix3D voxelMatrixTemp;
   CbArray3D<char> flagMatrix;

   std::vector<int> x1Nbr;
   std::vector<int> x2Nbr;
   std::vector<int> x3Nbr;

   std::vector<int> x1NbrTemp;
   std::vector<int> x2NbrTemp;
   std::vector<int> x3NbrTemp;

   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & minX1;
      ar & minX2;
      ar & minX3;
      ar & deltaX1;
      ar & deltaX2;
      ar & deltaX3;
      ar & nodesX1;
      ar & nodesX2;
      ar & nodesX3;
      ar & threshold;
      ar & transferViaFilename;
      ar & addSurfaceTriangleSetFlag;
      if(!transferViaFilename)
      {
         ar & voxelMatrix;
      }
      else
      {
         ar & filename;
         if(ArchiveTools::isReading(ar) ) 
         {
            this->readMatrixFromVtiASCIIFile(filename);
         }
      }

   }
#endif //CAB_RCF

protected:
   //for transfer
   std::string filename;
   bool transferViaFilename;

   bool addSurfaceTriangleSetFlag;

   int nodesX1;
   int nodesX2;
   int nodesX3;
   double lowerThreshold, upperThreshold;

   double minX1;
   double minX2;
   double minX3;
   double deltaX1;
   double deltaX2;
   double deltaX3;

   Matrix3D voxelMatrix;

   long numberOfSolid;
   long numberOfFluid;
};

//////////////////////////////////////////////////////////////////////////
template <class T>
void GbVoxelMatrix3D::readMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian)
{
   using namespace std;
   UBLOG(logINFO,"GbVoxelMatrix3D::readMatrixFromFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");
   ifstream in(filename.c_str(), ios::binary);
   if(!in) throw UbException(UB_EXARGS,"could not open file "+filename);

   in.seekg( 0, ios::end );     //Ende springen
   fstream::off_type length = in.tellg(); //Position abfragen
   in.seekg( 0, ios::beg );    //An den Anfang springen 

   //UBLOG(logINFO,"number of nodes = "<<nodesX1*nodesX2*nodesX3*sizeof(T)<<" file size = "<<(long)length);
   //if( (nodesX1*nodesX2*nodesX3)*sizeof(float) != (long)length )
   unsigned long long nofn = nodesX1*nodesX2*nodesX3*sizeof(T);
   if (nofn != (unsigned long long)length)
   {
      throw UbException(UB_EXARGS,"number of nodes("+UbSystem::toString(nofn)+") doesn't match file size("+UbSystem::toString((long)length)+")");
   }

   UBLOG(logINFO,"  - create GbVoxelMatrix3D");
   //GbVoxelMatrix3D* voxelGeo = new GbVoxelMatrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID);
   voxelMatrix = Matrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID);

   UBLOG(logINFO,"  - init values");
   //float val;
   T val;
   for(int x3=0; x3<nodesX3; x3++)
      for(int x2=0; x2<nodesX2; x2++)
         for(int x1=0; x1<nodesX1; x1++)
         {
            //in.read((char*)&val,sizeof(float));
            in.read((char*)&val,sizeof(T));
            if (endian == BigEndian)
               UbSystem::swapByteOrder((unsigned char*)(&(val)), sizeof(T));
            //if( UbMath::equal((double)val, threshold) ) 
            //if( UbMath::greater((double)val, threshold) )
            if( (double)val >= lowerThreshold && (double)val <= upperThreshold ) 
            {
               (voxelMatrix)(x1,x2,x3) = GbVoxelMatrix3D::SOLID;
            }
            //(voxelMatrix)(x1, x2, x3) = (float)val;
         }

         UBLOG(logINFO,"GbVoxelMatrix3D::readMatrixFromFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");
}

//////////////////////////////////////////////////////////////////////////
template <class T>
void GbVoxelMatrix3D::readBufferedMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian)
{
   using namespace std;
   UBLOG(logINFO, "GbVoxelMatrix3D::readMatrixFromRawFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");

   FILE *file;
   file = fopen(filename.c_str(), "rb");
   if (file==NULL)
   {
      throw UbException(UB_EXARGS, "Could not open file "+filename);
   }

   // obtain file size:
   fseek(file, 0, SEEK_END);
   unsigned long int length = ftell(file);
   rewind(file);

   UBLOG(logINFO, "number of nodes = "<<(long)nodesX1*(long)nodesX2*(long)nodesX3<<" file size = "<<length);
   
   unsigned long int nofn = (long)nodesX1*(long)nodesX2*(long)nodesX3*(long)sizeof(T);
   if (nofn != length)
   {
      //throw UbException(UB_EXARGS, "number of nodes("+UbSystem::toString(nofn)+") doesn't match file size("+UbSystem::toString(length)+")");
   }

   UBLOG(logINFO, "  - create GbVoxelMatrix3D");
   voxelMatrix = Matrix3D(nodesX1, nodesX2, nodesX3, GbVoxelMatrix3D::FLUID);

   CbArray3D<T> readMatrix(nodesX1, nodesX2, nodesX3);

   UBLOG(logINFO, "  - read file to matrix");
   fread(readMatrix.getStartAdressOfSortedArray(0, 0, 0), sizeof(T), readMatrix.getDataVector().size(), file);
   fclose(file);

   UBLOG(logINFO, "  - init values");

   numberOfSolid = 0;
   T val;
   for (int x3 = 0; x3<nodesX3; x3++)
      for (int x2 = 0; x2<nodesX2; x2++)
         for (int x1 = 0; x1<nodesX1; x1++)
         {
            val = readMatrix(x1, x2, x3);

            if (endian == BigEndian)
            {
               UbSystem::swapByteOrder((unsigned char*)(&(val)), sizeof(T));
            }

            if ((double)val >= lowerThreshold && (double)val <= upperThreshold)
            {
               voxelMatrix(x1, x2, x3) = GbVoxelMatrix3D::SOLID;
            }
         }

   UBLOG(logINFO, "GbVoxelMatrix3D::readMatrixFromRawFile \""<<filename<<"\" nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");
}


#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   UB_AUTO_RUN_NAMED(   SF::registerType<GbVoxelMatrix3D>("GbVoxelMatrix3D")        , SF_GbVoxelMatrix3D     );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbVoxelMatrix3D >()), SF_GbVoxelMatrix3D_BD1 );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< UbObserver, GbVoxelMatrix3D>() ), SF_GbVoxelMatrix3D_BD2 );
#endif //RCF_USE_SF_SERIALIZATION


#endif   
