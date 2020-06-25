#include "VtkGeometryReader.h"

#include "VirtualFluidsBasics/basics/utilities/UbFileInputASCII.h"

void VtkGeometryReader::readFile(const std::string& fileName, unsigned int* geoMat)
{
   UbFileInputASCII in(fileName);

   //in.readLine();
   //in.readLine();
   //in.readLine();
   //in.readLine();
   //in.readString();
   //int itsNX1 = in.readInteger();
   //int itsNX2 = in.readInteger();
   //int itsNX3 = in.readInteger();
   //int nn = itsNX1 * itsNX2 * itsNX3;
   //in.readLine();
   //in.readString();
   //in.readInteger();
   //in.readString();
   //in.readLine();

   //for(int k=0 ; k<itsNX3 ; k++){
   //   for(int j=0 ; j<itsNX2 ; j++){
   //      for(int i=0 ; i<itsNX1 ; i++){
   //         //float x =  i;//(j&0x1)^(k&0x1)  + i*2;
   //         //float y =  j;
   //         //float z =  k;
   //         in.readFloat();
   //         in.readFloat();
   //         in.readFloat();
   //         in.readLine();
   //      }
   //   }
   //}


   //in.readString();
   //in.readInteger();
   //in.readLine();

   //in.readLine();
   //in.readLine();

   //for(int k=0 ; k<itsNX3 ; k++){
   //   for(int j=0 ; j<itsNX2 ; j++){
   //      for(int i=0 ; i<itsNX1 ; i++){
   //         int m = itsNX1*(itsNX2*k + j) + i;
   //         geoMat[m] = in.readInteger();
   //         in.readLine();
   //      }
   //   }
   //} 
   //int nn = itsNX1 * itsNX2 * itsNX3;
   in.readLine();
   in.readLine();
   in.readLine();
   in.readLine();
   in.readString();
   unsigned int itsNX1 = in.readInteger();
   unsigned int itsNX2 = in.readInteger();
   unsigned int itsNX3 = in.readInteger();
      //int nn = itsNX1 * itsNX2 * itsNX3;

   std::cout << "NX = " << itsNX1 << "\n";
   std::cout << "NY = " << itsNX2 << "\n";
   std::cout << "NZ = " << itsNX3 << "\n";

   in.readLine();
   in.readString();
   in.readDouble();
   in.readDouble();
   in.readDouble();
   //in.readInteger();
   //in.readInteger();
   //in.readInteger();
   in.readLine();
   in.readString();
   in.readDouble();
   in.readDouble();
   in.readDouble();
   in.readLine();
   in.readString();
   in.readInteger();
   in.readLine();
   in.readLine();
   in.readLine();

   for(unsigned int k=0 ; k<itsNX3 ; k++){
      for(unsigned int j=0 ; j<itsNX2 ; j++){
         for(unsigned int i=0 ; i<itsNX1 ; i++){
            //int m = nx*(ny*k + j) + i;
            //in.readInteger( itsGeoMatrix(i,j,k) );
            //in.readLine();
                     unsigned int m = itsNX1*(itsNX2*k + j) + i;
                     geoMat[m] = in.readInteger();
                     //if (i==0 && j==0 && k>50 && k<60)
                     //{
                     //   std::cout << "Read = " << geoMat[m] << "\n";                        
                     //}
         }
      }
   }

}