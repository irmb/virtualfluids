//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MATRIX4DWRAPPER_H
#define MATRIX4DWRAPPER_H

//extension by CAB
#include <vector>
#include <string>
#include <fstream>

#include <basics/utilities/UbException.h>
#include <MarchingCubes/McTypes.h>

//neu: matrix muss lediglich double operator()(int x1, int x2, int x3, int x4) überladen, dann kann man sie verwenden!!

//////////////////////////////////////////////////////////////////////////
//Matrix4DWrapper-Wrapper
//example:
//   int indexForDataValue = 1;
//   CbUniformMatrix4D<double> data(10,8,5,2);
//   for(int x3=0; x3<data.getNX3(); x3++)
//    for(int x2=0; x2<data.getNX2(); x2++)
//      for(int x1=0; x1<data.getNX1(); x1++)
//        data(x1,x2,x3,indexForDataValue) = x1;
//
//   Matrix4DWrapper< CbUniformMatrix4D<double> > wrapper(&data,indexForDataValue);
//   MarchingCubes< Matrix4DWrapper< CbUniformMatrix4D<double> > > mc( wrapper );
//
//   mc.init_all();
//   mc.run(3.5);
//   mc.writeUCD("c:/temp/triangles.inp");
//   mc.clean_all();

namespace McCubes{

template< typename Matrix4D >
class Matrix4DWrapper
{
public:
   Matrix4DWrapper() 
      :  valIndex(-1), matrix(NULL)
       , minX1(-1), minX2(-1), minX3(-1)
       , maxX1(-1), maxX2(-1), maxX3(-1)
       , nx1(-1)  , nx2(-1)  , nx3(-1)
   {
      //wird benötigt, damit MarchingCubes generell eine membervariabel erstellen kann
   }
   /*==========================================================*/
   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex,const int& n1, const int& nx2, const int& nx3)
      :  valIndex(valIndex), matrix(matrix)
       , nx1(nx1), nx2(nx2), nx3(nx3)
   {
      minX1 = minX2 = minX3 = 0;
      maxX1 = nx1-1;
      maxX2 = nx2-1;
      maxX3 = nx3-1;
   }
   /*==========================================================*/
   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex)
      : valIndex(valIndex), matrix(matrix)
   {
      nx1 = matrix->getNX1();
      nx2 = matrix->getNX2();
      nx3 = matrix->getNX3();

      minX1 = minX2 = minX3 = 0;
      maxX1 = nx1-1;
      maxX2 = nx2-1;
      maxX3 = nx3-1;
   }
   /*==========================================================*/
   //wenn man z.B. matrixX1 von[0..10] geht und man nur den bereich 1..9 fuer MC
   //verwenden möchte -> minX1=1 und maxX2=2
   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex, const int& minX1, const int& minX2, const int& minX3,
                                                           const int& maxX1, const int& maxX2, const int& maxX3)
       :  valIndex(valIndex), matrix(matrix)
        , minX1(minX1), minX2(minX2), minX3(minX3)
        , maxX1(maxX1), maxX2(maxX2), maxX3(maxX3)
   {
      nx1 = matrix->getNX1();
      nx2 = matrix->getNX2();
      nx3 = matrix->getNX3();

      if(minX1<0 || minX2<0 || minX3<0 || maxX1>=nx1 || maxX2>=nx2 || maxX3>=nx3)
         throw UbException(UB_EXARGS,"range error");
   }
   /*==========================================================*/
   //wenn man z.B. matrixX1 von[0..10] geht und man nur den bereich 1..9 fuer MC
   //verwenden möchte -> minX1=1 und maxX2=2
   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex, const int& minX1, const int& minX2, const int& minX3
                                                         , const int& maxX1, const int& maxX2, const int& maxX3
                                                         , const int& n1   , const int& nx2  , const int& nx3   )
      :  valIndex(valIndex), matrix(matrix)
       , minX1(minX1), minX2(minX2), minX3(minX3)
       , maxX1(maxX1), maxX2(maxX2), maxX3(maxX3)
       , nx1(nx1)    , nx2(nx2)    , nx3(nx3)
   {
      if(minX1<0 || minX2<0 || minX3<0 || maxX1>=nx1 || maxX2>=nx2 || maxX3>=nx3)
         throw UbException(UB_EXARGS,"range error");
   }
   /*==========================================================*/
   inline real getData(const int& x1, const int& x2, const int& x3 ) const
   {
      return static_cast<real>( (*matrix)(x1, x2, x3, valIndex) );
   }
   /*==========================================================*/
   inline int getMinX1() const { return minX1; }  
   inline int getMinX2() const { return minX2; }  
   inline int getMinX3() const { return minX3; }  

   inline int getMaxX1() const { return maxX1; }  
   inline int getMaxX2() const { return maxX2; }  
   inline int getMaxX3() const { return maxX3; }  

   inline int getNX1()   const { return nx1;   }  
   inline int getNX2()   const { return nx2;   }  
   inline int getNX3()   const { return nx3;   }  

protected:
   int valIndex;
   Matrix4D* matrix;
   int minX1, minX2, minX3, maxX1, maxX2, maxX3, nx1, nx2, nx3;
};

} //namespace McCubes

#endif //MATRIX4DWRAPPER_H
