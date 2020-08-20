//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MATRIX3DWRAPPER_H
#define MATRIX3DWRAPPER_H

//extension by CAB
#include <vector>
#include <string>
#include <fstream>

#include <basics/utilities/UbException.h>
#include <MarchingCubes/McTypes.h>

//neu: matrix muss lediglich double operator()(int x1, int x2, int x3) ueberladen, dann kann man sie verwenden!!

//////////////////////////////////////////////////////////////////////////
//Matrix3DWrapper-Wrapper
//   CbUniformMatrix3D<double> data(10,8,5);
//   for(int x3=0; x3<data.getNX3(); x3++)
//    for(int x2=0; x2<data.getNX2(); x2++)
//      for(int x1=0; x1<data.getNX1(); x1++)
//        data(x1,x2,x3) = x1;
//
//   Matrix3DWrapper< CbUniformMatrix3D<double> > wrapper(&data);
//   MarchingCubes< Matrix3DWrapper< CbUniformMatrix3D<double> > > mc( wrapper );
//
//   mc.init_all();
//   mc.run(3.5);
//   mc.writeUCD("c:/temp/triangles.inp");
//   mc.clean_all();

namespace McCubes{

template< typename Matrix3D >
class Matrix3DWrapper
{
public:
   Matrix3DWrapper()
      :  matrix(NULL)
       , minX1(-1), minX2(-1), minX3(-1)
       , maxX1(-1), maxX2(-1), maxX3(-1)
       , nx1(-1)  , nx2(-1)  , nx3(-1)
   {
      //wird benoetigt, damit MarchingCubes generell eine membervariabel erstellen kann
   }
   /*==========================================================*/
   Matrix3DWrapper( Matrix3D* matrix)
      : matrix(matrix)
   {
      nx1 = (int)matrix->getNX1();
      nx2 = (int)matrix->getNX2();
      nx3 = (int)matrix->getNX3();

      minX1 = minX2 = minX3 = 0;
      maxX1 = nx1-1;
      maxX2 = nx2-1;
      maxX3 = nx3-1;
   }
   /*==========================================================*/
   Matrix3DWrapper( Matrix3D* matrix, const int& n1, const int& nx2, const int& nx3)
      : matrix(matrix)
      , nx1(nx1), nx2(nx2), nx3(nx3)
   {
      minX1 = minX2 = minX3 = 0;
      maxX1 = nx1-1;
      maxX2 = nx2-1;
      maxX3 = nx3-1;
   }
   /*==========================================================*/
   Matrix3DWrapper( Matrix3D* matrix, const int& minX1, const int& minX2, const int& minX3
                                    , const int& maxX1, const int& maxX2, const int& maxX3 )
      :  matrix(matrix)
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
   //verwenden moechte -> minX1=1 und maxX2=2
   Matrix3DWrapper( Matrix3D* matrix, const int& minX1, const int& minX2, const int& minX3
                                    , const int& maxX1, const int& maxX2, const int& maxX3
                                    , const int& n1   , const int& nx2  , const int& nx3   )
      :   matrix(matrix)
        , minX1(minX1), minX2(minX2), minX3(minX3)
        , maxX1(maxX1), maxX2(maxX2), maxX3(maxX3)
        , nx1(n1)     , nx2(nx2)    , nx3(nx3)
   {
      if(minX1<0 || minX2<0 || minX3<0 || maxX1>=nx1 || maxX2>=nx2 || maxX3>=nx3)
         throw UbException(UB_EXARGS,"range error");
   }
   /*==========================================================*/
   inline real getData(const int& x1, const int& x2, const int& x3 ) const
   {
      return static_cast<real>( (*matrix)(x1, x2, x3) );
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
   Matrix3D* matrix;
   int minX1, minX2, minX3, maxX1, maxX2, maxX3, nx1, nx2, nx3;
};

} //namespace McCubes

#endif //MATRIX3DWRAPPER_H
