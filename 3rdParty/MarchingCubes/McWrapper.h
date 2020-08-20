// _    ___      __              __________      _     __
//| |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
//| | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
//| |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
//|___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
//#ifndef MCWRAPPER_H
//#define MCWRAPPER_H
//
//extension by CAB
//#include <vector>
//#include <string>
//#include <fstream>
//
//#include <basics/utilities/UbException.h>
//#include <basics/container/CbUniformMatrix3D.h>
//#include <basics/container/CbUniformMatrix4D.h>
//
//#include <3rdParty/MarchingCubes/McTypes.h>
//
//namespace McCubes{
//
//////////////////////////////////////////////////////////////////////////
//Standard-Wrapper
//  MarchingCubes<DataWrapper> mc( 10,8,5 );
//  for(int z=0; z<mc.size_z(); z++)
//   for(int y=0; y<mc.size_y(); y++)
//     for(int x=0; x<mc.size_x(); x++)
//       mc.set_data(x,x,y,z);
//
//  //mc.set_method(false) //<-MC methode, STD=false
//  mc.init_all();
//
//  mc.run(3.5);
//  mc.writeUCD("c:/temp/triangles.inp");
//  mc.clean_all();

//template< typename T=real >
//class McDataWrapper
//{
//public:
//   typedef T                                                   value_type;
//   typedef typename std::vector< value_type >::reference       reference;
//   typedef typename std::vector< value_type >::const_reference const_reference;
//   typedef typename std::vector< value_type >::pointer         pointer;
//   typedef typename std::vector< value_type >::const_pointer   const_pointer;
//
//public:
//   McDataWrapper() : size_x1(-1), size_x2(-1), size_x3(-1)
//   {
//   }
//   /*==========================================================*/
//   McDataWrapper(const int& size_x1, const int& size_x2, const int& size_x3, const T& initVal=T())
//   {
//      this->resize(size_x1,size_x2,size_x3,initVal);
//   }
//   /*=======================================================================*/
//   reference operator() (const int& x1, const int& x2, const int& x3)
//   {
//      #ifdef _DEBUG
//         return this->data.at(x1 + x2*size_x1 + x3*size_x1*size_x2);
//      #else
//         return this->data[x1 + x2*size_x1 + x3*size_x1*size_x2];
//      #endif
//   }
//   /*=======================================================================*/
//   const_reference operator() (const int& x1, const int& x2, const int& x3)	const
//   {
//      #ifdef _DEBUG
//         return this->data.at(x1 + x2*size_x1 + x3*size_x1*size_x2);
//      #else
//         return this->data[x1 + x2*size_x1 + x3*size_x1*size_x2];
//      #endif
//   }
//   /*==========================================================*/
//   inline void setData( const T& val, const int& x1, const int& x2, const int& x3 )
//   {
//      #ifdef _DEBUG
//         this->data.at(x1 + x2*size_x1 + x3*size_x1*size_x2) = val;
//      #else
//         this->data[x1 + x2*size_x1 + x3*size_x1*size_x2] = val;
//      #endif
//   }
//   /*==========================================================*/
//   inline value_type getData(const int& x1, const int& x2, const int& x3 ) const
//   {
//      #ifdef _DEBUG
//         return this->data.at(x1 + x2*size_x1 + x3*size_x1*size_x2);
//      #else
//         return this->data[x1 + x2*size_x1 + x3*size_x1*size_x2];
//      #endif
//   }
//   /*==========================================================*/
//   inline void resize(const int& size_x1, const int& size_x2, const int& size_x3)
//   {
//      if(size_x1>0 && size_x2>0 && size_x3>0)
//      {
//         this->size_x1 = size_x1;
//         this->size_x2 = size_x2;
//         this->size_x3 = size_x3;
//         this->data.resize(size_x1*size_x2*size_x3);
//      }
//   }
//   /*==========================================================*/
//   inline void resize(const int& size_x1, const int& size_x2, const int& size_x3, const T& initVal)
//   {
//      if(size_x1>0 && size_x2>0 && size_x3>0)
//      {
//         this->size_x1 = size_x1;
//         this->size_x2 = size_x2;
//         this->size_x3 = size_x3;
//         this->data.resize(size_x1*size_x2*size_x3,initVal);
//      }
//   }
//   /*==========================================================*/
//   int getNX1() const { return size_x1; }
//   int getNX2() const { return size_x2; }
//   int getNX3() const { return size_x3; }
//
//protected:
//   std::vector<T> data;
//   int size_x1, size_x2, size_x3;
//};
//
//////////////////////////////////////////////////////////////////////////
//Matrix4DWrapper-Wrapper
//example:
//  int indexForDataValue = 1;
//  CbUniformMatrix4D<double> data(10,8,5,2);
//  for(int x3=0; x3<data.getNX3(); x3++)
//   for(int x2=0; x2<data.getNX2(); x2++)
//     for(int x1=0; x1<data.getNX1(); x1++)
//       data(x1,x2,x3,indexForDataValue) = x1;
//
//  Matrix4DWrapper< CbUniformMatrix4D<double> > wrapper(&data,indexForDataValue);
//  MarchingCubes< Matrix4DWrapper< CbUniformMatrix4D<double> > > mc( wrapper );
//
//  mc.init_all();
//  mc.run(3.5);
//  mc.writeUCD("c:/temp/triangles.inp");
//  mc.clean_all();
//template< typename Matrix4D >
//class Matrix4DWrapper
//{
//public:
//   typedef typename Matrix4D::value_type value_type;
//
//public:
//   Matrix4DWrapper() : valIndex(-1), matrix(NULL), minX1(-1), minX2(-1), minX3(-1), maxX1(-1), maxX2(-1), maxX3(-1)
//   {
//      //wird benoetigt, damit MarchingCubes generell eine membervariabel erstellen kann
//   }
//   /*==========================================================*/
//   //fuer beliebige matrizen
//   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex)
//      : valIndex(valIndex), matrix(matrix)
//   {
//
//   }
//   /*==========================================================*/
//   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex,const int& n1, const int& nx2, const int& nx3)
//      : valIndex(valIndex), matrix(matrix), nx1(nx1), nx2(nx2), nx3(nx3)
//   {
//      minX1 = minX2 = minX3 = 0;
//      maxX1 = nx1-1;
//      maxX2 = nx2-1;
//      maxX3 = nx3-1;
//   }
//   /*==========================================================*/
//   //wenn man z.B. matrixX1 von[0..10] geht und man nur den bereich 1..9 fuer MC
//   //verwenden moechte -> minX1=1 und maxX2=2
//   Matrix4DWrapper( Matrix4D* matrix, const int& valIndex, const int& minX1, const int& minX2, const int& minX3,
//                                                           const int& maxX1, const int& maxX2, const int& maxX3)
//       : valIndex(valIndex), matrix(matrix), minX1(minX1), minX2(minX2), minX3(minX3), maxX1(maxX1), maxX2(maxX2), maxX3(maxX3)
//   {
//      nx1 = matrix->getNX1()-1;
//      nx2 = matrix->getNX2()-1;
//      nx3 = matrix->getNX3()-1;
//
//      if(minX1<0 || minX2<0 || minX3<0 || maxX1>=nx1 || maxX2>=nx2 || maxX3>=nx3)
//         throw UbException(UB_EXARGS,"range error");
//   }
//   /*==========================================================*/
//   inline void setData( const real& val, const int& x1, const int& x2, const int& x3 )
//   {
//      (*matrix)(minX1+x1, minX2+x2, minX3+x3, valIndex) = static_cast<value_type>(val);
//   }
//   /*==========================================================*/
//   inline value_type getData(const int& x1, const int& x2, const int& x3 ) const
//   {
//      return (*matrix)(minX1+x1, minX2+x2, minX3+x3, valIndex);
//   }
//   /*==========================================================*/
//   inline void resize(const int& size_x1, const int& size_x2, const int& size_x3)
//   {
//      throw UbException("Matrix4DWrapper::resize(int,int,int) - mit diesem wrapper nicht erlaubt");
//   }
//   /*==========================================================*/
//   inline int getMinX1() const { return minX1; }  
//   inline int getMinX2() const { return minX2; }  
//   inline int getMinX3() const { return minX3; }  
//
//   inline int getMaxX1() const { return maxX1; }  
//   inline int getMaxX2() const { return maxX2; }  
//   inline int getMaxX3() const { return maxX3; }  
//
//   inline int getNX1()   const { nx1; }  
//   inline int getNX2()   const { nx2; }  
//   inline int getNX3()   const { nx3; }  
//
//protected:
//   int valIndex;
//   Matrix4D* matrix;
//   int minX1, minX2, minX3, maxX1, maxX2, maxX3, nx1, nx2, nx3;
//};
//
//////////////////////////////////////////////////////////////////////////
//Matrix3DWrapper-Wrapper
//  CbUniformMatrix3D<double> data(10,8,5);
//  for(int x3=0; x3<data.getNX3(); x3++)
//   for(int x2=0; x2<data.getNX2(); x2++)
//     for(int x1=0; x1<data.getNX1(); x1++)
//       data(x1,x2,x3) = x1;
//
//  Matrix3DWrapper< CbUniformMatrix3D<double> > wrapper(&data);
//  MarchingCubes< Matrix3DWrapper< CbUniformMatrix3D<double> > > mc( wrapper );
//
//  mc.init_all();
//  mc.run(3.5);
//  mc.writeUCD("c:/temp/triangles.inp");
//  mc.clean_all();
//template< typename Matrix3D >
//class Matrix3DWrapper
//{
//public:
//   typedef typename Matrix3D::value_type value_type;
//
//public:
//   Matrix3DWrapper() : matrix(NULL), minX1(-1), minX2(-1), minX3(-1), maxX1(-1), maxX2(-1), maxX3(-1)
//   {
//      //wird benoetigt, damit MarchingCubes generell eine membervariabel erstellen kann
//   }
//   /*==========================================================*/
//   Matrix3DWrapper( Matrix3D* matrix)
//      : matrix(matrix)
//   {
//      minX1 = minX2 = minX3 = 0;
//      maxX1 = matrix->getNX1();
//      maxX2 = matrix->getNX2();
//      maxX3 = matrix->getNX3();
//
//      minX1 = minX2 = minX3 = 0;
//      maxX1 = nx1-1;
//      maxX2 = nx2-1;
//      maxX3 = nx3-1;
//
//   }
//   /*==========================================================*/
//   Matrix3DWrapper( Matrix3D* matrix, const int& minX1, const int& minX2, const int& minX3
//                                    , const int& maxX1, const int& maxX2, const int& maxX3 )
//      : matrix(matrix), minX1(minX1), minX2(minX2), minX3(minX3), maxX1(maxX1), maxX2(maxX2), maxX3(maxX3)
//   {
//
//   }
//   /*==========================================================*/
//   Matrix3DWrapper(const int& size_x1, const int& size_x2, const int& size_x3)
//   {
//      throw UbException("Matrix3DWrapper(int,int,int) - mit diesem wrapper nicht erlaubt");
//   }
//   /*==========================================================*/
//   inline void setData( const real& val, const int& x1, const int& x2, const int& x3 )
//   {
//      (*matrix)(minX1+x1, minX2+x2, minX3+x3) = static_cast<value_type>(val);
//   }
//   /*==========================================================*/
//   inline value_type getData(const int& x1, const int& x2, const int& x3 ) const
//   {
//      return (*matrix)(minX1+x1, minX2+x2, minX3+x3);
//   }
//   /*==========================================================*/
//   inline void resize(const int& size_x1, const int& size_x2, const int& size_x3)
//   {
//      throw UbException("Matrix3DWrapper::resize(int,int,int) - mit diesem wrapper nicht erlaubt");
//   }
//   /*==========================================================*/
//   inline int getMinX1() const { return minX1; }  
//   inline int getMinX2() const { return minX2; }  
//   inline int getMinX3() const { return minX3; }  
//
//   inline int getMaxX1() const { return maxX1; }  
//   inline int getMaxX2() const { return maxX2; }  
//   inline int getMaxX3() const { return maxX3; }  
//
//   inline int getNX1()   const { nx1; }  
//   inline int getNX2()   const { nx2; }  
//   inline int getNX3()   const { nx3; }  
//
//protected:
//   Matrix3D* matrix;
//   int minX1, minX2, minX3, maxX1, maxX2, maxX3, nx1, nx2, nx3;
//};
//
//} //namespace McCubes
//
//#endif //MCWRAPPER_H
