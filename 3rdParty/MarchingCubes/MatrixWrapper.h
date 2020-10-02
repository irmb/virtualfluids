//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MATRIXWRAPPER_H
#define MATRIXWRAPPER_H

//extension by CAB
#include <vector>
#include <string>
#include <fstream>

#include <basics/utilities/UbException.h>

#include <MarchingCubes/McTypes.h>

//////////////////////////////////////////////////////////////////////////
//Standard-Wrapper
//   MarchingCubes<DataWrapper> mc( 10,8,5 );
//   for(int z=0; z<mc.size_z(); z++)
//    for(int y=0; y<mc.size_y(); y++)
//      for(int x=0; x<mc.size_x(); x++)
//        mc.set_data(x,x,y,z);
//
//   mc.init_all();
//   mc.run(3.5);
//   mc.writeUCD("c:/temp/triangles.inp");
//   mc.clean_all();

namespace McCubes{

template< typename T=real >
class MatrixWrapper
{
public:
   typedef T                                                   value_type;
   typedef typename std::vector< value_type >::reference       reference;
   typedef typename std::vector< value_type >::const_reference const_reference;
   typedef typename std::vector< value_type >::pointer         pointer;
   typedef typename std::vector< value_type >::const_pointer   const_pointer;

public:
   MatrixWrapper() : nx1(-1), nx2(-1), nx3(-1)
   {
   }
   /*==========================================================*/
   MatrixWrapper(const int& nx1, const int& nx2, const int& nx3, const T& initVal=T())
   {
      this->resize(nx1,nx2,nx3,initVal);
   }
   /*=======================================================================*/
   reference operator() (const int& x1, const int& x2, const int& x3)
   {
      #ifdef _DEBUG
         return this->data.at(x1 + x2*nx1 + x3*nx1*nx2);
      #else
         return this->data[x1 + x2*nx1 + x3*nx1*nx2];
      #endif
   }
   /*=======================================================================*/
   const_reference operator() (const int& x1, const int& x2, const int& x3)	const
   {
      #ifdef _DEBUG
         return this->data.at(x1 + x2*nx1 + x3*nx1*nx2);
      #else
         return this->data[x1 + x2*nx1 + x3*nx1*nx2];
      #endif
   }
   /*==========================================================*/
   inline void setData( const T& val, const int& x1, const int& x2, const int& x3 )
   {
      #ifdef _DEBUG
         this->data.at(x1 + x2*nx1 + x3*nx1*nx2) = val;
      #else
         this->data[x1 + x2*nx1 + x3*nx1*nx2] = val;
      #endif
   }
   /*==========================================================*/
   inline value_type getData(const int& x1, const int& x2, const int& x3 ) const
   {
      #ifdef _DEBUG
         return this->data.at(x1 + x2*nx1 + x3*nx1*nx2);
      #else
         return this->data[x1 + x2*nx1 + x3*nx1*nx2];
      #endif
   }
   /*==========================================================*/
   inline void resize(const int& nx1, const int& nx2, const int& nx3)
   {
      if(nx1>0 && nx2>0 && nx3>0)
      {
         this->nx1 = nx1;
         this->nx2 = nx2;
         this->nx3 = nx3;
         this->data.resize(nx1*nx2*nx3);
      }
   }
   /*==========================================================*/
   inline void resize(const int& nx1, const int& nx2, const int& nx3, const T& initVal)
   {
      if(nx1>0 && nx2>0 && nx3>0)
      {
         this->nx1 = nx1;
         this->nx2 = nx2;
         this->nx3 = nx3;
         this->data.resize(nx1*nx2*nx3,initVal);
      }
   }
   /*==========================================================*/
   inline int getMinX1() const { return 0;     }  
   inline int getMinX2() const { return 0;     }  
   inline int getMinX3() const { return 0;     }  

   inline int getMaxX1() const { return nx1-1; }  
   inline int getMaxX2() const { return nx2-1; }  
   inline int getMaxX3() const { return nx3-1; }  

   inline int getNX1()   const { return nx1;   }  
   inline int getNX2()   const { return nx2;   }  
   inline int getNX3()   const { return nx3;   }  

protected:
   std::vector<T> data;
   int nx1, nx2, nx3;
};

} //namespace McCubes

#endif //MATRIXWRAPPER_H
