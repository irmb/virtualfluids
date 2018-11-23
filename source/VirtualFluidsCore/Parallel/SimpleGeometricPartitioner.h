/**
* @file SimpleGeometricPartitioner.h
* @author Kostyantyn Kucher
* @date 06/06/2011
*
* @section DESCRIPTION
*
* This class make simple geometric partitioning.
*/

#ifndef SIMPLEGEOMETRICPARTITIONER_H 
#define SIMPLEGEOMETRICPARTITIONER_H

#include "basics/utilities/UbTuple.h"
#include "basics/utilities/UbException.h"
#include "MathUtil.hpp"

class SimpleGeometricPartitioner
{
public:
   static UbTupleInt3 createDimensions(const int& x, const int& y, const int& z, const int& numberOfProcess)
   {
      int xyz = x*y*z;

      int p = numberOfProcess;

      if (p == 1)
         return UbTupleInt3(1, 1, 1);

      double a = pow(p*pow(x,3.0)/xyz,1.0/3.0);
      double b = pow(p*pow(y,3.0)/xyz,1.0/3.0);
      double c = pow(p*pow(z,3.0)/xyz,1.0/3.0);

      MaxDim maxDim;
 
      if(c >= a && c >= b)
         maxDim = cDim;
      if(b >= a && b >= c)
         maxDim = bDim;
      if(a >= b && a >= c)
         maxDim = aDim;

      int dim1, dim2, dim3;
      dim1 = (int)Utilities::cint(a);
      dim2 = (int)Utilities::cint(b);
      dim3 = (int)Utilities::cint(c);
      if(dim1 <= 0) dim1 = 1;
      if(dim2 <= 0) dim2 = 1;
      if(dim3 <= 0) dim3 = 1;

      switch (maxDim)
      {
      case aDim: 
         dim1 = p/(dim2*dim3);
         if (dim1*dim2*dim3 != p)
         {
            dim2 = 1;
            dim3 = 1;
            dim1 = p;
         }
         break;
      case bDim: 
         dim2 = p/(dim1*dim3);
         if (dim1*dim2*dim3 != p)
         {
            dim1 = 1;
            dim3 = 1;
            dim2 = p;
         }
         break;
      case cDim: 
         dim3 = p/(dim1*dim2);
         if (dim1*dim2*dim3 != p)
         {
            dim1 = 1;
            dim2 = 1;
            dim3 = p;
         }
         break;
      }

      if (dim1>x || dim2>y || dim3>z)
      {
         UB_THROW( UbException(UB_EXARGS,"SimpleGeometricPartitioner::createDimensions: Segmentation fault - bad number of prozess") );
      }

      UbTupleInt3 dims(dim1, dim2, dim3);
      
      return dims;
   }

protected:
private:
   enum MaxDim {aDim,bDim,cDim};
};

#endif
