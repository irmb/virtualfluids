#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <math.h>
#include "muParser.h"

namespace Utilities
{
   static bool isEven( int integer )
   {
      if ( integer % 2 == 0 )
         return true;
      else
         return false;
   }

   static bool isOdd( int integer )
   {
      if ( integer % 2 != 0 )
         return true;
      else
         return false;
   }

   //convert from real to int
   static int cint(real x)
   {
      real intpart;
      if (modf(x,&intpart)>=.5)
         return static_cast<int> (floor(x)+1);
      else
         return static_cast<int> (floor(x));
   }

   //create new mu parser for duct parabolic profile
   //inflow in X
   static mu::Parser getDuctParaboloidX(real Cy, real Hy, real Cz, real Hz, real V)
   {
      mu::Parser fct;
      fct.SetExpr("V*(((-(x2-Cy)^2.0+(Hy/2.0)^2.0)/(Hy/2.0)^2.0)*((-(x3-Cz)^2.0+(Hz/2.0)^2.0)/(Hz/2.0)^2.0))" );
      fct.DefineConst("Cy", Cy);
      fct.DefineConst("Hy", Hy);
      fct.DefineConst("Cz", Cz);
      fct.DefineConst("Hz", Hz);
      fct.DefineConst("V" , V );
      return fct;
   }
   //inflow in Y
   static mu::Parser getDuctParaboloidY(real Cx, real Hx, real Cz, real Hz, real V)
   {
      mu::Parser fct;
      fct.SetExpr("V*(((-(x1-Cx)^2.0+(Hx/2.0)^2.0)/(Hx/2.0)^2.0)*((-(x3-Cz)^2.0+(Hz/2.0)^2.0)/(Hz/2.0)^2.0))" );
      fct.DefineConst("Cx", Cx);
      fct.DefineConst("Hx", Hx);
      fct.DefineConst("Cz", Cz);
      fct.DefineConst("Hz", Hz);
      fct.DefineConst("V" , V );
      return fct;
   }
   //inflow in Z
   static mu::Parser getDuctParaboloidZ(real Cx, real Hx, real Cy, real Hy, real V)
   {
      mu::Parser fct;
      fct.SetExpr("V*(((-(x1-Cx)^2.0+(Hx/2.0)^2.0)/(Hx/2.0)^2.0)*((-(x2-Cy)^2.0+(Hy/2.0)^2.0)/(Hy/2.0)^2.0))" );
      fct.DefineConst("Cx", Cx);
      fct.DefineConst("Hx", Hx);
      fct.DefineConst("Cy", Cy);
      fct.DefineConst("Hy", Hy);
      fct.DefineConst("V" , V );
      return fct;
   }
   //hash function
   static unsigned int RSHash(const std::string& str)
   {
      unsigned int b    = 378551;
      unsigned int a    = 63689;
      unsigned int hash = 0;

      for(std::size_t i = 0; i < str.length(); i++)
      {
         hash = hash * a + str[i];
         a    = a * b;
      }

      return hash;
   }
   //linear interpolation
   static real linear_interpolation1D(real x0, real y0, real x1, real y1, real x)
   {
      real a = (y1 - y0) / (x1 - x0);
      real b = -a*x0 + y0;
      real y = a * x + b;
      return y;
   }
}

#endif
