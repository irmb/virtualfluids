//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Utilities Utilities
//! \ingroup cpu_core core
//! \{
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

//! \}
