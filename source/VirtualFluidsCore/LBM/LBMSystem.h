#ifndef LBMSYSTEM_H
#define LBMSYSTEM_H

#include <cmath>
#include <string>
#include <iostream>

#ifdef RCF_USE_SF_SERIALIZATION
#include <SF/Serializer.hpp>

#if CAB_RCF <= 903
#include <SF/SerializeEnum.hpp>   
#endif
#endif //RCF_USE_SF_SERIALIZATION

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbTuple.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbSystem.h>

/*=========================================================================*/
/*  LBMSystem                                                            */
/*                                                                         */
/**
namespace for global system-functions
<BR><BR>
@author <A HREF="mailto:geller@irmb.tu-bs.de">S. Geller</A>
@version 1.0 - 07.01.11
*/ 

/*
usage: ...
*/

namespace LBMSystem
{
#ifndef SWIG
   using namespace UbMath;
#endif

//#define SINGLEPRECISION

#ifdef SINGLEPRECISION
   typedef float real;
   #define REAL_CAST(x) ( (LBMSystem::real)(x) )
#else
   typedef double real;
   #define REAL_CAST(x) ( x )
#endif

   extern real SMAG_CONST;

   //////////////////////////////////////////////////////////////////////////
   //!get LBM deltaT is equal LBM DeltaX
   //!deltaT is dependent from grid level 
   //!for first grid level is deltaT = 1.0
   //!for next grid level 1/2 etc.
   static real getDeltaT(int level)
   {
      return REAL_CAST(1.0/REAL_CAST(1<<level)); 
   }

   //////////////////////////////////////////////////////////////////////////
   //!calculate collision factor omega = 1.0/(3.0*viscosity/deltaT+0.5)
   //!deltaT is dependent from grid level 
   //!for first grid level is deltaT = 1.0
   //!for next grid level 1/2 etc.
   static real calcCollisionFactor(real viscosity, int level)
   {
      //return REAL_CAST(1.0/(3.0*viscosity/deltaT+0.5));
      return REAL_CAST(1.0/(3.0*viscosity/(1.0/REAL_CAST(1<<level))+0.5));
   }
}

//some typedefs for global namespace
typedef LBMSystem::real LBMReal;

//#define LBMSystem::real LBMReal



#endif

