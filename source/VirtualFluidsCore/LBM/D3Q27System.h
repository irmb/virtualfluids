#ifndef D3Q27SYSTEM_H
#define D3Q27SYSTEM_H

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
//#include "Patch3DSystem.h"
#include "LBMSystem.h"
	
/*=========================================================================*/
/*  D3Q27System                                                            */
/*                                                                         */
/**
class for global system-functions
<BR><BR>
@author <A HREF="mailto:kucher@irmb.tu-bs.de">K. Kucher</A>
@version 1.0 - 22.10.09
*/ 

/*
usage: ...
*/


#ifndef SWIG
   using namespace UbMath;
#endif

namespace D3Q27System
{
   //enum COLLISIONMODEL { UNDEFINED, INCOMPLBGKMODEL,   COMPLBGKMODEL,   COMPLBGKWTMODEL,   INCOMPLBGKLESMODEL, INCOMPLBGKNONNEWTONIANMODEL    
   //                               , INCOMPGLBEJTMODEL, COMPGLBEJTMODEL, COMPGLBEJTWTMODEL, INCOMPGLBEJTLESMODEL, INCOMPGLBEJTWALEMODEL  
   //                               , CASCADEDMODEL};
   //
   // #if defined(RCF_USE_SF_SERIALIZATION) && (CAB_RCF <= 903)
   //    SF_SERIALIZE_ENUM(COLLISIONMODEL) //muss im namespace stehen, sonst funzt es nicht!
   // #endif

   ///*=====================================================================*/
   //std::string toString(const COLLISIONMODEL& model);
   ///*=====================================================================*/
   //COLLISIONMODEL getCollModelByString(const std::string& str);
   ///*=====================================================================*/

   ///*=====================================================================*/
   //static bool isCompModel(const COLLISIONMODEL& model) 
   //{
   //   switch(model)
   //   {
   //   case COMPLBGKMODEL               : return true; 
   //   case COMPLBGKWTMODEL             : return true; 
   //   case COMPGLBEJTWTMODEL           : return true;
   //   case COMPGLBEJTMODEL             : return true;
   //   case CASCADEDMODEL               : return true;
   //   
   //   case INCOMPLBGKMODEL             : return false;
   //   case INCOMPGLBEJTMODEL           : return false;
   //   case INCOMPLBGKLESMODEL          : return false;
   //   case INCOMPGLBEJTLESMODEL        : return false;
   //   case INCOMPGLBEJTWALEMODEL       : return false;
   //   case INCOMPLBGKNONNEWTONIANMODEL : return false;

   //   default: throw UbException(UB_EXARGS,"unknown model");
   //   }
   //}
   ///*=====================================================================*/
   //static bool isGLBEModel(const COLLISIONMODEL& model) 
   //{
   //   switch(model)
   //   {
   //   case COMPGLBEJTWTMODEL           : return true;
   //   case COMPGLBEJTMODEL             : return true;
   //   case INCOMPGLBEJTMODEL           : return true;
   //   case INCOMPGLBEJTLESMODEL        : return true;
   //   case INCOMPGLBEJTWALEMODEL       : return false;

   //   case COMPLBGKMODEL               : return false; 
   //   case COMPLBGKWTMODEL             : return false; 
   //   case INCOMPLBGKMODEL             : return false;
   //   case INCOMPLBGKLESMODEL          : return false;
   //   case INCOMPLBGKNONNEWTONIANMODEL : return false;

   //   default: throw UbException(UB_EXARGS,"unknown model");
   //   }
   //}
   //static bool isLESModel(const COLLISIONMODEL& model) 
   //{
   //   switch(model)
   //   {
   //   case INCOMPGLBEJTLESMODEL        : return true;
   //   case INCOMPLBGKLESMODEL          : return true;
   //   case INCOMPGLBEJTWALEMODEL       : return true;
   //   
   //   case COMPGLBEJTWTMODEL           : return false;
   //   case COMPGLBEJTMODEL             : return false;
   //   case INCOMPGLBEJTMODEL           : return false;
   //   case COMPLBGKMODEL               : return false; 
   //   case COMPLBGKWTMODEL             : return false; 
   //   case INCOMPLBGKMODEL             : return false;
   //   case INCOMPLBGKNONNEWTONIANMODEL : return false;

   //   default: throw UbException(UB_EXARGS,"unknown model");
   //   }
   //}

   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   //DIRECTION STUFF
   static const int FSTARTDIR = 0;
   static const int FENDDIR   = 25;   //D3Q27

   //static const int FSTARTDIR = 1;
   //static const int FENDDIR   = 26;   //D3Q27

   static const int STARTF = 0;
   static const int ENDF   = 26;   //D3Q27

   //extern const int EX1[ENDF+1];
   //extern const int EX2[ENDF+1];
   //extern const int EX3[ENDF+1];

   static const int STARTDIR = 0;
   static const int ENDDIR   = 26; //alle geometrischen richtungen

   extern const int DX1[ENDDIR+1];
   extern const int DX2[ENDDIR+1];
   extern const int DX3[ENDDIR+1];
   extern const double WEIGTH[ENDDIR+1];

   extern const double cNorm[3][ENDDIR];
   
   //static const int ZERO /*f0 */ = 0;
   //static const int E    /*f1 */ = 1;
   //static const int W    /*f2 */ = 2;
   //static const int N    /*f3 */ = 3;
   //static const int S    /*f4 */ = 4;
   //static const int T    /*f5 */ = 5;
   //static const int B    /*f6 */ = 6;
   //static const int NE   /*f7 */ = 7;
   //static const int SW   /*f8 */ = 8;
   //static const int SE   /*f9 */ = 9;
   //static const int NW   /*f10*/ = 10;
   //static const int TE   /*f11*/ = 11;
   //static const int BW   /*f12*/ = 12;
   //static const int BE   /*f13*/ = 13;
   //static const int TW   /*f14*/ = 14;
   //static const int TN   /*f15*/ = 15;
   //static const int BS   /*f16*/ = 16;
   //static const int BN   /*f17*/ = 17;
   //static const int TS   /*f18*/ = 18;
   //static const int TNE          = 19;
   //static const int TNW          = 20;
   //static const int TSE          = 21;
   //static const int TSW          = 22;
   //static const int BNE          = 23;
   //static const int BNW          = 24;
   //static const int BSE          = 25;
   //static const int BSW          = 26;

   static const int E    /*f1 */ = 0;
   static const int W    /*f2 */ = 1;
   static const int N    /*f3 */ = 2;
   static const int S    /*f4 */ = 3;
   static const int T    /*f5 */ = 4;
   static const int B    /*f6 */ = 5;
   static const int NE   /*f7 */ = 6;
   static const int SW   /*f8 */ = 7;
   static const int SE   /*f9 */ = 8;
   static const int NW   /*f10*/ = 9;
   static const int TE   /*f11*/ = 10;
   static const int BW   /*f12*/ = 11;
   static const int BE   /*f13*/ = 12;
   static const int TW   /*f14*/ = 13;
   static const int TN   /*f15*/ = 14;
   static const int BS   /*f16*/ = 15;
   static const int BN   /*f17*/ = 16;
   static const int TS   /*f18*/ = 17;
   static const int TNE          = 18;
   static const int TNW          = 19;
   static const int TSE          = 20;
   static const int TSW          = 21;
   static const int BNE          = 22;
   static const int BNW          = 23;
   static const int BSE          = 24;
   static const int BSW          = 25;
   static const int ZERO /*f0 */ = 26;

   static const int INV_E   = W;  
   static const int INV_W   = E;  
   static const int INV_N   = S;  
   static const int INV_S   = N;  
   static const int INV_T   = B;  
   static const int INV_B   = T;  
   static const int INV_NE  = SW; 
   static const int INV_SW  = NE; 
   static const int INV_SE  = NW; 
   static const int INV_NW  = SE; 
   static const int INV_TE  = BW; 
   static const int INV_BW  = TE; 
   static const int INV_BE  = TW; 
   static const int INV_TW  = BE; 
   static const int INV_TN  = BS; 
   static const int INV_BS  = TN; 
   static const int INV_BN  = TS; 
   static const int INV_TS  = BN; 
   static const int INV_TNE = BSW;
   static const int INV_TNW = BSE;
   static const int INV_TSE = BNW;
   static const int INV_TSW = BNE;
   static const int INV_BNE = TSW;
   static const int INV_BNW = TSE;
   static const int INV_BSE = TNW;
   static const int INV_BSW = TNE;
                                       
   extern const int INVDIR[ENDDIR+1];

   static const int ET_E   = 0;
   static const int ET_W   = 0;
   static const int ET_N   = 1;
   static const int ET_S   = 1;
   static const int ET_T   = 2;
   static const int ET_B   = 2;
   static const int ET_NE  = 3;
   static const int ET_SW  = 3;
   static const int ET_SE  = 4;
   static const int ET_NW  = 4;
   static const int ET_TE  = 5;
   static const int ET_BW  = 5;
   static const int ET_BE  = 6;
   static const int ET_TW  = 6;
   static const int ET_TN  = 7;
   static const int ET_BS  = 7;
   static const int ET_BN  = 8;
   static const int ET_TS  = 8;
   static const int ET_TNE = 9;
   static const int ET_BSW = 9;
   static const int ET_TNW = 10;
   static const int ET_BSE = 10;
   static const int ET_TSE = 11;
   static const int ET_BNW = 11;
   static const int ET_TSW = 12;
   static const int ET_BNE = 12;

   static const int M_RHO     = 0;  
   static const int M_EN      = 1;  
   static const int M_EPS     = 2;  
   static const int M_JX1     = 3;  
   static const int M_QX1     = 4;  
   static const int M_JX2     = 5;  
   static const int M_QX2     = 6;  
   static const int M_JX3     = 7;  
   static const int M_QX3     = 8;  
   static const int M_3PX1X1  = 9;  
   static const int M_3PIX1X1 = 10; 
   static const int M_PWW     = 11; 
   static const int M_PIWW    = 12; 
   static const int M_PX1X2   = 13; 
   static const int M_PX2X3   = 14; 
   static const int M_PX1X3   = 15; 
   static const int M_MX1     = 16; 
   static const int M_MX2     = 17; 
   static const int M_MX3     = 18; 
   
   static const int STARTM = 0;
   static const int ENDM   = 18;   //D3Q27


   
   //////////////////////////////////////////////////////////////////////////
   //MACROSCOPIC VALUES                  
   /*=====================================================================*/
   static LBMReal getDensity(const LBMReal* const& f/*[27]*/)
   {
      return  ((f[TNE] + f[BSW])+(f[TSE]+f[BNW]))+((f[BSE]+f[TNW])+ (f[TSW]+f[BNE]))
             +(((f[NE] + f[SW]) + (f[SE] + f[NW]))+((f[TE] + f[BW])+(f[BE]+ f[TW]))
             +((f[BN] + f[TS]) + (f[TN] + f[BS])))+((f[E] + f[W])+(f[N] + f[S])
             +(f[T] + f[B]))+f[ZERO];
   }
   /*=====================================================================*/
   //ACHTUNG: gilt nicht fuer alle modelle -> praedikat verwenden anstelle static! toDo
   static LBMReal getPressure(const LBMReal* const& f/*[27]*/)
   {
      return  REAL_CAST( c1o3 )*getDensity(f);
   }
   /*=====================================================================*/
   static LBMReal getIncompVelocityX1(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[BSE]-f[TNW]) + (f[BNE]-f[TSW]))) +
             (((f[BE]-f[TW]) + (f[TE]-f[BW])) + ((f[SE]-f[NW]) + (f[NE]-f[SW]))) +
             (f[E]-f[W]));  
   }
   /*=====================================================================*/
   static LBMReal getIncompVelocityX2(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE]-f[BSW]) + (f[BNW]-f[TSE])) + ((f[TNW]-f[BSE]) + (f[BNE]-f[TSW]))) +
             (((f[BN]-f[TS]) + (f[TN]-f[BS])) + ((f[NW]-f[SE]) + (f[NE]-f[SW]))) +
             (f[N]-f[S]));  
   }
   /*=====================================================================*/
   static LBMReal getIncompVelocityX3(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
             (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) +
             (f[T] - f[B]));
   }
   /*=====================================================================*/
   static void calcDensity(const LBMReal* const& f/*[27]*/, LBMReal& rho)
   {
      rho = ((f[TNE] + f[BSW])+(f[TSE]+f[BNW]))+((f[BSE]+f[TNW])+ (f[TSW]+f[BNE]))
         +(((f[NE] + f[SW]) + (f[SE] + f[NW]))+((f[TE] + f[BW])+(f[BE]+ f[TW]))
         +((f[BN] + f[TS]) + (f[TN] + f[BS])))+((f[E] + f[W])+(f[N] + f[S])
         +(f[T] + f[B]))+f[ZERO];
         
   }
   /*=====================================================================*/
   static void calcIncompVelocityX1(const LBMReal* const& f/*[27]*/, LBMReal& vx1)
   {
      vx1 = ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[BSE]-f[TNW]) + (f[BNE]-f[TSW]))) +
             (((f[BE]-f[TW]) + (f[TE]-f[BW])) + ((f[SE]-f[NW]) + (f[NE]-f[SW]))) +
             (f[E]-f[W]));
   }
   /*=====================================================================*/
   static void calcIncompVelocityX2(const LBMReal* const& f/*[27]*/, LBMReal& vx2)
   {
      vx2 = ((((f[TNE]-f[BSW]) + (f[BNW]-f[TSE])) + ((f[TNW]-f[BSE]) + (f[BNE]-f[TSW]))) +
             (((f[BN]-f[TS]) + (f[TN]-f[BS])) + ((f[NW]-f[SE]) + (f[NE]-f[SW]))) +
             (f[N]-f[S]));
   }
   /*=====================================================================*/
   static void calcIncompVelocityX3(const LBMReal* const& f/*[27]*/, LBMReal& vx3)
   {
      vx3 =((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[TNW]-f[BSE]) + (f[TSW]-f[BNE]))) +
             (((f[TS]-f[BN]) + (f[TN]-f[BS])) + ((f[TW]-f[BE]) + (f[TE]-f[BW]))) +
             (f[T]-f[B]));
   }
   /*=====================================================================*/
   static LBMReal getCompVelocityX1(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[BSE]-f[TNW]) + (f[BNE]-f[TSW]))) +
             (((f[BE]-f[TW]) + (f[TE]-f[BW])) + ((f[SE]-f[NW]) + (f[NE]-f[SW]))) +
             (f[E]-f[W]))/getDensity(f);  
   }
   /*=====================================================================*/
   static LBMReal getCompVelocityX2(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE]-f[BSW]) + (f[BNW]-f[TSE])) + ((f[TNW]-f[BSE]) + (f[BNE]-f[TSW]))) +
             (((f[BN]-f[TS]) + (f[TN]-f[BS])) + ((f[NW]-f[SE]) + (f[NE]-f[SW]))) +
             (f[N]-f[S]))/getDensity(f);  
  }
   /*=====================================================================*/
   static LBMReal getCompVelocityX3(const LBMReal* const& f/*[27]*/)
   {
      return ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[TNW]-f[BSE]) + (f[TSW]-f[BNE]))) +
             (((f[TS]-f[BN]) + (f[TN]-f[BS])) + ((f[TW]-f[BE]) + (f[TE]-f[BW]))) +
             (f[T]-f[B]))/getDensity(f);
   }
   /*=====================================================================*/
   static void calcCompVelocityX1(const LBMReal* const& f/*[27]*/, LBMReal& vx1)
   {
      vx1 = ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[BSE]-f[TNW]) + (f[BNE]-f[TSW]))) +
            (((f[BE]-f[TW]) + (f[TE]-f[BW])) + ((f[SE]-f[NW]) + (f[NE]-f[SW]))) +
            (f[E]-f[W]))/getDensity(f);  
   }
   /*=====================================================================*/
   static void calcCompVelocityX2(const LBMReal* const& f/*[27]*/, LBMReal& vx2)
   {
      vx2 = ((((f[TNE]-f[BSW]) + (f[BNW]-f[TSE])) + ((f[TNW]-f[BSE]) + (f[BNE]-f[TSW]))) +
            (((f[BN]-f[TS]) + (f[TN]-f[BS])) + ((f[NW]-f[SE]) + (f[NE]-f[SW]))) +
            (f[N]-f[S]))/getDensity(f);  
   }
   /*=====================================================================*/
   static void calcCompVelocityX3(const LBMReal* const& f/*[27]*/, LBMReal& vx3)
   {
      vx3 = ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[TNW]-f[BSE]) + (f[TSW]-f[BNE]))) +
            (((f[TS]-f[BN]) + (f[TN]-f[BS])) + ((f[TW]-f[BE]) + (f[TE]-f[BW]))) +
            (f[T]-f[B]))/getDensity(f);
   }
   /*=====================================================================*/
   static void calcIncompMacroscopicValues(const LBMReal* const& f/*[27]*/, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3)
   {
      D3Q27System::calcDensity(f, rho);
      D3Q27System::calcIncompVelocityX1(f, vx1);
      D3Q27System::calcIncompVelocityX2(f, vx2);
      D3Q27System::calcIncompVelocityX3(f, vx3);
   }

   /*=====================================================================*/
   static void calcCompMacroscopicValues(const LBMReal* const& f/*[27]*/, LBMReal& drho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3)
   {
      D3Q27System::calcDensity(f, drho);
      D3Q27System::calcIncompVelocityX1(f, vx1);
      D3Q27System::calcIncompVelocityX2(f, vx2);
      D3Q27System::calcIncompVelocityX3(f, vx3);
      LBMReal rho = drho+one;
      vx1/=rho;
      vx2/=rho;
      vx3/=rho;
   }
   //////////////////////////////////////////////////////////////////////////
   static LBMReal getCompFeqForDirection(const int& direction, const LBMReal& drho,const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)
   {
      LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

      //switch(direction)    
      //{
      //   case ZERO : return REAL_CAST( c8o27*rho*(1.0-cu_sq));
      //   case E : return REAL_CAST(  c2o27*rho*(1.0+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq));
      //   case W : return REAL_CAST(  c2o27*rho*(1.0+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq));
      //   case N : return REAL_CAST(  c2o27*rho*(1.0+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq));
      //   case S : return REAL_CAST(  c2o27*rho*(1.0+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq));
      //   case T : return REAL_CAST(  c2o27*rho*(1.0+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq));
      //   case B : return REAL_CAST(  c2o27*rho*(1.0+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq));
      //   case NE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq));
      //   case SW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq));
      //   case SE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq));
      //   case NW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq));
      //   case TE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq));
      //   case BW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq));
      //   case BE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq));
      //   case TW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq));
      //   case TN : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq));
      //   case BS : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq));
      //   case BN : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq));
      //   case TS : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq));
      //   case TNE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      //   case BSW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      //   case BNE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      //   case TSW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      //   case TSE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      //   case BNW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      //   case BSE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      //   case TNW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      //   default: throw UbException(UB_EXARGS,"unknown dir");
      //}


      ////-----
      LBMReal rho = drho+one;
      switch (direction)
      {
      case ZERO: return REAL_CAST(c8o27*(drho+rho*(-cu_sq)));
      case E: return REAL_CAST(c2o27*(drho+rho*(3.0*(vx1)+c9o2*(vx1)*(vx1)-cu_sq)));
      case W: return REAL_CAST(c2o27*(drho+rho*(3.0*(-vx1)+c9o2*(-vx1)*(-vx1)-cu_sq)));
      case N: return REAL_CAST(c2o27*(drho+rho*(3.0*(vx2)+c9o2*(vx2)*(vx2)-cu_sq)));
      case S: return REAL_CAST(c2o27*(drho+rho*(3.0*(-vx2)+c9o2*(-vx2)*(-vx2)-cu_sq)));
      case T: return REAL_CAST(c2o27*(drho+rho*(3.0*(vx3)+c9o2*(vx3)*(vx3)-cu_sq)));
      case B: return REAL_CAST(c2o27*(drho+rho*(3.0*(-vx3)+c9o2*(-vx3)*(-vx3)-cu_sq)));
      case NE: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx1+vx2)+c9o2*(vx1+vx2)*(vx1+vx2)-cu_sq)));
      case SW: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq)));
      case SE: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx1-vx2)+c9o2*(vx1-vx2)*(vx1-vx2)-cu_sq)));
      case NW: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq)));
      case TE: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx1+vx3)+c9o2*(vx1+vx3)*(vx1+vx3)-cu_sq)));
      case BW: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq)));
      case BE: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx1-vx3)+c9o2*(vx1-vx3)*(vx1-vx3)-cu_sq)));
      case TW: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq)));
      case TN: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx2+vx3)+c9o2*(vx2+vx3)*(vx2+vx3)-cu_sq)));
      case BS: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq)));
      case BN: return REAL_CAST(c1o54*(drho+rho*(3.0*(vx2-vx3)+c9o2*(vx2-vx3)*(vx2-vx3)-cu_sq)));
      case TS: return REAL_CAST(c1o54*(drho+rho*(3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq)));
      case TNE: return REAL_CAST(c1o216*(drho+rho*(3.0*(vx1+vx2+vx3)+c9o2*(vx1+vx2+vx3)*(vx1+vx2+vx3)-cu_sq)));
      case BSW: return REAL_CAST(c1o216*(drho+rho*(3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq)));
      case BNE: return REAL_CAST(c1o216*(drho+rho*(3.0*(vx1+vx2-vx3)+c9o2*(vx1+vx2-vx3)*(vx1+vx2-vx3)-cu_sq)));
      case TSW: return REAL_CAST(c1o216*(drho+rho*(3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq)));
      case TSE: return REAL_CAST(c1o216*(drho+rho*(3.0*(vx1-vx2+vx3)+c9o2*(vx1-vx2+vx3)*(vx1-vx2+vx3)-cu_sq)));
      case BNW: return REAL_CAST(c1o216*(drho+rho*(3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq)));
      case BSE: return REAL_CAST(c1o216*(drho+rho*(3.0*(vx1-vx2-vx3)+c9o2*(vx1-vx2-vx3)*(vx1-vx2-vx3)-cu_sq)));
      case TNW: return REAL_CAST(c1o216*(drho+rho*(3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq)));
      default: throw UbException(UB_EXARGS, "unknown dir");
      }

   }
   //////////////////////////////////////////////////////////////////////////
   static void calcCompFeq(LBMReal* const& feq/*[27]*/,const LBMReal& drho,const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)	
   {
      //LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

      //feq[ZERO] =  c8o27*rho*(1.0-cu_sq);
      //feq[E] =   c2o27*rho*(1.0+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq);
      //feq[W] =   c2o27*rho*(1.0+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq);
      //feq[N] =   c2o27*rho*(1.0+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq);
      //feq[S] =   c2o27*rho*(1.0+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq);
      //feq[T] =   c2o27*rho*(1.0+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq);
      //feq[B] =   c2o27*rho*(1.0+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq);
      //feq[NE] =  c1o54*rho*(1.0+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq);
      //feq[SW] =  c1o54*rho*(1.0+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
      //feq[SE] =  c1o54*rho*(1.0+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq);
      //feq[NW] =  c1o54*rho*(1.0+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
      //feq[TE] =  c1o54*rho*(1.0+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq);
      //feq[BW] =  c1o54*rho*(1.0+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
      //feq[BE] =  c1o54*rho*(1.0+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq);
      //feq[TW] =  c1o54*rho*(1.0+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
      //feq[TN] =  c1o54*rho*(1.0+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq);
      //feq[BS] =  c1o54*rho*(1.0+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
      //feq[BN] =  c1o54*rho*(1.0+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq);
      //feq[TS] =  c1o54*rho*(1.0+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
      //feq[TNE] = c1o216*rho*(1.0+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      //feq[BSW] = c1o216*rho*(1.0+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      //feq[BNE] = c1o216*rho*(1.0+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      //feq[TSW] = c1o216*rho*(1.0+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      //feq[TSE] = c1o216*rho*(1.0+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      //feq[BNW] = c1o216*rho*(1.0+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      //feq[BSE] = c1o216*rho*(1.0+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      //feq[TNW] = c1o216*rho*(1.0+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

      //////////////////////////////////////////////////////////////////////////

      LBMReal cu_sq = 1.5*(vx1*vx1+vx2*vx2+vx3*vx3);
      LBMReal rho = drho+one;

      feq[ZERO] = c8o27*(drho+rho*(-cu_sq));
      feq[E] = c2o27*(drho+rho*(3.0*(vx1)+c9o2*(vx1)*(vx1)-cu_sq));
      feq[W] = c2o27*(drho+rho*(3.0*(-vx1)+c9o2*(-vx1)*(-vx1)-cu_sq));
      feq[N] = c2o27*(drho+rho*(3.0*(vx2)+c9o2*(vx2)*(vx2)-cu_sq));
      feq[S] = c2o27*(drho+rho*(3.0*(-vx2)+c9o2*(-vx2)*(-vx2)-cu_sq));
      feq[T] = c2o27*(drho+rho*(3.0*(vx3)+c9o2*(vx3)*(vx3)-cu_sq));
      feq[B] = c2o27*(drho+rho*(3.0*(-vx3)+c9o2*(-vx3)*(-vx3)-cu_sq));
      feq[NE] = c1o54*(drho+rho*(3.0*(vx1+vx2)+c9o2*(vx1+vx2)*(vx1+vx2)-cu_sq));
      feq[SW] = c1o54*(drho+rho*(3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq));
      feq[SE] = c1o54*(drho+rho*(3.0*(vx1-vx2)+c9o2*(vx1-vx2)*(vx1-vx2)-cu_sq));
      feq[NW] = c1o54*(drho+rho*(3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq));
      feq[TE] = c1o54*(drho+rho*(3.0*(vx1+vx3)+c9o2*(vx1+vx3)*(vx1+vx3)-cu_sq));
      feq[BW] = c1o54*(drho+rho*(3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq));
      feq[BE] = c1o54*(drho+rho*(3.0*(vx1-vx3)+c9o2*(vx1-vx3)*(vx1-vx3)-cu_sq));
      feq[TW] = c1o54*(drho+rho*(3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq));
      feq[TN] = c1o54*(drho+rho*(3.0*(vx2+vx3)+c9o2*(vx2+vx3)*(vx2+vx3)-cu_sq));
      feq[BS] = c1o54*(drho+rho*(3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq));
      feq[BN] = c1o54*(drho+rho*(3.0*(vx2-vx3)+c9o2*(vx2-vx3)*(vx2-vx3)-cu_sq));
      feq[TS] = c1o54*(drho+rho*(3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq));
      feq[TNE] = c1o216*(drho+rho*(3.0*(vx1+vx2+vx3)+c9o2*(vx1+vx2+vx3)*(vx1+vx2+vx3)-cu_sq));
      feq[BSW] = c1o216*(drho+rho*(3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      feq[BNE] = c1o216*(drho+rho*(3.0*(vx1+vx2-vx3)+c9o2*(vx1+vx2-vx3)*(vx1+vx2-vx3)-cu_sq));
      feq[TSW] = c1o216*(drho+rho*(3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      feq[TSE] = c1o216*(drho+rho*(3.0*(vx1-vx2+vx3)+c9o2*(vx1-vx2+vx3)*(vx1-vx2+vx3)-cu_sq));
      feq[BNW] = c1o216*(drho+rho*(3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      feq[BSE] = c1o216*(drho+rho*(3.0*(vx1-vx2-vx3)+c9o2*(vx1-vx2-vx3)*(vx1-vx2-vx3)-cu_sq));
      feq[TNW] = c1o216*(drho+rho*(3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
   }
   //////////////////////////////////////////////////////////////////////////
   static LBMReal getIncompFeqForDirection(const int& direction,const LBMReal& drho, const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)	
   {	
      LBMReal cu_sq=1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);

      switch(direction)    
      {		 
         case ZERO : return REAL_CAST( c8o27*(drho-cu_sq));
         case E : return REAL_CAST( c2o27*(drho+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq));
         case W : return REAL_CAST( c2o27*(drho+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq));
         case N : return REAL_CAST( c2o27*(drho+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq));
         case S : return REAL_CAST( c2o27*(drho+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq));
         case T : return REAL_CAST( c2o27*(drho+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq));
         case B : return REAL_CAST( c2o27*(drho+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq));
         case NE : return REAL_CAST( c1o54*(drho+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq));
         case SW : return REAL_CAST( c1o54*(drho+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq));
         case SE : return REAL_CAST( c1o54*(drho+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq));
         case NW : return REAL_CAST( c1o54*(drho+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq));
         case TE : return REAL_CAST( c1o54*(drho+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq));
         case BW : return REAL_CAST( c1o54*(drho+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq));
         case BE : return REAL_CAST( c1o54*(drho+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq));
         case TW : return REAL_CAST( c1o54*(drho+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq));
         case TN : return REAL_CAST( c1o54*(drho+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq));
         case BS : return REAL_CAST( c1o54*(drho+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq));
         case BN : return REAL_CAST( c1o54*(drho+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq));
         case TS : return REAL_CAST( c1o54*(drho+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq));
         case TNE : return REAL_CAST(c1o216*(drho+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
         case BSW : return REAL_CAST(c1o216*(drho+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
         case BNE : return REAL_CAST(c1o216*(drho+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
         case TSW : return REAL_CAST(c1o216*(drho+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
         case TSE : return REAL_CAST(c1o216*(drho+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
         case BNW : return REAL_CAST(c1o216*(drho+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
         case BSE : return REAL_CAST(c1o216*(drho+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
         case TNW : return REAL_CAST(c1o216*(drho+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
         default: throw UbException(UB_EXARGS,"unknown dir");
      }
   }
   //////////////////////////////////////////////////////////////////////////
   static void calcIncompFeq(LBMReal* const& feq/*[27]*/,const LBMReal& drho,const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)	
   {
      LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

      feq[ZERO] =  c8o27*(drho-cu_sq);
      feq[E] =  c2o27*(drho+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq);
      feq[W] =  c2o27*(drho+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq);
      feq[N] =  c2o27*(drho+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq);
      feq[S] =  c2o27*(drho+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq);
      feq[T] =  c2o27*(drho+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq);
      feq[B] =  c2o27*(drho+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq);
      feq[NE] =  c1o54*(drho+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq);
      feq[SW] =  c1o54*(drho+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
      feq[SE] =  c1o54*(drho+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq);
      feq[NW] =  c1o54*(drho+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
      feq[TE] =  c1o54*(drho+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq);
      feq[BW] =  c1o54*(drho+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
      feq[BE] =  c1o54*(drho+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq);
      feq[TW] =  c1o54*(drho+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
      feq[TN] =  c1o54*(drho+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq);
      feq[BS] =  c1o54*(drho+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
      feq[BN] =  c1o54*(drho+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq);
      feq[TS] =  c1o54*(drho+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
      feq[TNE] = c1o216*(drho+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      feq[BSW] = c1o216*(drho+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      feq[BNE] = c1o216*(drho+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      feq[TSW] = c1o216*(drho+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      feq[TSE] = c1o216*(drho+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      feq[BNW] = c1o216*(drho+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      feq[BSE] = c1o216*(drho+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      feq[TNW] = c1o216*(drho+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);   
   }
   //////////////////////////////////////////////////////////////////////////
   static inline float getBoundaryVelocityForDirection(const int& direction, const float& bcVelocityX1,const float& bcVelocityX2,const float& bcVelocityX3)
   {
      switch(direction) 
      {          
      case E:   return (float)( UbMath::c4o9*(+bcVelocityX1) );
      case W:   return (float)( UbMath::c4o9*(-bcVelocityX1) );
      case N:   return (float)( UbMath::c4o9*(+bcVelocityX2) );
      case S:   return (float)( UbMath::c4o9*(-bcVelocityX2) );
      case T:   return (float)( UbMath::c4o9*(+bcVelocityX3) );
      case B:   return (float)( UbMath::c4o9*(-bcVelocityX3) );
      case NE:  return (float)( UbMath::c1o9*(+bcVelocityX1+bcVelocityX2             ) );
      case SW:  return (float)( UbMath::c1o9*(-bcVelocityX1-bcVelocityX2             ) );
      case SE:  return (float)( UbMath::c1o9*(+bcVelocityX1-bcVelocityX2             ) );
      case NW:  return (float)( UbMath::c1o9*(-bcVelocityX1+bcVelocityX2             ) );
      case TE:  return (float)( UbMath::c1o9*(+bcVelocityX1             +bcVelocityX3) );
      case BW:  return (float)( UbMath::c1o9*(-bcVelocityX1             -bcVelocityX3) );
      case BE:  return (float)( UbMath::c1o9*(+bcVelocityX1             -bcVelocityX3) );
      case TW:  return (float)( UbMath::c1o9*(-bcVelocityX1             +bcVelocityX3) );
      case TN:  return (float)( UbMath::c1o9*(             +bcVelocityX2+bcVelocityX3) );
      case BS:  return (float)( UbMath::c1o9*(             -bcVelocityX2-bcVelocityX3) );
      case BN:  return (float)( UbMath::c1o9*(             +bcVelocityX2-bcVelocityX3) );
      case TS:  return (float)( UbMath::c1o9*(             -bcVelocityX2+bcVelocityX3) );
      case TNE: return (float)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2+bcVelocityX3) );
      case BSW: return (float)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case BNE: return (float)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case TSW: return (float)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case TSE: return (float)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case BNW: return (float)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case BSE: return (float)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case TNW: return (float)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2+bcVelocityX3) );
      default: throw UbException(UB_EXARGS,"unknown direction"); 
      }
   }
   /*=====================================================================*/
   static const int& getInvertDirection(const int& direction)
   {  
   #ifdef _DEBUG
      if(direction<STARTDIR || direction>ENDDIR) 
         throw UbException(UB_EXARGS,"unknown direction");
   #endif
      return INVDIR[direction];
   }
   /*=====================================================================*/
   static void getLBMDirections(std::vector<int>& dirs, bool onlyLBdirs = false)
   {
      std::vector<int> D3Q27Dirs;
      if(onlyLBdirs) /*FSTARTDIR->FENDDIR*/
      {
         dirs.resize(FENDDIR+1);
         for(int dir=FSTARTDIR; dir<=FENDDIR; ++dir)
            dirs[dir] = dir;
      }
      else /*STARTDIR->ENDDIR*/
      {
         dirs.resize(ENDDIR+1);
         for(int dir=STARTDIR; dir<=ENDDIR; ++dir)
            dirs[dir] = dir;
      }
   }
//////////////////////////////////////////////////////////////////////////
   static std::vector<int> getEX(const int& exn)
   {
      std::vector<int> ex;
      ex.resize(ENDDIR+1);
      switch (exn)
      {
      case 1:
         for(int dir=STARTDIR; dir<ENDDIR; ++dir)
            ex[dir] = DX1[dir];
      	break;
      case 2:
         for(int dir=STARTDIR; dir<ENDDIR; ++dir)
            ex[dir] = DX2[dir];
         break;
      case 3:
         for(int dir=STARTDIR; dir<ENDDIR; ++dir)
            ex[dir] = DX3[dir];
         break;
      }
      return ex;
   }
//////////////////////////////////////////////////////////////////////////
   static inline void calcDistanceToNeighbors(std::vector<double>& distNeigh, const double& deltaX1)
   {
      //distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);
      double sqrt3 = UbMath::sqrt3;
      double sqrt2 = UbMath::sqrt2;
      distNeigh[E] = distNeigh[W] = distNeigh[N] = deltaX1;
      distNeigh[S] = distNeigh[T] = distNeigh[B] = deltaX1;
      distNeigh[NE] = distNeigh[NW] = distNeigh[SW] = distNeigh[SE] = sqrt2*deltaX1;
      distNeigh[TE] = distNeigh[TN] = distNeigh[TW] = distNeigh[TS] = sqrt2*deltaX1;
      distNeigh[BE] = distNeigh[BN] = distNeigh[BW] = distNeigh[BS] = sqrt2*deltaX1;
      distNeigh[TNE] = distNeigh[TNW] = distNeigh[TSE] = distNeigh[TSW] = sqrt3*deltaX1;
      distNeigh[BNE] = distNeigh[BNW] = distNeigh[BSE] = distNeigh[BSW] = sqrt3*deltaX1;
   }
//////////////////////////////////////////////////////////////////////////
   static inline void calcDistanceToNeighbors(std::vector<double>& distNeigh, const double& deltaX1,const double& deltaX2,const double& deltaX3)
   {
      //distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);
      double sqrt3 = UbMath::sqrt3;
      double sqrt2 = UbMath::sqrt2;
      distNeigh[E] = distNeigh[W] =  deltaX1;
      distNeigh[N] = distNeigh[S] =  deltaX2;
      distNeigh[T] = distNeigh[B] = deltaX3;
      distNeigh[NE] = distNeigh[NW] = distNeigh[SW] = distNeigh[SE] = sqrt(deltaX1*deltaX1+deltaX2*deltaX2);
      distNeigh[TE] = distNeigh[TN] = distNeigh[TW] = distNeigh[TS] = sqrt(deltaX1*deltaX1+deltaX3*deltaX3);
      distNeigh[BE] = distNeigh[BN] = distNeigh[BW] = distNeigh[BS] = sqrt(deltaX2*deltaX2+deltaX3*deltaX3);
      distNeigh[TNE] = distNeigh[TNW] = distNeigh[TSE] = distNeigh[TSW] = sqrt(deltaX1*deltaX1+deltaX2*deltaX2+deltaX3*deltaX3);
      distNeigh[BNE] = distNeigh[BNW] = distNeigh[BSE] = distNeigh[BSW] = sqrt(deltaX1*deltaX1+deltaX2*deltaX2+deltaX3*deltaX3);
   }
//////////////////////////////////////////////////////////////////////////
   static inline LBMReal calcPress(const LBMReal* const f, LBMReal rho, LBMReal vx1, LBMReal vx2, LBMReal vx3)
   {
      LBMReal op=1.0;
      return ((f[E]+f[W]+f[N]+f[S]+f[T]+f[B]+2.*(f[NE]+f[SW]+f[SE]+f[NW]+f[TE]+f[BW]+f[BE]+f[TW]+f[TN]+f[BS]+f[BN]+f[TS])+
         3.*(f[TNE]+f[TSW]+f[TSE]+f[TNW]+f[BNE]+f[BSW]+f[BSE]+f[BNW])-(vx1*vx1+vx2*vx2+vx3*vx3))*(1-0.5*op)+op*0.5*(rho))*c1o3;

   }
}

#endif



