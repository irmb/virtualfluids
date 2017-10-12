#ifndef D3Q27SYSTEM_H
#define D3Q27SYSTEM_H

#include <cmath>
#include <string>
#include <iostream>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbTuple.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbSystem.h>
#include "Patch3DSystem.h"
#include "LBMSystem.h"


class LBMSystems
{
   //////////////////////////////////////////////////////////////////////////


   //////////////////////////////////////////////////////////////////////////
   //MACROSCOPIC VALUES                  
   /*=====================================================================*/
   static LBMReal getDensity(const LBMReal* const& f/*[27]*/)
   {
      return  f[ZERO] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
      + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
      + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
      + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];
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
      return f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] 
      + f[BE] - f[BW] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
      + f[BSE] - f[BNW];  
   }
   /*=====================================================================*/
   static LBMReal getIncompVelocityX2(const LBMReal* const& f/*[27]*/)
   {
      return f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
      - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
      + f[BNW];  
   }
   /*=====================================================================*/
   static LBMReal getIncompVelocityX3(const LBMReal* const& f/*[27]*/)
   {
      return f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
      + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
      - f[BNW];
   }
   /*=====================================================================*/
   static void calcDensity(const LBMReal* const& f/*[27]*/, LBMReal& rho)
   {
      rho = f[ZERO] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
      + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
      + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
      + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];
   }
   /*=====================================================================*/
   static void calcIncompVelocityX1(const LBMReal* const& f/*[27]*/, LBMReal& vx1)
   {
      vx1 = f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
      + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
      + f[BSE] - f[BNW];  
   }
   /*=====================================================================*/
   static void calcIncompVelocityX2(const LBMReal* const& f/*[27]*/, LBMReal& vx2)
   {
      vx2 = f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
      - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
      + f[BNW];  
   }
   /*=====================================================================*/
   static void calcIncompVelocityX3(const LBMReal* const& f/*[27]*/, LBMReal& vx3)
   {
      vx3 = f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
      + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
      - f[BNW];
   }
   /*=====================================================================*/
   static LBMReal getCompVelocityX1(const LBMReal* const& f/*[27]*/)
   {
      return (f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
      + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
      + f[BSE] - f[BNW])/getDensity(f);  
   }
   /*=====================================================================*/
   static LBMReal getCompVelocityX2(const LBMReal* const& f/*[27]*/)
   {
      return (f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
      - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
      + f[BNW])/getDensity(f);  
   }
   /*=====================================================================*/
   static LBMReal getCompVelocityX3(const LBMReal* const& f/*[27]*/)
   {
      return (f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
      + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
      - f[BNW])/getDensity(f);
   }
   /*=====================================================================*/
   static void calcCompVelocityX1(const LBMReal* const& f/*[27]*/, LBMReal& vx1)
   {
      vx1 = (f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
      + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
      + f[BSE] - f[BNW])/getDensity(f);  
   }
   /*=====================================================================*/
   static void calcCompVelocityX2(const LBMReal* const& f/*[27]*/, LBMReal& vx2)
   {
      vx2 = (f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
      - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
      + f[BNW])/getDensity(f);  
   }
   /*=====================================================================*/
   static void calcCompVelocityX3(const LBMReal* const& f/*[27]*/, LBMReal& vx3)
   {
      vx3 = (f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
      + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
      - f[BNW])/getDensity(f);
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
   static void calcCompMacroscopicValues(const LBMReal* const& f/*[27]*/, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3)
   {
      D3Q27System::calcDensity(f, rho);
      D3Q27System::calcIncompVelocityX1(f, vx1);
      D3Q27System::calcIncompVelocityX2(f, vx2);
      D3Q27System::calcIncompVelocityX3(f, vx3);
      vx1/=rho;
      vx2/=rho;
      vx3/=rho;
   }
   //////////////////////////////////////////////////////////////////////////
   static LBMReal getCompFeqForDirection(const int& direction, const LBMReal& rho,const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)
   {
      LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

      switch(direction)    
      {
      case ZERO : return REAL_CAST( c8o27*rho*(1.0-cu_sq));
      case E : return REAL_CAST(  c2o27*rho*(1.0+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq));
      case W : return REAL_CAST(  c2o27*rho*(1.0+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq));
      case N : return REAL_CAST(  c2o27*rho*(1.0+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq));
      case S : return REAL_CAST(  c2o27*rho*(1.0+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq));
      case T : return REAL_CAST(  c2o27*rho*(1.0+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq));
      case B : return REAL_CAST(  c2o27*rho*(1.0+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq));
      case NE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq));
      case SW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq));
      case SE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq));
      case NW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq));
      case TE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq));
      case BW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq));
      case BE : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq));
      case TW : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq));
      case TN : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq));
      case BS : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq));
      case BN : return REAL_CAST( c1o54*rho*(1.0+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq));
      case TS : return REAL_CAST( c1o54*rho*(1.0+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq));
      case TNE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      case BSW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      case BNE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      case TSW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      case TSE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      case BNW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      case BSE : return REAL_CAST(c1o216*rho*(1.0+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      case TNW : return REAL_CAST(c1o216*rho*(1.0+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      default: throw UbException(UB_EXARGS,"unknown dir");
      }
   }
   //////////////////////////////////////////////////////////////////////////
   static void calcCompFeq(LBMReal* const& feq/*[27]*/,const LBMReal& rho,const LBMReal& vx1,const LBMReal& vx2,const LBMReal& vx3)	
   {
      LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

      feq[ZERO] =  c8o27*rho*(1.0-cu_sq);
      feq[E] =   c2o27*rho*(1.0+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq);
      feq[W] =   c2o27*rho*(1.0+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq);
      feq[N] =   c2o27*rho*(1.0+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq);
      feq[S] =   c2o27*rho*(1.0+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq);
      feq[T] =   c2o27*rho*(1.0+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq);
      feq[B] =   c2o27*rho*(1.0+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq);
      feq[NE] =  c1o54*rho*(1.0+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq);
      feq[SW] =  c1o54*rho*(1.0+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
      feq[SE] =  c1o54*rho*(1.0+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq);
      feq[NW] =  c1o54*rho*(1.0+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
      feq[TE] =  c1o54*rho*(1.0+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq);
      feq[BW] =  c1o54*rho*(1.0+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
      feq[BE] =  c1o54*rho*(1.0+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq);
      feq[TW] =  c1o54*rho*(1.0+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
      feq[TN] =  c1o54*rho*(1.0+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq);
      feq[BS] =  c1o54*rho*(1.0+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
      feq[BN] =  c1o54*rho*(1.0+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq);
      feq[TS] =  c1o54*rho*(1.0+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
      feq[TNE] = c1o216*rho*(1.0+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      feq[BSW] = c1o216*rho*(1.0+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      feq[BNE] = c1o216*rho*(1.0+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      feq[TSW] = c1o216*rho*(1.0+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      feq[TSE] = c1o216*rho*(1.0+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      feq[BNW] = c1o216*rho*(1.0+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      feq[BSE] = c1o216*rho*(1.0+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      feq[TNW] = c1o216*rho*(1.0+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
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
      LBMReal cu_sq=1.5f*(vx1*vx1+vx2*vx2+vx3*vx3);

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
   static std::vector<int> getFDirs()
   {
      std::vector<int> D3Q27_f_Dirs;
      D3Q27_f_Dirs.resize(ENDF+1);
      for(int dir=STARTF; dir<=ENDF; ++dir)
         D3Q27_f_Dirs[dir] = dir;

      return D3Q27_f_Dirs;
   }
   /*=====================================================================*/
   static std::vector<int> getDirs(const bool& onlyLBdirs)
   {
      std::vector<int> D3Q27Dirs;
      if(onlyLBdirs) /*FSTARTDIR->FENDDIR*/
      {
         D3Q27Dirs.resize(FENDDIR+1);
         for(int dir=FSTARTDIR; dir<=FENDDIR; ++dir)
            D3Q27Dirs[dir] = dir;
      }
      else /*STARTDIR->ENDDIR*/
      {
         D3Q27Dirs.resize(ENDDIR+1);
         for(int dir=STARTDIR; dir<=ENDDIR; ++dir)
            D3Q27Dirs[dir] = dir;
      }
      return D3Q27Dirs;
   }
   //////////////////////////////////////////////////////////////////////////
   static LBMReal calcCollisionFactor(LBMReal viscosity, LBMReal deltaX)
   {
      return 1.0/(3.0*viscosity/deltaX+0.5);
   }




};

#endif

