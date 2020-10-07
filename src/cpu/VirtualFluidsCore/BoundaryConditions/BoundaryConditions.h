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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file BoundaryConditions.h
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================
#ifndef BoundaryConditions_H
#define BoundaryConditions_H

#include <vector>
#include <string>

#include "Vector3D.h"
#include "UbException.h"                  
#include "UbSystem.h"
#include "UbTuple.h"
#include "D3Q27System.h"
#include <PointerDefinitions.h>

//! Difenition of baundary conditions in grid generation
class BoundaryConditions 
{
public:
   BoundaryConditions() 
      : noslipBoundaryFlags(0)		
      , slipBoundaryFlags(0)		
      , velocityBoundaryFlags(0)		
      , densityBoundaryFlags(0)		
      , wallModelBoundaryFlags(0)
      , bcVelocityX1(0.0f)
      , bcVelocityX2(0.0f)
      , bcVelocityX3(0.0f)
      , bcDensity(0.0f)
      , bcLodiDensity(0.0f)
      , bcLodiVelocityX1(0.0f)
      , bcLodiVelocityX2(0.0f)
      , bcLodiVelocityX3(0.0f)
      , bcLodiLentgh(0.0f)
      , nx1(0.0f)
      , nx2(0.0f)
      , nx3(0.0f)
      , algorithmType(-1)
   {
      UB_STATIC_ASSERT( sizeof(long long) >= 8);
      UB_STATIC_ASSERT( (sizeof(long long)*8) >= (D3Q27System::FENDDIR+1)*BoundaryConditions::optionDigits );

      for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++) 
         q[fdir] = -999.; 
   }
   virtual ~BoundaryConditions() = default;

   virtual bool isEmpty() { return (noslipBoundaryFlags&slipBoundaryFlags&velocityBoundaryFlags&densityBoundaryFlags)==0;}
   virtual bool hasBoundaryCondition()
   {
      return (  hasNoSlipBoundary() || hasSlipBoundary() 
             || hasDensityBoundary() || hasVelocityBoundary() || hasWallModelBoundary() );
   }

   virtual bool hasBoundaryConditionFlag(const int& direction)
   {
      assert( direction >= D3Q27System::FSTARTDIR && direction <= D3Q27System::FENDDIR );

      return (   hasNoSlipBoundaryFlag(direction) || hasSlipBoundaryFlag(direction) 
              || hasDensityBoundaryFlag(direction) || hasVelocityBoundaryFlag(direction)  || hasWallModelBoundaryFlag(direction));
   }
protected:
   void setFlagBits(long long& flag, const int& direction, const short& secOpt)
   {
      if( (secOpt+1)>maxOptionVal ) 
         throw UbException(UB_EXARGS,"error: option > "+UbSystem::toString(maxOptionVal-1));
      
      //all digits at the respective positions to "0"
      flag &= ~( maxOptionVal<<(direction*optionDigits) );
      //set all digits according to the flag at the respective positions
      flag |= ((long long)(secOpt+1)<<(direction*optionDigits));
   }
public:
   /*===================== NoSlip Boundary ==================================================*/	
   void       setNoSlipBoundaryFlag(const int& direction, const short& secOpt=0)    { this->setFlagBits(noslipBoundaryFlags,direction,secOpt);                                     }  
   void       unsetNoSlipBoundaryFlag(const int& direction)                         { this->noslipBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetNoSlipBoundary()                                                 { this->noslipBoundaryFlags = 0;                                                               }
   long long  getNoSlipBoundary()	                                                { return this->noslipBoundaryFlags;                                                            }
   bool       hasNoSlipBoundary()					                                    { return (noslipBoundaryFlags!=0);                                                             }
   bool       hasNoSlipBoundaryFlag(const int& direction)                           { return ( ( ( noslipBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getNoSlipSecondaryOption(const int& direction)                        { return (short)( (  ( noslipBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }
   /*===================== WallModel Boundary ==================================================*/	
   void       setWallModelBoundaryFlag(const int& direction, const short& secOpt=0) { this->setFlagBits(wallModelBoundaryFlags,direction,secOpt);                                     }  
   void       unsetWallModelBoundaryFlag(const int& direction)                      { this->wallModelBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetWallModelBoundary()                                              { this->wallModelBoundaryFlags = 0;                                                               }
   long long  getWallModelBoundary()	                                             { return this->wallModelBoundaryFlags;                                                            }
   bool       hasWallModelBoundary()					                                 { return (wallModelBoundaryFlags!=0);                                                             }
   bool       hasWallModelBoundaryFlag(const int& direction)                        { return ( ( ( wallModelBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getWallModelSecondaryOption(const int& direction)                     { return (short)( (  ( wallModelBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }
   /*===================== Slip-Solid Boundary ==================================================*/	
   void       setSlipBoundaryFlag(const int& direction, const short& secOpt=0)      { this->setFlagBits(slipBoundaryFlags,direction,secOpt);                                     }  
   void       unsetSlipBoundaryFlag(const int& direction)                           { this->slipBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetSlipBoundary()                                                   { this->slipBoundaryFlags = 0;                                                               }
   long long  getSlipBoundary()	                                                   { return this->slipBoundaryFlags;                                                            }
   bool       hasSlipBoundary()					                                       { return (slipBoundaryFlags!=0);                                                             }
   bool       hasSlipBoundaryFlag(const int& direction)	                           { return ( ( ( slipBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getSlipSecondaryOption(const int& direction)                          { return (short)( (  ( slipBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }
   void       setNormalVector(const LBMReal& nx1,const LBMReal& nx2,const LBMReal& nx3)   { this->nx1 = nx1; this->nx2 = nx2;  this->nx3 = nx3;}
   UbTupleDouble3 getNormalVector()                                                  { return makeUbTuple(nx1,nx2,nx3); }

   /*============== Velocity Boundary ========================*/
   void       setVelocityBoundaryFlag(const int& direction, const short& secOpt=0)  { this->setFlagBits(velocityBoundaryFlags,direction,secOpt);                                     }  
   void       unsetVelocityBoundaryFlag(const int& direction)                       { this->velocityBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetVelocityBoundary()                 		                        { this->velocityBoundaryFlags = 0;                                                               }
   long long  getVelocityBoundary()	               		                           { return this->velocityBoundaryFlags;                                                            }
   bool       hasVelocityBoundary()   					 		                           { return this->velocityBoundaryFlags!=0;                                                         }
   bool       hasVelocityBoundaryFlag(const int& direction)                         { return ( ( ( velocityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getVelocitySecondaryOption(const int& direction)                      { return (short)( (  ( velocityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }

   void setBoundaryVelocity(const Vector3D& vx) 
    {
       setBoundaryVelocityX1((LBMReal)vx[0]); 
       setBoundaryVelocityX2((LBMReal)vx[1]);
       setBoundaryVelocityX3((LBMReal)vx[2]);
   }
   void  setBoundaryVelocityX1(const LBMReal& vx1) { this->bcVelocityX1 = vx1;  } 
   void  setBoundaryVelocityX2(const LBMReal& vx2) { this->bcVelocityX2 = vx2;  } 
   void  setBoundaryVelocityX3(const LBMReal& vx3) { this->bcVelocityX3 = vx3;  } 
   LBMReal getBoundaryVelocityX1()                 { return this->bcVelocityX1; }
   LBMReal getBoundaryVelocityX2()                 { return this->bcVelocityX2; }
   LBMReal getBoundaryVelocityX3()                 { return this->bcVelocityX3; }
   LBMReal getBoundaryVelocity(const int& direction) 
   {                   
      switch(direction)
      {
      case D3Q27System::E : return (LBMReal)( UbMath::c4o9*(+bcVelocityX1) );      //(2/cs^2)(=6)*rho_0(=1 for incompressible)*wi*u*ei with cs=1/sqrt(3)
      case D3Q27System::W : return (LBMReal)( UbMath::c4o9*(-bcVelocityX1) );         
      case D3Q27System::N : return (LBMReal)( UbMath::c4o9*(+bcVelocityX2) );   
      case D3Q27System::S : return (LBMReal)( UbMath::c4o9*(-bcVelocityX2) );
      case D3Q27System::T : return (LBMReal)( UbMath::c4o9*(+bcVelocityX3) );
      case D3Q27System::B : return (LBMReal)( UbMath::c4o9*(-bcVelocityX3) );
      case D3Q27System::NE: return (LBMReal)( UbMath::c1o9*(+bcVelocityX1+bcVelocityX2             ) );
      case D3Q27System::SW: return (LBMReal)( UbMath::c1o9*(-bcVelocityX1-bcVelocityX2             ) );
      case D3Q27System::SE: return (LBMReal)( UbMath::c1o9*(+bcVelocityX1-bcVelocityX2             ) );
      case D3Q27System::NW: return (LBMReal)( UbMath::c1o9*(-bcVelocityX1+bcVelocityX2             ) );
      case D3Q27System::TE: return (LBMReal)( UbMath::c1o9*(+bcVelocityX1             +bcVelocityX3) );
      case D3Q27System::BW: return (LBMReal)( UbMath::c1o9*(-bcVelocityX1             -bcVelocityX3) );
      case D3Q27System::BE: return (LBMReal)( UbMath::c1o9*(+bcVelocityX1             -bcVelocityX3) );
      case D3Q27System::TW: return (LBMReal)( UbMath::c1o9*(-bcVelocityX1             +bcVelocityX3) );
      case D3Q27System::TN: return (LBMReal)( UbMath::c1o9*(             +bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BS: return (LBMReal)( UbMath::c1o9*(             -bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BN: return (LBMReal)( UbMath::c1o9*(             +bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TS: return (LBMReal)( UbMath::c1o9*(             -bcVelocityX2+bcVelocityX3) );
      case D3Q27System::TNE: return (LBMReal)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BSW: return (LBMReal)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BNE: return (LBMReal)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TSW: return (LBMReal)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case D3Q27System::TSE: return (LBMReal)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BNW: return (LBMReal)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BSE: return (LBMReal)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TNW: return (LBMReal)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2+bcVelocityX3) ); 
      default: throw UbException(UB_EXARGS,"unknown error");
      }
   }

   /*============== Density Boundary ========================*/
   void       setDensityBoundaryFlag(const int& direction, const short& secOpt=0) { this->setFlagBits(densityBoundaryFlags,direction,secOpt);                                     }  
   void       unsetDensityBoundaryFlag(const int& direction)                      { this->densityBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetDensityBoundary()                                              { this->densityBoundaryFlags = 0;                                                               }
   long long  getDensityBoundary()	                                              { return this->densityBoundaryFlags;                                                            }
   bool       hasDensityBoundary()					                                  { return (this->densityBoundaryFlags!=0);                                                       }
   bool       hasDensityBoundaryFlag(const int& direction)	                      { return ( ( ( densityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getDensitySecondaryOption(const int& direction)                     { return (short)( (  ( densityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }

   void  setBoundaryDensity(LBMReal density) { this->bcDensity = density; } 
   LBMReal getBoundaryDensity()              { return this->bcDensity;    }

   //Lodi extension
   void  setDensityLodiDensity(const LBMReal& bcLodiDensity)       { this->bcLodiDensity    = bcLodiDensity;    } 
   void  setDensityLodiVelocityX1(const LBMReal& bcLodiVelocityX1) { this->bcLodiVelocityX1 = bcLodiVelocityX1; } 
   void  setDensityLodiVelocityX2(const LBMReal& bcLodiVelocityX2) { this->bcLodiVelocityX2 = bcLodiVelocityX2; } 
   void  setDensityLodiVelocityX3(const LBMReal& bcLodiVelocityX3) { this->bcLodiVelocityX3 = bcLodiVelocityX3; } 
   void  setDensityLodiLength(const LBMReal& bcLodiLentgh)         { this->bcLodiLentgh     = bcLodiLentgh;     } 
   LBMReal getDensityLodiDensity() const                           { return this->bcLodiDensity;    } 
   LBMReal getDensityLodiVelocityX1() const                        { return this->bcLodiVelocityX1; }
   LBMReal getDensityLodiVelocityX2() const                        { return this->bcLodiVelocityX2; }
   LBMReal getDensityLodiVelocityX3() const                        { return this->bcLodiVelocityX3; }
   LBMReal getDensityLodiLength() const                            { return this->bcLodiLentgh;     }

   LBMReal& densityLodiDensity()                                   { return this->bcLodiDensity;    } 
   LBMReal& densityLodiVelocityX1()                                { return this->bcLodiVelocityX1; }
   LBMReal& densityLodiVelocityX2()                                { return this->bcLodiVelocityX2; }
   LBMReal& densityLodiVelocityX3()                                { return this->bcLodiVelocityX3; }
   LBMReal& densityLodiLentgh()                                    { return this->bcLodiLentgh;     }

   const LBMReal& densityLodiDensity()  const                      { return this->bcLodiDensity;    } 
   const LBMReal& densityLodiVelocityX1() const                    { return this->bcLodiVelocityX1; }
   const LBMReal& densityLodiVelocityX2() const                    { return this->bcLodiVelocityX2; }
   const LBMReal& densityLodiVelocityX3() const                    { return this->bcLodiVelocityX3; }
   const LBMReal& densityLodiLentgh()  const                       { return this->bcLodiLentgh;     }


   /*======================= Qs =============================*/
   void  setQ(const LBMReal& val, const int& direction) { q[direction] = val; }
   LBMReal getQ(const int& direction)                   { return q[direction]; }
   
   virtual std::vector< std::string > getBCNames()
   {
      std::vector< std::string > tmp;
      tmp.push_back( "NoSlipBC"   );
      tmp.push_back( "SlipBC"     );
      tmp.push_back( "VelocityBC" );
      tmp.push_back( "DensityBC"  );
      return tmp;
   }
   virtual std::vector< long long > getBCFlags()
   {
      std::vector< long long > tmp;
      tmp.push_back( noslipBoundaryFlags   );
      tmp.push_back( slipBoundaryFlags     );
      tmp.push_back( velocityBoundaryFlags );
      tmp.push_back( densityBoundaryFlags  );
      return tmp;
   }

   static bool hasFlagForDirection(const long long& flag, const int& direction)
   {
      return ( ( ( flag>>(optionDigits*direction) ) & maxOptionVal ) != 0);
   }

   void setBcAlgorithmType(char alg) { algorithmType = alg; }
   char getBcAlgorithmType() { return algorithmType; }

public:
   static const int       optionDigits = 2;  //--> 2 bits for secondary Option --> maxOptionVal = 7
   static const long long maxOptionVal;// = ( 1<<optionDigits ) - 1; //2^3-1 -> 7

protected:
   LBMReal q[D3Q27System::FENDDIR+1];

   long long noslipBoundaryFlags;		
   long long slipBoundaryFlags;		
   long long velocityBoundaryFlags;		
   long long densityBoundaryFlags;		
   long long wallModelBoundaryFlags;

   LBMReal  bcVelocityX1;
   LBMReal  bcVelocityX2;
   LBMReal  bcVelocityX3;
   LBMReal  bcDensity;

   LBMReal  bcLodiDensity;
   LBMReal  bcLodiVelocityX1;
   LBMReal  bcLodiVelocityX2;
   LBMReal  bcLodiVelocityX3;
   LBMReal  bcLodiLentgh;

   LBMReal  nx1,nx2,nx3;

   char algorithmType;
};

#endif
