//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef BoundaryConditions_H
#define BoundaryConditions_H

#include <vector>
#include <string>

#include <basics/utilities/UbException.h>                  
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbTuple.h>
#include "D3Q27System.h"

#include <boost/serialization/serialization.hpp>

class BoundaryConditions;
typedef boost::shared_ptr<BoundaryConditions> BoundaryConditionsPtr;

class BoundaryConditions 
{
//public:
//   enum BcAlgorithm{VelocityBC, SlipBC, NoSlipBC, ThinWallNoSlipBC, HighViscosityNoSlipBC, EqDensityBC, NonEqDensityBC, NonReflectingVelocityBC, NonReflectingDensityBC};
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
      //wenn folgendes nicht geht, dann hat man weiter unten bei der bit-geschichte ein ernstes problem!!!
      UB_STATIC_ASSERT( sizeof(long long) >= 8);
      //UB_STATIC_ASSERT( sizeof(double) >= 16);
      //UB_STATIC_ASSERT( sizeof(long long) == 32);
      UB_STATIC_ASSERT( (sizeof(long long)*8) >= (D3Q27System::FENDDIR+1)*BoundaryConditions::optionDigits );

      for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++) 
         q[fdir] = -999.; 
   }
   virtual ~BoundaryConditions() {}

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
      
      //alle digits an den betreffenden postionen auf "0"
      flag &= ~( maxOptionVal<<(direction*optionDigits) );
      //alle digitsan den betreffenden postionen entsprechend der marke setzen
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
   void       setNormalVector(const float& nx1,const float& nx2,const float& nx3)   { this->nx1 = nx1; this->nx2 = nx2;  this->nx3 = nx3;}
   UbTupleFloat3 getNormalVector()                                                  { return makeUbTuple(nx1,nx2,nx3); }

   /*============== Velocity Boundary ========================*/
   void       setVelocityBoundaryFlag(const int& direction, const short& secOpt=0)  { this->setFlagBits(velocityBoundaryFlags,direction,secOpt);                                     }  
   void       unsetVelocityBoundaryFlag(const int& direction)                       { this->velocityBoundaryFlags &= ~( maxOptionVal<<(direction*optionDigits) );                    }
   void       unsetVelocityBoundary()                 		                        { this->velocityBoundaryFlags = 0;                                                               }
   long long  getVelocityBoundary()	               		                           { return this->velocityBoundaryFlags;                                                            }
   bool       hasVelocityBoundary()   					 		                           { return this->velocityBoundaryFlags!=0;                                                         }
   bool       hasVelocityBoundaryFlag(const int& direction)                         { return ( ( ( velocityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) != 0);         }
   short      getVelocitySecondaryOption(const int& direction)                      { return (short)( (  ( velocityBoundaryFlags>>(optionDigits*direction) ) & maxOptionVal ) - 1 ); }

   void  setBoundaryVelocityX1(const float& vx1) { this->bcVelocityX1 = vx1;  } 
   void  setBoundaryVelocityX2(const float& vx2) { this->bcVelocityX2 = vx2;  } 
   void  setBoundaryVelocityX3(const float& vx3) { this->bcVelocityX3 = vx3;  } 
   float getBoundaryVelocityX1()                 { return this->bcVelocityX1; }
   float getBoundaryVelocityX2()                 { return this->bcVelocityX2; }
   float getBoundaryVelocityX3()                 { return this->bcVelocityX3; }
   float getBoundaryVelocity(const int& direction) 
   {                   
      switch(direction)
      {
      case D3Q27System::E : return (float)( UbMath::c4o9*(+bcVelocityX1) );      //(2/cs^2)(=6)*rho_0(=1 bei inkompr)*wi*u*ei mit cs=1/sqrt(3)
      case D3Q27System::W : return (float)( UbMath::c4o9*(-bcVelocityX1) );      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
      case D3Q27System::N : return (float)( UbMath::c4o9*(+bcVelocityX2) );   
      case D3Q27System::S : return (float)( UbMath::c4o9*(-bcVelocityX2) );
      case D3Q27System::T : return (float)( UbMath::c4o9*(+bcVelocityX3) );
      case D3Q27System::B : return (float)( UbMath::c4o9*(-bcVelocityX3) );
      case D3Q27System::NE: return (float)( UbMath::c1o9*(+bcVelocityX1+bcVelocityX2             ) );
      case D3Q27System::SW: return (float)( UbMath::c1o9*(-bcVelocityX1-bcVelocityX2             ) );
      case D3Q27System::SE: return (float)( UbMath::c1o9*(+bcVelocityX1-bcVelocityX2             ) );
      case D3Q27System::NW: return (float)( UbMath::c1o9*(-bcVelocityX1+bcVelocityX2             ) );
      case D3Q27System::TE: return (float)( UbMath::c1o9*(+bcVelocityX1             +bcVelocityX3) );
      case D3Q27System::BW: return (float)( UbMath::c1o9*(-bcVelocityX1             -bcVelocityX3) );
      case D3Q27System::BE: return (float)( UbMath::c1o9*(+bcVelocityX1             -bcVelocityX3) );
      case D3Q27System::TW: return (float)( UbMath::c1o9*(-bcVelocityX1             +bcVelocityX3) );
      case D3Q27System::TN: return (float)( UbMath::c1o9*(             +bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BS: return (float)( UbMath::c1o9*(             -bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BN: return (float)( UbMath::c1o9*(             +bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TS: return (float)( UbMath::c1o9*(             -bcVelocityX2+bcVelocityX3) );
      case D3Q27System::TNE: return (float)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BSW: return (float)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BNE: return (float)( UbMath::c1o36*(+bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TSW: return (float)( UbMath::c1o36*(-bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case D3Q27System::TSE: return (float)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2+bcVelocityX3) );
      case D3Q27System::BNW: return (float)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2-bcVelocityX3) );
      case D3Q27System::BSE: return (float)( UbMath::c1o36*(+bcVelocityX1-bcVelocityX2-bcVelocityX3) );
      case D3Q27System::TNW: return (float)( UbMath::c1o36*(-bcVelocityX1+bcVelocityX2+bcVelocityX3) ); 
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

   void  setBoundaryDensity(float density) { this->bcDensity = density; } 
   float getBoundaryDensity()              { return this->bcDensity;    }

   //Lodi extension
   void  setDensityLodiDensity(const float& bcLodiDensity)       { this->bcLodiDensity    = bcLodiDensity;    } 
   void  setDensityLodiVelocityX1(const float& bcLodiVelocityX1) { this->bcLodiVelocityX1 = bcLodiVelocityX1; } 
   void  setDensityLodiVelocityX2(const float& bcLodiVelocityX2) { this->bcLodiVelocityX2 = bcLodiVelocityX2; } 
   void  setDensityLodiVelocityX3(const float& bcLodiVelocityX3) { this->bcLodiVelocityX3 = bcLodiVelocityX3; } 
   void  setDensityLodiLength(const float& bcLodiLentgh)         { this->bcLodiLentgh     = bcLodiLentgh;     } 
   float getDensityLodiDensity() const                           { return this->bcLodiDensity;    } 
   float getDensityLodiVelocityX1() const                        { return this->bcLodiVelocityX1; }
   float getDensityLodiVelocityX2() const                        { return this->bcLodiVelocityX2; }
   float getDensityLodiVelocityX3() const                        { return this->bcLodiVelocityX3; }
   float getDensityLodiLength() const                            { return this->bcLodiLentgh;     }

   float& densityLodiDensity()                                   { return this->bcLodiDensity;    } 
   float& densityLodiVelocityX1()                                { return this->bcLodiVelocityX1; }
   float& densityLodiVelocityX2()                                { return this->bcLodiVelocityX2; }
   float& densityLodiVelocityX3()                                { return this->bcLodiVelocityX3; }
   float& densityLodiLentgh()                                    { return this->bcLodiLentgh;     }

   const float& densityLodiDensity()  const                      { return this->bcLodiDensity;    } 
   const float& densityLodiVelocityX1() const                    { return this->bcLodiVelocityX1; }
   const float& densityLodiVelocityX2() const                    { return this->bcLodiVelocityX2; }
   const float& densityLodiVelocityX3() const                    { return this->bcLodiVelocityX3; }
   const float& densityLodiLentgh()  const                       { return this->bcLodiLentgh;     }


   /*======================= Qs =============================*/
   void  setQ(const float& val, const int& direction) { q[direction] = val; }
   float getQ(const int& direction)                   { return q[direction]; }
   
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
   static const int       optionDigits = 2;  //--> 3 bits für secondary Option --> maxOptionVal = 7, da man mit drei Digits max die 7 darstellen kann
   static const long long maxOptionVal;// = ( 1<<optionDigits ) - 1; //2^3-1 -> 7

protected:
   float q[D3Q27System::FENDDIR+1];
   //float q[D3Q27System::STARTF+1];

   long long noslipBoundaryFlags;		
   long long slipBoundaryFlags;		
   long long velocityBoundaryFlags;		
   long long densityBoundaryFlags;		
   long long wallModelBoundaryFlags;

   float  bcVelocityX1;
   float  bcVelocityX2;
   float  bcVelocityX3;
   float  bcDensity;

   float  bcLodiDensity;
   float  bcLodiVelocityX1;
   float  bcLodiVelocityX2;
   float  bcLodiVelocityX3;
   float  bcLodiLentgh;

   float  nx1,nx2,nx3;

   char algorithmType;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & q; 

      ar & noslipBoundaryFlags;		
      ar & slipBoundaryFlags;		
      ar & velocityBoundaryFlags;		
      ar & densityBoundaryFlags;		

      ar & bcVelocityX1;
      ar & bcVelocityX2;
      ar & bcVelocityX3;
      ar & bcDensity;

      ar & bcLodiDensity;
      ar & bcLodiVelocityX1;
      ar & bcLodiVelocityX2;
      ar & bcLodiVelocityX3;
      ar & bcLodiLentgh;

      ar & wallModelBoundaryFlags;

      ar & nx1;
      ar & nx2;
      ar & nx3;

      ar & algorithmType;
   }

};

#endif //D3Q27BOUNDARYCONDITION_H
