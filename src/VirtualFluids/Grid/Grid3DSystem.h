//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef Grid3DSystem_H
#define Grid3DSystem_H

#include <cmath>
#include <iostream>
#include <string>

#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbException.h>


namespace Grid3DSystem
{
   static const int STARTDIR = 0;

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

   static const int ENDDIR = 25; 

   static const int INV_E   = W;  
   static const int INV_W   = E;  
   static const int INV_N   = S;  
   static const int INV_S   = N;  
   static const int INV_T   = B;  
   static const int INV_B   = T;  
   static const int INV_NE  = SW; 
   static const int INV_NW  = SE; 
   static const int INV_SE  = NW; 
   static const int INV_SW  = NE; 
   static const int INV_TE  = BW; 
   static const int INV_TW  = BE; 
   static const int INV_BE  = TW; 
   static const int INV_BW  = TE; 
   static const int INV_TN  = BS; 
   static const int INV_TS  = BN; 
   static const int INV_BN  = TS; 
   static const int INV_BS  = TN; 
   static const int INV_TNE = BSW;
   static const int INV_TNW = BSE;
   static const int INV_TSE = BNW;
   static const int INV_TSW = BNE;
   static const int INV_BNE = TSW;
   static const int INV_BNW = TSE;
   static const int INV_BSE = TNW;
   static const int INV_BSW = TNE;

   extern const int INVDIR[ENDDIR+1];

   static const int MAXLEVEL  = 25;

   extern const int EX1[ENDDIR+1];
   extern const int EX2[ENDDIR+1];
   extern const int EX3[ENDDIR+1];

   inline std::string getDirectionString(int direction)
   {
      switch(direction)
      {
      case E   : return "E"; 
      case W   : return "W"; 
      case N   : return "N"; 
      case S   : return "S"; 
      case T   : return "T";
      case B   : return "B"; 
      case NE  : return "NE";
      case NW  : return "NW";
      case SE  : return "SE";
      case SW  : return "SW";
      case TE  : return "TE";
      case TW  : return "TW";
      case BE  : return "BE";
      case BW  : return "BW";
      case TN  : return "TN";
      case TS  : return "TS";
      case BN  : return "BN";
      case BS  : return "BS";
      case TNE : return "TNE";
      case TNW : return "TNW";
      case TSE : return "TSE";
      case TSW : return "TSW";
      case BNE : return "BNE";
      case BNW : return "BNW";
      case BSE : return "BSE";
      case BSW : return "BSW";
      default  : return "Cell3DSystem::getDrectionString(...) - unknown dir";
      }
   }
   static const int&       getInvertDirection(const int& direction);

//////////////////////////////////////////////////////////////////////////
   static inline void setNeighborCoordinatesForDirection(int &x1, int &x2,int &x3, const int& direction)
   {
      switch(direction)
      {
      case Grid3DSystem::E  :  x1++;             break;
      case Grid3DSystem::N  :  x2++;             break;
      case Grid3DSystem::T  :  x3++;             break;
      case Grid3DSystem::W  :  x1--;             break;
      case Grid3DSystem::S  :  x2--;             break;
      case Grid3DSystem::B  :  x3--;             break;
      case Grid3DSystem::NE :  x1++; x2++;       break;
      case Grid3DSystem::NW :  x1--; x2++;       break;
      case Grid3DSystem::SW :  x1--; x2--;       break;
      case Grid3DSystem::SE :  x1++; x2--;       break;
      case Grid3DSystem::TE :  x1++; x3++;       break;
      case Grid3DSystem::BW :  x1--; x3--;       break;
      case Grid3DSystem::BE :  x1++; x3--;       break;
      case Grid3DSystem::TW :  x1--; x3++;       break;
      case Grid3DSystem::TN :  x2++; x3++;       break;
      case Grid3DSystem::BS :  x2--; x3--;       break;
      case Grid3DSystem::BN :  x2++; x3--;       break;
      case Grid3DSystem::TS :  x2--; x3++;       break;
      case Grid3DSystem::TNE:  x1++; x2++; x3++; break;
      case Grid3DSystem::TNW:  x1--; x2++; x3++; break;
      case Grid3DSystem::TSE:  x1++; x2--; x3++; break;
      case Grid3DSystem::TSW:  x1--; x2--; x3++; break;
      case Grid3DSystem::BNE:  x1++; x2++; x3--; break;
      case Grid3DSystem::BNW:  x1--; x2++; x3--; break;
      case Grid3DSystem::BSE:  x1++; x2--; x3--; break;
      case Grid3DSystem::BSW:  x1--; x2--; x3--; break;
      default: throw UbException(UB_EXARGS,"no direction ...");
      }
   }
}

#endif 
