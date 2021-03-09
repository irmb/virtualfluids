//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"

#include "lbm/CalcMac.h"


__device__ Distributions27 getDistributions(real* DD, unsigned int size_Mat, bool evenOrOdd)
{
    Distributions27 D;
    if (evenOrOdd)
    {
        D.f[dirE   ] = &DD[dirE   *size_Mat];
        D.f[dirW   ] = &DD[dirW   *size_Mat];
        D.f[dirN   ] = &DD[dirN   *size_Mat];
        D.f[dirS   ] = &DD[dirS   *size_Mat];
        D.f[dirT   ] = &DD[dirT   *size_Mat];
        D.f[dirB   ] = &DD[dirB   *size_Mat];
        D.f[dirNE  ] = &DD[dirNE  *size_Mat];
        D.f[dirSW  ] = &DD[dirSW  *size_Mat];
        D.f[dirSE  ] = &DD[dirSE  *size_Mat];
        D.f[dirNW  ] = &DD[dirNW  *size_Mat];
        D.f[dirTE  ] = &DD[dirTE  *size_Mat];
        D.f[dirBW  ] = &DD[dirBW  *size_Mat];
        D.f[dirBE  ] = &DD[dirBE  *size_Mat];
        D.f[dirTW  ] = &DD[dirTW  *size_Mat];
        D.f[dirTN  ] = &DD[dirTN  *size_Mat];
        D.f[dirBS  ] = &DD[dirBS  *size_Mat];
        D.f[dirBN  ] = &DD[dirBN  *size_Mat];
        D.f[dirTS  ] = &DD[dirTS  *size_Mat];
        D.f[dirZERO] = &DD[dirZERO*size_Mat];
        D.f[dirTNE ] = &DD[dirTNE *size_Mat];
        D.f[dirTSW ] = &DD[dirTSW *size_Mat];
        D.f[dirTSE ] = &DD[dirTSE *size_Mat];
        D.f[dirTNW ] = &DD[dirTNW *size_Mat];
        D.f[dirBNE ] = &DD[dirBNE *size_Mat];
        D.f[dirBSW ] = &DD[dirBSW *size_Mat];
        D.f[dirBSE ] = &DD[dirBSE *size_Mat];
        D.f[dirBNW ] = &DD[dirBNW *size_Mat];
    }
    else
    {
        D.f[dirW   ] = &DD[dirE   *size_Mat];
        D.f[dirE   ] = &DD[dirW   *size_Mat];
        D.f[dirS   ] = &DD[dirN   *size_Mat];
        D.f[dirN   ] = &DD[dirS   *size_Mat];
        D.f[dirB   ] = &DD[dirT   *size_Mat];
        D.f[dirT   ] = &DD[dirB   *size_Mat];
        D.f[dirSW  ] = &DD[dirNE  *size_Mat];
        D.f[dirNE  ] = &DD[dirSW  *size_Mat];
        D.f[dirNW  ] = &DD[dirSE  *size_Mat];
        D.f[dirSE  ] = &DD[dirNW  *size_Mat];
        D.f[dirBW  ] = &DD[dirTE  *size_Mat];
        D.f[dirTE  ] = &DD[dirBW  *size_Mat];
        D.f[dirTW  ] = &DD[dirBE  *size_Mat];
        D.f[dirBE  ] = &DD[dirTW  *size_Mat];
        D.f[dirBS  ] = &DD[dirTN  *size_Mat];
        D.f[dirTN  ] = &DD[dirBS  *size_Mat];
        D.f[dirTS  ] = &DD[dirBN  *size_Mat];
        D.f[dirBN  ] = &DD[dirTS  *size_Mat];
        D.f[dirZERO] = &DD[dirZERO*size_Mat];
        D.f[dirTNE ] = &DD[dirBSW *size_Mat];
        D.f[dirTSW ] = &DD[dirBNE *size_Mat];
        D.f[dirTSE ] = &DD[dirBNW *size_Mat];
        D.f[dirTNW ] = &DD[dirBSE *size_Mat];
        D.f[dirBNE ] = &DD[dirTSW *size_Mat];
        D.f[dirBSW ] = &DD[dirTNE *size_Mat];
        D.f[dirBSE ] = &DD[dirTNW *size_Mat];
        D.f[dirBNW ] = &DD[dirTSE *size_Mat];
    }
    return D;
}

struct Distribution27
{
   real f[27];
};

__device__ Distribution27 getDistribution(real* DD, unsigned int size_Mat, bool evenOrOdd, unsigned int k,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ)
{
    unsigned int kzero= k;
    unsigned int ke   = k;
    unsigned int kw   = neighborX[k];
    unsigned int kn   = k;
    unsigned int ks   = neighborY[k];
    unsigned int kt   = k;
    unsigned int kb   = neighborZ[k];
    unsigned int ksw  = neighborY[kw];
    unsigned int kne  = k;
    unsigned int kse  = ks;
    unsigned int knw  = kw;
    unsigned int kbw  = neighborZ[kw];
    unsigned int kte  = k;
    unsigned int kbe  = kb;
    unsigned int ktw  = kw;
    unsigned int kbs  = neighborZ[ks];
    unsigned int ktn  = k;
    unsigned int kbn  = kb;
    unsigned int kts  = ks;
    unsigned int ktse = ks;
    unsigned int kbnw = kbw;
    unsigned int ktnw = kw;
    unsigned int kbse = kbs;
    unsigned int ktsw = ksw;
    unsigned int kbne = kb;
    unsigned int ktne = k;
    unsigned int kbsw = neighborZ[ksw];

    Distributions27 D = getDistributions(DD, size_Mat, evenOrOdd);

//       real f_dirE = (D.f[dirE])[ke];
//       real f_dirW = (D.f[dirW])[kw];
//       real f_dirN = (D.f[dirN])[kn];
//       real f_dirS = (D.f[dirS])[ks];
//       real f_dirT = (D.f[dirT])[kt];
//       real f_dirB = (D.f[dirB])[kb];
//
//       real f_dirNE = (D.f[dirNE])[kne];
//       real f_dirSW = (D.f[dirSW])[ksw];
//       real f_dirSE = (D.f[dirSE])[kse];
//       real f_dirNW = (D.f[dirNW])[knw];
//       real f_dirTE = (D.f[dirTE])[kte];
//       real f_dirBW = (D.f[dirBW])[kbw];
//
//       real f_dirBE = (D.f[dirBE])[kbe];
//       real f_dirTW = (D.f[dirTW])[ktw];
//       real f_dirTN = (D.f[dirTN])[ktn];
//       real f_dirBS = (D.f[dirBS])[kbs];
//       real f_dirBN = (D.f[dirBN])[kbn];
//       real f_dirTS = (D.f[dirTS])[kts];
//
//       real f_dirZERO = (D.f[dirZERO])[kzero];
//
//       real f_dirTNE = (D.f[dirTNE ])[ktne];
//       real f_dirTSW = (D.f[dirTSW ])[ktsw];
//       real f_dirTSE = (D.f[dirTSE ])[ktse];
//       real f_dirTNW = (D.f[dirTNW ])[ktnw];
//       real f_dirBNE = (D.f[dirBNE ])[kbne];
//       real f_dirBSW = (D.f[dirBSW ])[kbsw];
//       real f_dirBSE = (D.f[dirBSE ])[kbse];
//       real f_dirBNW = (D.f[dirBNW ])[kbnw];


return {
        (D.f[dirE])[ke],
        (D.f[dirW])[kw],
        (D.f[dirN])[kn],
        (D.f[dirS])[ks],
        (D.f[dirT])[kt],
        (D.f[dirB])[kb],
        (D.f[dirNE])[kne],
        (D.f[dirSW])[ksw],
        (D.f[dirSE])[kse],
        (D.f[dirNW])[knw],
        (D.f[dirTE])[kte],
        (D.f[dirBW])[kbw],
        (D.f[dirBE])[kbe],
        (D.f[dirTW])[ktw],
        (D.f[dirTN])[ktn],
        (D.f[dirBS])[kbs],
        (D.f[dirBN])[kbn],
        (D.f[dirTS])[kts],
        (D.f[dirZERO])[kzero],
        (D.f[dirTNE ])[ktne],
        (D.f[dirTSW ])[ktsw],
        (D.f[dirTSE ])[ktse],
        (D.f[dirTNW ])[ktnw],
        (D.f[dirBNE ])[kbne],
        (D.f[dirBSW ])[kbsw],
        (D.f[dirBSE ])[kbse],
        (D.f[dirBNW ])[kbnw]
    };
}

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMac27( real* vxD,
                                        real* vyD,
                                        real* vzD,
                                        real* rhoD,
                                        unsigned int* geoD,
                                        unsigned int* neighborX,
                                        unsigned int* neighborY,
                                        unsigned int* neighborZ,
                                        unsigned int size_Mat,
                                        real* DD,
                                        bool evenOrOdd)
{
   const unsigned int tx = threadIdx.x;    // Thread index = lokaler i index
   const unsigned int by = blockIdx.x;     // Block index x
   const unsigned int bz = blockIdx.y;     // Block index y
   const unsigned int x = tx + STARTOFFX;  // Globaler x-Index 
   const unsigned int y = by + STARTOFFY;  // Globaler y-Index 
   const unsigned int z = bz + STARTOFFZ;  // Globaler z-Index 

   const unsigned sizeX = blockDim.x;
   const unsigned sizeY = gridDim.x;
   const unsigned nx = sizeX + 2 * STARTOFFX;
   const unsigned ny = sizeY + 2 * STARTOFFY;

   const unsigned int k = nx*(ny*z + y) + x; // Zugriff auf arrays im device

   rhoD[k] = c0o1;
   vxD[k]  = c0o1;
   vyD[k]  = c0o1;
   vzD[k]  = c0o1;

   if(geoD[k] == GEO_FLUID)
   {
       const auto distribution = getDistribution(DD, size_Mat, evenOrOdd, k, neighborX, neighborY, neighborZ);

       rhoD[k] = LBM::getDensity(distribution.f);
       vxD[k] = LBM::getIncompVelocityX1(distribution.f);
       vyD[k] = LBM::getIncompVelocityX2(distribution.f);
       vzD[k] = LBM::getIncompVelocityX3(distribution.f);
   }
}





////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacSP27( real* vxD,
                                          real* vyD,
                                          real* vzD,
                                          real* rhoD,
                                          real* pressD,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD,
                                          bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
   {
      D.f[dirE   ] = &DD[dirE   *size_Mat];
      D.f[dirW   ] = &DD[dirW   *size_Mat];
      D.f[dirN   ] = &DD[dirN   *size_Mat];
      D.f[dirS   ] = &DD[dirS   *size_Mat];
      D.f[dirT   ] = &DD[dirT   *size_Mat];
      D.f[dirB   ] = &DD[dirB   *size_Mat];
      D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      D.f[dirBNW ] = &DD[dirBNW *size_Mat];
   } 
   else
   {
      D.f[dirW   ] = &DD[dirE   *size_Mat];
      D.f[dirE   ] = &DD[dirW   *size_Mat];
      D.f[dirS   ] = &DD[dirN   *size_Mat];
      D.f[dirN   ] = &DD[dirS   *size_Mat];
      D.f[dirB   ] = &DD[dirT   *size_Mat];
      D.f[dirT   ] = &DD[dirB   *size_Mat];
      D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirZERO])[kzero]+ 
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vxD[k]     =   (D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

         vyD[k]     =   (D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vzD[k]     =   (D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

         pressD[k]  =  ((D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        2.f*(
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ])+
                        3.f*(
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+c0o1*rhoD[k])) * c1o2+rhoD[k]; // times zero for incompressible case   
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5

      }
   }
}




























////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacCompSP27( real* vxD,
											  real* vyD,
											  real* vzD,
											  real* rhoD,
											  real* pressD,
											  unsigned int* geoD,
											  unsigned int* neighborX,
											  unsigned int* neighborY,
											  unsigned int* neighborZ,
											  unsigned int size_Mat,
											  real* DD,
											  bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
   {
      D.f[dirE   ] = &DD[dirE   *size_Mat];
      D.f[dirW   ] = &DD[dirW   *size_Mat];
      D.f[dirN   ] = &DD[dirN   *size_Mat];
      D.f[dirS   ] = &DD[dirS   *size_Mat];
      D.f[dirT   ] = &DD[dirT   *size_Mat];
      D.f[dirB   ] = &DD[dirB   *size_Mat];
      D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      D.f[dirBNW ] = &DD[dirBNW *size_Mat];
   } 
   else
   {
      D.f[dirW   ] = &DD[dirE   *size_Mat];
      D.f[dirE   ] = &DD[dirW   *size_Mat];
      D.f[dirS   ] = &DD[dirN   *size_Mat];
      D.f[dirN   ] = &DD[dirS   *size_Mat];
      D.f[dirB   ] = &DD[dirT   *size_Mat];
      D.f[dirT   ] = &DD[dirB   *size_Mat];
      D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID || geoD[k] == GEO_PM_0 || geoD[k] == GEO_PM_1 || geoD[k] == GEO_PM_2)
      {
         rhoD[k]    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirZERO])[kzero]+ 
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vxD[k]     =  ((D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
						(D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (c1o1 + rhoD[k]);

         vyD[k]     =  ((D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw]) / (c1o1 + rhoD[k]);

         vzD[k]     =  ((D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (c1o1 + rhoD[k]);

         pressD[k]  =  ((D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        c2o1*(
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ])+
                        c3o1*(
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]; // times zero for incompressible case   
         //achtung op hart gesetzt Annahme op = 1 ;                                                      ^^^^(1.0/op-0.5)=0.5

      }
   }
}



























////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacThS7( real* Conc,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD7,
                                          bool evenOrOdd)
{
   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   } 
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D7.f[1])[ke   ]+ (D7.f[2])[kw  ]+ 
                        (D7.f[3])[kn   ]+ (D7.f[4])[ks  ]+
                        (D7.f[5])[kt   ]+ (D7.f[6])[kb  ]+
                        (D7.f[0])[kzero];  
      }
   }
}





























////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void GetPlaneConcThS7(real* Conc,
								            int* kPC,
								            unsigned int numberOfPointskPC,
											unsigned int* geoD,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat,
											real* DD7,
											bool evenOrOdd)
{
   Distributions7 D7;
   if (evenOrOdd==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   } 
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfPointskPC)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kPC[k];
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D7.f[1])[ke   ]+ (D7.f[2])[kw  ]+ 
                        (D7.f[3])[kn   ]+ (D7.f[4])[ks  ]+
                        (D7.f[5])[kt   ]+ (D7.f[6])[kb  ]+
                        (D7.f[0])[kzero];  
      }
   }
}




































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void GetPlaneConcThS27(real* Conc,
								             int* kPC,
								             unsigned int numberOfPointskPC,
											 unsigned int* geoD,
											 unsigned int* neighborX,
											 unsigned int* neighborY,
											 unsigned int* neighborZ,
											 unsigned int size_Mat,
											 real* DD27,
											 bool evenOrOdd)
{
   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   }
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfPointskPC)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kPC[k];
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = kzero;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = kzero;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = kzero;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = kzero;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D27.f[dirE   ])[ke  ]+ (D27.f[dirW   ])[kw  ]+ 
                        (D27.f[dirN   ])[kn  ]+ (D27.f[dirS   ])[ks  ]+
                        (D27.f[dirT   ])[kt  ]+ (D27.f[dirB   ])[kb  ]+
                        (D27.f[dirNE  ])[kne ]+ (D27.f[dirSW  ])[ksw ]+
                        (D27.f[dirSE  ])[kse ]+ (D27.f[dirNW  ])[knw ]+
                        (D27.f[dirTE  ])[kte ]+ (D27.f[dirBW  ])[kbw ]+
                        (D27.f[dirBE  ])[kbe ]+ (D27.f[dirTW  ])[ktw ]+
                        (D27.f[dirTN  ])[ktn ]+ (D27.f[dirBS  ])[kbs ]+
                        (D27.f[dirBN  ])[kbn ]+ (D27.f[dirTS  ])[kts ]+
                        (D27.f[dirZERO])[kzero]+ 
                        (D27.f[dirTNE ])[ktne]+ (D27.f[dirTSW ])[ktsw]+
                        (D27.f[dirTSE ])[ktse]+ (D27.f[dirTNW ])[ktnw]+
                        (D27.f[dirBNE ])[kbne]+ (D27.f[dirBSW ])[kbsw]+
                        (D27.f[dirBSE ])[kbse]+ (D27.f[dirBNW ])[kbnw];
      }
   }   
}




































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacThS27(real* Conc,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD27,
                                          bool evenOrOdd)
{
   Distributions27 D27;
   if (evenOrOdd==true)
   {
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   }
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D27.f[dirE   ])[ke  ]+ (D27.f[dirW   ])[kw  ]+ 
                        (D27.f[dirN   ])[kn  ]+ (D27.f[dirS   ])[ks  ]+
                        (D27.f[dirT   ])[kt  ]+ (D27.f[dirB   ])[kb  ]+
                        (D27.f[dirNE  ])[kne ]+ (D27.f[dirSW  ])[ksw ]+
                        (D27.f[dirSE  ])[kse ]+ (D27.f[dirNW  ])[knw ]+
                        (D27.f[dirTE  ])[kte ]+ (D27.f[dirBW  ])[kbw ]+
                        (D27.f[dirBE  ])[kbe ]+ (D27.f[dirTW  ])[ktw ]+
                        (D27.f[dirTN  ])[ktn ]+ (D27.f[dirBS  ])[kbs ]+
                        (D27.f[dirBN  ])[kbn ]+ (D27.f[dirTS  ])[kts ]+
                        (D27.f[dirZERO])[kzero]+ 
                        (D27.f[dirTNE ])[ktne]+ (D27.f[dirTSW ])[ktsw]+
                        (D27.f[dirTSE ])[ktse]+ (D27.f[dirTNW ])[ktnw]+
                        (D27.f[dirBNE ])[kbne]+ (D27.f[dirBSW ])[kbsw]+
                        (D27.f[dirBSE ])[kbse]+ (D27.f[dirBNW ])[kbnw];
      }
   }   
}




















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMedSP27( real* vxD,
                                          real* vyD,
                                          real* vzD,
                                          real* rhoD,
                                          real* pressD,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD,
                                          bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
   {
      D.f[dirE   ] = &DD[dirE   *size_Mat];
      D.f[dirW   ] = &DD[dirW   *size_Mat];
      D.f[dirN   ] = &DD[dirN   *size_Mat];
      D.f[dirS   ] = &DD[dirS   *size_Mat];
      D.f[dirT   ] = &DD[dirT   *size_Mat];
      D.f[dirB   ] = &DD[dirB   *size_Mat];
      D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      D.f[dirBNW ] = &DD[dirBNW *size_Mat];
   } 
   else
   {
      D.f[dirW   ] = &DD[dirE   *size_Mat];
      D.f[dirE   ] = &DD[dirW   *size_Mat];
      D.f[dirS   ] = &DD[dirN   *size_Mat];
      D.f[dirN   ] = &DD[dirS   *size_Mat];
      D.f[dirB   ] = &DD[dirT   *size_Mat];
      D.f[dirT   ] = &DD[dirB   *size_Mat];
      D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirZERO])[kzero]+ 
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw]+
                        RHO;

         vxD[k]     =   (D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]+
                        VX;

         vyD[k]     =   (D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw]+
                        VY;

         vzD[k]     =   (D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]+
                        VZ;

         pressD[k]  =   ((D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        c2o1*(
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ])+
                        c3o1*(
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMedCompSP27( real* vxD,
											  real* vyD,
											  real* vzD,
											  real* rhoD,
											  real* pressD,
											  unsigned int* geoD,
											  unsigned int* neighborX,
											  unsigned int* neighborY,
											  unsigned int* neighborZ,
											  unsigned int size_Mat,
											  real* DD,
											  bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
   {
      D.f[dirE   ] = &DD[dirE   *size_Mat];
      D.f[dirW   ] = &DD[dirW   *size_Mat];
      D.f[dirN   ] = &DD[dirN   *size_Mat];
      D.f[dirS   ] = &DD[dirS   *size_Mat];
      D.f[dirT   ] = &DD[dirT   *size_Mat];
      D.f[dirB   ] = &DD[dirB   *size_Mat];
      D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      D.f[dirBNW ] = &DD[dirBNW *size_Mat];
   } 
   else
   {
      D.f[dirW   ] = &DD[dirE   *size_Mat];
      D.f[dirE   ] = &DD[dirW   *size_Mat];
      D.f[dirS   ] = &DD[dirN   *size_Mat];
      D.f[dirN   ] = &DD[dirS   *size_Mat];
      D.f[dirB   ] = &DD[dirT   *size_Mat];
      D.f[dirT   ] = &DD[dirB   *size_Mat];
      D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      D.f[dirZERO] = &DD[dirZERO*size_Mat];
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      //unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
	  rhoD[k]   = c0o1;
	  vxD[k]    = c0o1;
	  vyD[k]    = c0o1;
	  vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
		  real mfcbb = (D.f[dirE])[k];//[ke   ];
		  real mfabb = (D.f[dirW])[kw];//[kw   ];  
		  real mfbcb = (D.f[dirN])[k];//[kn   ];
		  real mfbab = (D.f[dirS])[ks];//[ks   ];  
		  real mfbbc = (D.f[dirT])[k];//[kt   ];
		  real mfbba = (D.f[dirB])[kb];//[kb   ];  
		  real mfccb = (D.f[dirNE])[k];//[kne  ];  
		  real mfaab = (D.f[dirSW])[ksw];//[ksw  ];
		  real mfcab = (D.f[dirSE])[ks];//[kse  ]; 
		  real mfacb = (D.f[dirNW])[kw];//[knw  ]; 
		  real mfcbc = (D.f[dirTE])[k];//[kte  ];  
		  real mfaba = (D.f[dirBW])[kbw];//[kbw  ];
		  real mfcba = (D.f[dirBE])[kb];//[kbe  ]; 
		  real mfabc = (D.f[dirTW])[kw];//[ktw  ]; 
		  real mfbcc = (D.f[dirTN])[k];//[ktn  ];  
		  real mfbaa = (D.f[dirBS])[kbs];//[kbs  ];
		  real mfbca = (D.f[dirBN])[kb];//[kbn  ]; 
		  real mfbac = (D.f[dirTS])[ks];//[kts  ]; 
		  real mfbbb = (D.f[dirZERO])[k];//[kzero];
		  real mfccc = (D.f[dirTNE])[k];//[ktne ]; 
		  real mfaac = (D.f[dirTSW])[ksw];//[ktsw ]; 
		  real mfcac = (D.f[dirTSE])[ks];//[ktse ];
		  real mfacc = (D.f[dirTNW])[kw];//[ktnw ];
		  real mfcca = (D.f[dirBNE])[kb];//[kbne ];
		  real mfaaa = (D.f[dirBSW])[kbsw];//[kbsw ];
		  real mfcaa = (D.f[dirBSE])[kbs];//[kbse ]; 
		  real mfaca = (D.f[dirBNW])[kbw];//[kbnw ]; 
		  ////////////////////////////////////////////////////////////////////////////////////
		  real drho = 
			  ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
			  (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
			  ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

		  real rho = c1o1 + drho;
		  
		  rhoD[k] = drho + RHO;

		  vxD[k] = 
			  (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
			  (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
			  (mfcbb - mfabb)) / rho) + VX;
		  vyD[k] = 
			  (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
			  (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
			  (mfbcb - mfbab)) / rho) + VY;
		  vzD[k] = 
			  (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
			  (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
			  (mfbbc - mfbba)) / rho) + VZ;

		  //rhoD[k] =
			 // (D.f[dirE])[ke] + (D.f[dirW])[kw] +
			 // (D.f[dirN])[kn] + (D.f[dirS])[ks] +
			 // (D.f[dirT])[kt] + (D.f[dirB])[kb] +
			 // (D.f[dirNE])[kne] + (D.f[dirSW])[ksw] +
			 // (D.f[dirSE])[kse] + (D.f[dirNW])[knw] +
			 // (D.f[dirTE])[kte] + (D.f[dirBW])[kbw] +
			 // (D.f[dirBE])[kbe] + (D.f[dirTW])[ktw] +
			 // (D.f[dirTN])[ktn] + (D.f[dirBS])[kbs] +
			 // (D.f[dirBN])[kbn] + (D.f[dirTS])[kts] +
			 // (D.f[dirZERO])[kzero] +
			 // (D.f[dirTNE])[ktne] + (D.f[dirTSW])[ktsw] +
			 // (D.f[dirTSE])[ktse] + (D.f[dirTNW])[ktnw] +
			 // (D.f[dirBNE])[kbne] + (D.f[dirBSW])[kbsw] +
			 // (D.f[dirBSE])[kbse] + (D.f[dirBNW])[kbnw];// +RHO;

    //     vxD[k] =  
			 //((D.f[dirE  ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
    //         (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
    //         (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
    //         (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
    //         (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
    //         (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
    //         (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
    //         (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
    //         (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (one + rhoD[k])+
    //         VX;

    //     vyD[k] =  
			 //((D.f[dirN  ])[kn  ]- (D.f[dirS   ])[ks  ]+
    //         (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
    //         (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
    //         (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
    //         (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
    //         (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
    //         (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
    //         (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
    //         (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw]) / (one + rhoD[k])+
    //         VY;

    //     vzD[k] =  
			 //((D.f[dirT  ])[kt  ]- (D.f[dirB   ])[kb  ]+
    //         (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
    //         (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
    //         (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
    //         (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
    //         (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
    //         (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
    //         (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
    //         (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (one + rhoD[k])+
    //         VZ;

         pressD[k]  =  ((D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        c2o1*(
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ])+
                        c3o1*(
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1+rhoD[k])) * c1o2+rhoD[k]+
                        PRESS;    
         //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMedCompAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int* geoD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int size_Mat,
	real* DD,
	real* DD_AD,
	bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd == true)
	{
		D.f[dirE] = &DD[dirE   *size_Mat];
		D.f[dirW] = &DD[dirW   *size_Mat];
		D.f[dirN] = &DD[dirN   *size_Mat];
		D.f[dirS] = &DD[dirS   *size_Mat];
		D.f[dirT] = &DD[dirT   *size_Mat];
		D.f[dirB] = &DD[dirB   *size_Mat];
		D.f[dirNE] = &DD[dirNE  *size_Mat];
		D.f[dirSW] = &DD[dirSW  *size_Mat];
		D.f[dirSE] = &DD[dirSE  *size_Mat];
		D.f[dirNW] = &DD[dirNW  *size_Mat];
		D.f[dirTE] = &DD[dirTE  *size_Mat];
		D.f[dirBW] = &DD[dirBW  *size_Mat];
		D.f[dirBE] = &DD[dirBE  *size_Mat];
		D.f[dirTW] = &DD[dirTW  *size_Mat];
		D.f[dirTN] = &DD[dirTN  *size_Mat];
		D.f[dirBS] = &DD[dirBS  *size_Mat];
		D.f[dirBN] = &DD[dirBN  *size_Mat];
		D.f[dirTS] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE] = &DD[dirTNE *size_Mat];
		D.f[dirTSW] = &DD[dirTSW *size_Mat];
		D.f[dirTSE] = &DD[dirTSE *size_Mat];
		D.f[dirTNW] = &DD[dirTNW *size_Mat];
		D.f[dirBNE] = &DD[dirBNE *size_Mat];
		D.f[dirBSW] = &DD[dirBSW *size_Mat];
		D.f[dirBSE] = &DD[dirBSE *size_Mat];
		D.f[dirBNW] = &DD[dirBNW *size_Mat];
	}
	else
	{
		D.f[dirW] = &DD[dirE   *size_Mat];
		D.f[dirE] = &DD[dirW   *size_Mat];
		D.f[dirS] = &DD[dirN   *size_Mat];
		D.f[dirN] = &DD[dirS   *size_Mat];
		D.f[dirB] = &DD[dirT   *size_Mat];
		D.f[dirT] = &DD[dirB   *size_Mat];
		D.f[dirSW] = &DD[dirNE  *size_Mat];
		D.f[dirNE] = &DD[dirSW  *size_Mat];
		D.f[dirNW] = &DD[dirSE  *size_Mat];
		D.f[dirSE] = &DD[dirNW  *size_Mat];
		D.f[dirBW] = &DD[dirTE  *size_Mat];
		D.f[dirTE] = &DD[dirBW  *size_Mat];
		D.f[dirTW] = &DD[dirBE  *size_Mat];
		D.f[dirBE] = &DD[dirTW  *size_Mat];
		D.f[dirBS] = &DD[dirTN  *size_Mat];
		D.f[dirTN] = &DD[dirBS  *size_Mat];
		D.f[dirTS] = &DD[dirBN  *size_Mat];
		D.f[dirBN] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE] = &DD[dirBSW *size_Mat];
		D.f[dirTSW] = &DD[dirBNE *size_Mat];
		D.f[dirTSE] = &DD[dirBNW *size_Mat];
		D.f[dirTNW] = &DD[dirBSE *size_Mat];
		D.f[dirBNE] = &DD[dirTSW *size_Mat];
		D.f[dirBSW] = &DD[dirTNE *size_Mat];
		D.f[dirBSE] = &DD[dirTNW *size_Mat];
		D.f[dirBNW] = &DD[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	Distributions27 Dad;
	if (evenOrOdd == true)
	{
		Dad.f[dirE]    = &DD_AD[dirE   *size_Mat];
		Dad.f[dirW]    = &DD_AD[dirW   *size_Mat];
		Dad.f[dirN]    = &DD_AD[dirN   *size_Mat];
		Dad.f[dirS]    = &DD_AD[dirS   *size_Mat];
		Dad.f[dirT]    = &DD_AD[dirT   *size_Mat];
		Dad.f[dirB]    = &DD_AD[dirB   *size_Mat];
		Dad.f[dirNE]   = &DD_AD[dirNE  *size_Mat];
		Dad.f[dirSW]   = &DD_AD[dirSW  *size_Mat];
		Dad.f[dirSE]   = &DD_AD[dirSE  *size_Mat];
		Dad.f[dirNW]   = &DD_AD[dirNW  *size_Mat];
		Dad.f[dirTE]   = &DD_AD[dirTE  *size_Mat];
		Dad.f[dirBW]   = &DD_AD[dirBW  *size_Mat];
		Dad.f[dirBE]   = &DD_AD[dirBE  *size_Mat];
		Dad.f[dirTW]   = &DD_AD[dirTW  *size_Mat];
		Dad.f[dirTN]   = &DD_AD[dirTN  *size_Mat];
		Dad.f[dirBS]   = &DD_AD[dirBS  *size_Mat];
		Dad.f[dirBN]   = &DD_AD[dirBN  *size_Mat];
		Dad.f[dirTS]   = &DD_AD[dirTS  *size_Mat];
		Dad.f[dirZERO] = &DD_AD[dirZERO*size_Mat];
		Dad.f[dirTNE]  = &DD_AD[dirTNE *size_Mat];
		Dad.f[dirTSW]  = &DD_AD[dirTSW *size_Mat];
		Dad.f[dirTSE]  = &DD_AD[dirTSE *size_Mat];
		Dad.f[dirTNW]  = &DD_AD[dirTNW *size_Mat];
		Dad.f[dirBNE]  = &DD_AD[dirBNE *size_Mat];
		Dad.f[dirBSW]  = &DD_AD[dirBSW *size_Mat];
		Dad.f[dirBSE]  = &DD_AD[dirBSE *size_Mat];
		Dad.f[dirBNW]  = &DD_AD[dirBNW *size_Mat];
	}						
	else					
	{						
		Dad.f[dirW]    = &DD_AD[dirE   *size_Mat];
		Dad.f[dirE]    = &DD_AD[dirW   *size_Mat];
		Dad.f[dirS]    = &DD_AD[dirN   *size_Mat];
		Dad.f[dirN]    = &DD_AD[dirS   *size_Mat];
		Dad.f[dirB]    = &DD_AD[dirT   *size_Mat];
		Dad.f[dirT]    = &DD_AD[dirB   *size_Mat];
		Dad.f[dirSW]   = &DD_AD[dirNE  *size_Mat];
		Dad.f[dirNE]   = &DD_AD[dirSW  *size_Mat];
		Dad.f[dirNW]   = &DD_AD[dirSE  *size_Mat];
		Dad.f[dirSE]   = &DD_AD[dirNW  *size_Mat];
		Dad.f[dirBW]   = &DD_AD[dirTE  *size_Mat];
		Dad.f[dirTE]   = &DD_AD[dirBW  *size_Mat];
		Dad.f[dirTW]   = &DD_AD[dirBE  *size_Mat];
		Dad.f[dirBE]   = &DD_AD[dirTW  *size_Mat];
		Dad.f[dirBS]   = &DD_AD[dirTN  *size_Mat];
		Dad.f[dirTN]   = &DD_AD[dirBS  *size_Mat];
		Dad.f[dirTS]   = &DD_AD[dirBN  *size_Mat];
		Dad.f[dirBN]   = &DD_AD[dirTS  *size_Mat];
		Dad.f[dirZERO] = &DD_AD[dirZERO*size_Mat];
		Dad.f[dirTNE]  = &DD_AD[dirBSW *size_Mat];
		Dad.f[dirTSW]  = &DD_AD[dirBNE *size_Mat];
		Dad.f[dirTSE]  = &DD_AD[dirBNW *size_Mat];
		Dad.f[dirTNW]  = &DD_AD[dirBSE *size_Mat];
		Dad.f[dirBNE]  = &DD_AD[dirTSW *size_Mat];
		Dad.f[dirBSW]  = &DD_AD[dirTNE *size_Mat];
		Dad.f[dirBSE]  = &DD_AD[dirTNW *size_Mat];
		Dad.f[dirBNW]  = &DD_AD[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		//index
		//unsigned int kzero = k;
		unsigned int ke = k;
		unsigned int kw = neighborX[k];
		unsigned int kn = k;
		unsigned int ks = neighborY[k];
		unsigned int kt = k;
		unsigned int kb = neighborZ[k];
		unsigned int ksw = neighborY[kw];
		unsigned int kne = k;
		unsigned int kse = ks;
		unsigned int knw = kw;
		unsigned int kbw = neighborZ[kw];
		unsigned int kte = k;
		unsigned int kbe = kb;
		unsigned int ktw = kw;
		unsigned int kbs = neighborZ[ks];
		unsigned int ktn = k;
		unsigned int kbn = kb;
		unsigned int kts = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = k;
		unsigned int kbsw = neighborZ[ksw];
		//////////////////////////////////////////////////////////////////////////
		real CONC  = concD[k];
		real PRESS = pressD[k];
		real RHO   = rhoD[k];
		real VX    = vxD[k];
		real VY    = vyD[k];
		real VZ    = vzD[k];
		//////////////////////////////////////////////////////////////////////////
		concD[k] = c0o1;
		pressD[k] = c0o1;
		rhoD[k] = c0o1;
		vxD[k] = c0o1;
		vyD[k] = c0o1;
		vzD[k] = c0o1;

		if (geoD[k] == GEO_FLUID)
		{
			real mfcbb = (D.f[dirE])[k];//[ke   ];
			real mfabb = (D.f[dirW])[kw];//[kw   ];  
			real mfbcb = (D.f[dirN])[k];//[kn   ];
			real mfbab = (D.f[dirS])[ks];//[ks   ];  
			real mfbbc = (D.f[dirT])[k];//[kt   ];
			real mfbba = (D.f[dirB])[kb];//[kb   ];  
			real mfccb = (D.f[dirNE])[k];//[kne  ];  
			real mfaab = (D.f[dirSW])[ksw];//[ksw  ];
			real mfcab = (D.f[dirSE])[ks];//[kse  ]; 
			real mfacb = (D.f[dirNW])[kw];//[knw  ]; 
			real mfcbc = (D.f[dirTE])[k];//[kte  ];  
			real mfaba = (D.f[dirBW])[kbw];//[kbw  ];
			real mfcba = (D.f[dirBE])[kb];//[kbe  ]; 
			real mfabc = (D.f[dirTW])[kw];//[ktw  ]; 
			real mfbcc = (D.f[dirTN])[k];//[ktn  ];  
			real mfbaa = (D.f[dirBS])[kbs];//[kbs  ];
			real mfbca = (D.f[dirBN])[kb];//[kbn  ]; 
			real mfbac = (D.f[dirTS])[ks];//[kts  ]; 
			real mfbbb = (D.f[dirZERO])[k];//[kzero];
			real mfccc = (D.f[dirTNE])[k];//[ktne ]; 
			real mfaac = (D.f[dirTSW])[ksw];//[ktsw ]; 
			real mfcac = (D.f[dirTSE])[ks];//[ktse ];
			real mfacc = (D.f[dirTNW])[kw];//[ktnw ];
			real mfcca = (D.f[dirBNE])[kb];//[kbne ];
			real mfaaa = (D.f[dirBSW])[kbsw];//[kbsw ];
			real mfcaa = (D.f[dirBSE])[kbs];//[kbse ]; 
			real mfaca = (D.f[dirBNW])[kbw];//[kbnw ]; 
			////////////////////////////////////////////////////////////////////////////////////
			real drho =
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) + mfbbb;
			real rho = c1o1 + drho;
			////////////////////////////////////////////////////////////////////////////////////

			rhoD[k] = drho + RHO;

			vxD[k] =
				(((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					(mfcbb - mfabb)) / rho) + VX;
			
			vyD[k] =
				(((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					(mfbcb - mfbab)) / rho) + VY;
			
			vzD[k] =
				(((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					(mfbbc - mfbba)) / rho) + VZ;

			pressD[k] = 
				((D.f[dirE])[ke] + (D.f[dirW])[kw] +
				 (D.f[dirN])[kn] + (D.f[dirS])[ks] +
				 (D.f[dirT])[kt] + (D.f[dirB])[kb] +
				 c2o1*(
				 (D.f[dirNE])[kne] + (D.f[dirSW])[ksw] +
				 (D.f[dirSE])[kse] + (D.f[dirNW])[knw] +
				 (D.f[dirTE])[kte] + (D.f[dirBW])[kbw] +
				 (D.f[dirBE])[kbe] + (D.f[dirTW])[ktw] +
				 (D.f[dirTN])[ktn] + (D.f[dirBS])[kbs] +
				 (D.f[dirBN])[kbn] + (D.f[dirTS])[kts]) +
				 c3o1*(
				 (D.f[dirTNE])[ktne] + (D.f[dirTSW])[ktsw] +
				 (D.f[dirTSE])[ktse] + (D.f[dirTNW])[ktnw] +
				 (D.f[dirBNE])[kbne] + (D.f[dirBSW])[kbsw] +
				 (D.f[dirBSE])[kbse] + (D.f[dirBNW])[kbnw]) -
				 rhoD[k] - (vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (c1o1 + rhoD[k])) * c1o2 + rhoD[k] +
				 PRESS;
				 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
			//////////////////////////////////////////////////////////////////////////
			mfcbb = (Dad.f[dirE   ])[k   ];
			mfabb = (Dad.f[dirW   ])[kw  ];
			mfbcb = (Dad.f[dirN   ])[k   ];
			mfbab = (Dad.f[dirS   ])[ks  ];
			mfbbc = (Dad.f[dirT   ])[k   ];
			mfbba = (Dad.f[dirB   ])[kb  ];
			mfccb = (Dad.f[dirNE  ])[k   ];
			mfaab = (Dad.f[dirSW  ])[ksw ];
			mfcab = (Dad.f[dirSE  ])[ks  ];
			mfacb = (Dad.f[dirNW  ])[kw  ];
			mfcbc = (Dad.f[dirTE  ])[k   ];
			mfaba = (Dad.f[dirBW  ])[kbw ];
			mfcba = (Dad.f[dirBE  ])[kb  ];
			mfabc = (Dad.f[dirTW  ])[kw  ];
			mfbcc = (Dad.f[dirTN  ])[k   ];
			mfbaa = (Dad.f[dirBS  ])[kbs ];
			mfbca = (Dad.f[dirBN  ])[kb  ];
			mfbac = (Dad.f[dirTS  ])[ks  ];
			mfbbb = (Dad.f[dirZERO])[k   ];
			mfccc = (Dad.f[dirTNE ])[k   ];
			mfaac = (Dad.f[dirTSW ])[ksw ];
			mfcac = (Dad.f[dirTSE ])[ks  ];
			mfacc = (Dad.f[dirTNW ])[kw  ];
			mfcca = (Dad.f[dirBNE ])[kb  ];
			mfaaa = (Dad.f[dirBSW ])[kbsw];
			mfcaa = (Dad.f[dirBSE ])[kbs ];
			mfaca = (Dad.f[dirBNW ])[kbw ];
			//////////////////////////////////////////////////////////////////////////
			concD[k] = 
				((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
				 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb + CONC;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacMedSP27( real* vxD,
                                             real* vyD,
                                             real* vzD,
                                             real* rhoD,
                                             real* pressD,
                                             unsigned int* geoD,
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int tdiff,
                                             unsigned int size_Mat,
                                             bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      real PRESS = pressD[k];
      real RHO   = rhoD[k];
      real VX    = vxD[k];
      real VY    = vyD[k];
      real VZ    = vzD[k];
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = c0o1;
      rhoD[k]   = c0o1;
      vxD[k]    = c0o1;
      vyD[k]    = c0o1;
      vzD[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   RHO   / tdiff;
         vxD[k]     =   VX    / tdiff;
         vyD[k]     =   VY    / tdiff;
         vzD[k]     =   VZ    / tdiff;
         pressD[k]  =   PRESS / tdiff;    
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBResetMedianValuesSP27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	unsigned int size_Mat,
	bool evenOrOdd)
{
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k<size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		pressD[k] = c0o1;
		rhoD[k] = c0o1;
		vxD[k] = c0o1;
		vyD[k] = c0o1;
		vzD[k] = c0o1;
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBResetMedianValuesAD27(
	real* vxD,
	real* vyD,
	real* vzD,
	real* rhoD,
	real* pressD,
	real* concD,
	unsigned int size_Mat,
	bool evenOrOdd)
{
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < size_Mat)
	{
		//////////////////////////////////////////////////////////////////////////
		concD[k]  = c0o1;
		pressD[k] = c0o1;
		rhoD[k]   = c0o1;
		vxD[k]    = c0o1;
		vyD[k]    = c0o1;
		vzD[k]    = c0o1;
	}
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMeasurePoints( real* vxMP,
												real* vyMP,
												real* vzMP,
												real* rhoMP,
												unsigned int* kMP,
												unsigned int numberOfPointskMP,
												unsigned int MPClockCycle,
												unsigned int t,
												unsigned int* geoD,
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned int size_Mat,
												real* DD,
												bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd==true)
	{
		D.f[dirE   ] = &DD[dirE   *size_Mat];
		D.f[dirW   ] = &DD[dirW   *size_Mat];
		D.f[dirN   ] = &DD[dirN   *size_Mat];
		D.f[dirS   ] = &DD[dirS   *size_Mat];
		D.f[dirT   ] = &DD[dirT   *size_Mat];
		D.f[dirB   ] = &DD[dirB   *size_Mat];
		D.f[dirNE  ] = &DD[dirNE  *size_Mat];
		D.f[dirSW  ] = &DD[dirSW  *size_Mat];
		D.f[dirSE  ] = &DD[dirSE  *size_Mat];
		D.f[dirNW  ] = &DD[dirNW  *size_Mat];
		D.f[dirTE  ] = &DD[dirTE  *size_Mat];
		D.f[dirBW  ] = &DD[dirBW  *size_Mat];
		D.f[dirBE  ] = &DD[dirBE  *size_Mat];
		D.f[dirTW  ] = &DD[dirTW  *size_Mat];
		D.f[dirTN  ] = &DD[dirTN  *size_Mat];
		D.f[dirBS  ] = &DD[dirBS  *size_Mat];
		D.f[dirBN  ] = &DD[dirBN  *size_Mat];
		D.f[dirTS  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirTNE *size_Mat];
		D.f[dirTSW ] = &DD[dirTSW *size_Mat];
		D.f[dirTSE ] = &DD[dirTSE *size_Mat];
		D.f[dirTNW ] = &DD[dirTNW *size_Mat];
		D.f[dirBNE ] = &DD[dirBNE *size_Mat];
		D.f[dirBSW ] = &DD[dirBSW *size_Mat];
		D.f[dirBSE ] = &DD[dirBSE *size_Mat];
		D.f[dirBNW ] = &DD[dirBNW *size_Mat];
	} 
	else
	{
		D.f[dirW   ] = &DD[dirE   *size_Mat];
		D.f[dirE   ] = &DD[dirW   *size_Mat];
		D.f[dirS   ] = &DD[dirN   *size_Mat];
		D.f[dirN   ] = &DD[dirS   *size_Mat];
		D.f[dirB   ] = &DD[dirT   *size_Mat];
		D.f[dirT   ] = &DD[dirB   *size_Mat];
		D.f[dirSW  ] = &DD[dirNE  *size_Mat];
		D.f[dirNE  ] = &DD[dirSW  *size_Mat];
		D.f[dirNW  ] = &DD[dirSE  *size_Mat];
		D.f[dirSE  ] = &DD[dirNW  *size_Mat];
		D.f[dirBW  ] = &DD[dirTE  *size_Mat];
		D.f[dirTE  ] = &DD[dirBW  *size_Mat];
		D.f[dirTW  ] = &DD[dirBE  *size_Mat];
		D.f[dirBE  ] = &DD[dirTW  *size_Mat];
		D.f[dirBS  ] = &DD[dirTN  *size_Mat];
		D.f[dirTN  ] = &DD[dirBS  *size_Mat];
		D.f[dirTS  ] = &DD[dirBN  *size_Mat];
		D.f[dirBN  ] = &DD[dirTS  *size_Mat];
		D.f[dirZERO] = &DD[dirZERO*size_Mat];
		D.f[dirTNE ] = &DD[dirBSW *size_Mat];
		D.f[dirTSW ] = &DD[dirBNE *size_Mat];
		D.f[dirTSE ] = &DD[dirBNW *size_Mat];
		D.f[dirTNW ] = &DD[dirBSE *size_Mat];
		D.f[dirBNE ] = &DD[dirTSW *size_Mat];
		D.f[dirBSW ] = &DD[dirTNE *size_Mat];
		D.f[dirBSE ] = &DD[dirTNW *size_Mat];
		D.f[dirBNW ] = &DD[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<numberOfPointskMP)
	{
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kMP[k];//k;
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = kzero;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = kzero;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = kzero;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = kzero;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
	  unsigned int kMac = k*MPClockCycle + t;
	  //////////////////////////////////////////////////////////////////////////

      if(geoD[kzero] == GEO_FLUID)
      {
         rhoMP[kMac]=   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirZERO])[kzero]+ 
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vxMP[kMac] =   (D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];

         vyMP[kMac] =   (D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vzMP[kMac] =   (D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw];
      }
   }
}
////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBSetOutputWallVelocitySP27( real* vxD,
														real* vyD,
														real* vzD,
														real* vxWall,
														real* vyWall,
														real* vzWall,
														int numberOfWallNodes, 
														int* kWallNodes, 
														real* rhoD,
														real* pressD,
														unsigned int* geoD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														unsigned int size_Mat,
														real* DD,
														bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfWallNodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KWN  = kWallNodes[k];
      //////////////////////////////////////////////////////////////////////////
      vxD[KWN] = 0.0;//vxWall[k];
      vyD[KWN] = 0.0;//vyWall[k];
      vzD[KWN] = 0.0;//vzWall[k];
   }
}





























