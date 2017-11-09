/* Device code */
#include "LBM/D3Q27.h"
//#include "LBM/LB.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitF3(unsigned int* neighborX,
									unsigned int* neighborY,
									unsigned int* neighborZ,
									unsigned int* geoD,
									doubflo* rho,
									doubflo* ux,
									doubflo* uy,
									doubflo* uz,
									unsigned int size_Mat,
									doubflo* G6,
									bool EvenOrOdd)
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
		////////////////////////////////////////////////////////////////////////////////
		unsigned int BC;
		BC = geoD[k];

		if (BC != GEO_SOLID &&  BC != GEO_VOID)
		{
			Distributions6 D;
			if (EvenOrOdd == true)
			{
				D.g[dirE] = &G6[dirE   *size_Mat];
				D.g[dirW] = &G6[dirW   *size_Mat];
				D.g[dirN] = &G6[dirN   *size_Mat];
				D.g[dirS] = &G6[dirS   *size_Mat];
				D.g[dirT] = &G6[dirT   *size_Mat];
				D.g[dirB] = &G6[dirB   *size_Mat];
			}
			else
			{
				D.g[dirW] = &G6[dirE   *size_Mat];
				D.g[dirE] = &G6[dirW   *size_Mat];
				D.g[dirS] = &G6[dirN   *size_Mat];
				D.g[dirN] = &G6[dirS   *size_Mat];
				D.g[dirB] = &G6[dirT   *size_Mat];
				D.g[dirT] = &G6[dirB   *size_Mat];
			}
			//////////////////////////////////////////////////////////////////////////
			//index
			//////////////////////////////////////////////////////////////////////////
			unsigned int kzero = k;
			unsigned int ke = k;
			unsigned int kw = neighborX[k];
			unsigned int kn = k;
			unsigned int ks = neighborY[k];
			unsigned int kt = k;
			unsigned int kb = neighborZ[k];
			//////////////////////////////////////////////////////////////////////////

			(D.g[dirE])[ke] = 0.0f;
			(D.g[dirW])[kw] = 0.0f;
			(D.g[dirN])[kn] = 0.0f;
			(D.g[dirS])[ks] = 0.0f;
			(D.g[dirT])[kt] = 0.0f;
			(D.g[dirB])[kb] = 0.0f;
		}
	}
}







////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInit27( int myid,
                                     int numprocs,
                                     doubflo u0,
                                     unsigned int* geoD,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     doubflo* vParabel,
                                     unsigned int size_Mat,
                                     unsigned int grid_nx, 
                                     unsigned int grid_ny, 
                                     unsigned int grid_nz, 
                                     doubflo* DD,
                                     int lev,
                                     int maxlev)
{
   Distributions27 D;
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
   ////////////////////////////////////////////////////////////////////////////////
   unsigned int  k;                   // Zugriff auf arrays im device
   //
   unsigned int tx = threadIdx.x;     // Thread index = lokaler i index
   unsigned int by = blockIdx.x;      // Block index x
   unsigned int bz = blockIdx.y;      // Block index y
   unsigned int  x = tx + STARTOFFX;  // Globaler x-Index 
   unsigned int  y = by + STARTOFFY;  // Globaler y-Index 
   unsigned int  z = bz + STARTOFFZ;  // Globaler z-Index 

   const unsigned sizeX = blockDim.x;
   const unsigned sizeY = gridDim.x;
   const unsigned nx = sizeX + 2 * STARTOFFX;
   const unsigned ny = sizeY + 2 * STARTOFFY;

   k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////
   geoD[k] = GEO_FLUID;
   if (lev==0)
   {
      if( by == 0 || by == grid_ny-1 || tx == 0 || tx == grid_nx-1 )             
         geoD[k] = GEO_SOLID;
      else if( bz == grid_nz-1 && myid == numprocs - 1 && geoD[k] != GEO_SOLID )
         geoD[k] = GEO_PRESS;				 
      else if( bz == 0 && myid == 0 && geoD[k] != GEO_SOLID)
         geoD[k] = GEO_SOLID;//GEO_VELO;
   }
   else if (lev==maxlev-1)
   {
      unsigned int centerX = grid_nx / 2;
      unsigned int centerY = grid_ny / 2;
      unsigned int centerZ = grid_nz / 2;
      doubflo        radius  = grid_ny / 2.56;

      unsigned int distSq = (centerX-tx)*(centerX-tx)+(centerY-by)*(centerY-by)+(centerZ-bz)*(centerZ-bz);
      doubflo radiSq = radius*radius;

      if( distSq < radiSq)        geoD[k] = GEO_SOLID;
   }
   //////////////////////////////////////////////////////////////////////////
   doubflo drho = zero;
   doubflo  vx1 = zero;
   doubflo  vx2 = zero;
   doubflo  vx3 = u0;
   vParabel[k] = vx3;
   ////////////////////////////////////////////////////////////////////////////////
   //index
   unsigned int nxny = nx*ny;
   ////////////////////////////////////////////////////////////////////////////////
   //neighborX[k]      = k+1;
   //neighborY[k+1]    = k+nx+1;
   //neighborZ[k+1]    = k+nxny+1;
   //neighborY[k]      = k+nx;
   //neighborX[k+nx]   = k+nx+1;
   //neighborZ[k+nx]   = k+nx+nxny;
   //neighborZ[k]      = k+nxny;
   //neighborX[k+nxny] = k+nxny+1;
   //neighborY[k+nxny] = k+nxny+nx;
   ////////////////////////////////////////////////////////////////////////////////
   unsigned int kzero= k;
   unsigned int ke   = k;
   unsigned int kw   = k + 1;
   unsigned int kn   = k;
   unsigned int ks   = k + nx;
   unsigned int kt   = k;
   unsigned int kb   = k + nxny;
   unsigned int ksw  = k + nx + 1;
   unsigned int kne  = k;
   unsigned int kse  = k + nx;
   unsigned int knw  = k + 1;
   unsigned int kbw  = k + nxny + 1;
   unsigned int kte  = k;
   unsigned int kbe  = k + nxny;
   unsigned int ktw  = k + 1;
   unsigned int kbs  = k + nxny + nx;
   unsigned int ktn  = k;
   unsigned int kbn  = k + nxny;
   unsigned int kts  = k + nx;
   unsigned int ktse = k + nx;
   unsigned int kbnw = k + nxny + 1;
   unsigned int ktnw = k + 1;
   unsigned int kbse = k + nxny + nx;
   unsigned int ktsw = k + nx + 1;
   unsigned int kbne = k + nxny;
   unsigned int ktne = k;
   unsigned int kbsw = k + nxny + nx + 1;
   //////////////////////////////////////////////////////////////////////////

   doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

   (D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq);
   (D.f[dirE   ])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
   (D.f[dirW   ])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
   (D.f[dirN   ])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
   (D.f[dirS   ])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
   (D.f[dirT   ])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
   (D.f[dirB   ])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
   (D.f[dirNE  ])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
   (D.f[dirSW  ])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
   (D.f[dirSE  ])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
   (D.f[dirNW  ])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
   (D.f[dirTE  ])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
   (D.f[dirBW  ])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
   (D.f[dirBE  ])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
   (D.f[dirTW  ])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
   (D.f[dirTN  ])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
   (D.f[dirBS  ])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
   (D.f[dirBN  ])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
   (D.f[dirTS  ])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
   (D.f[dirTNE ])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
   (D.f[dirBSW ])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
   (D.f[dirBNE ])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
   (D.f[dirTSW ])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
   (D.f[dirTSE ])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
   (D.f[dirBNW ])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
   (D.f[dirBSE ])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
   (D.f[dirTNW ])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

}
////////////////////////////////////////////////////////////////////////////////










////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitSP27( unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int* geoD,
                                       doubflo* rho,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       unsigned int size_Mat,
                                       doubflo* DD,
                                       bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID &&  BC != GEO_VOID)
      {
         Distributions27 D;
         if (EvenOrOdd==true)
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
            D.f[dirBSW ] = &DD[dirTNE *size_Mat];
            D.f[dirBNE ] = &DD[dirTSW *size_Mat];
            D.f[dirBNW ] = &DD[dirTSE *size_Mat];
            D.f[dirBSE ] = &DD[dirTNW *size_Mat];
            D.f[dirTSW ] = &DD[dirBNE *size_Mat];
            D.f[dirTNE ] = &DD[dirBSW *size_Mat];
            D.f[dirTNW ] = &DD[dirBSE *size_Mat];
            D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         doubflo drho = rho[k];//0.0f;//
         doubflo  vx1 = ux[k]; //0.0f;//
         doubflo  vx2 = uy[k]; //0.0f;//
         doubflo  vx3 = uz[k]; //0.0f;//
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
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
         doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq);
         (D.f[dirE   ])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
         (D.f[dirW   ])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[dirN   ])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[dirS   ])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[dirT   ])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
         (D.f[dirB   ])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[dirNE  ])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[dirSW  ])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[dirSE  ])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[dirNW  ])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[dirTE  ])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[dirBW  ])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[dirBE  ])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[dirTW  ])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[dirTN  ])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[dirBS  ])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[dirBN  ])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[dirTS  ])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[dirTNE ])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[dirBSW ])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[dirBNE ])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[dirTSW ])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[dirTSE ])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[dirBNW ])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[dirBSE ])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[dirTNW ])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
	  else
	  {
		  //////////////////////////////////////////////////////////////////////////
		  Distributions27 D;
		  D.f[dirZERO] = &DD[dirZERO*size_Mat];
		  //////////////////////////////////////////////////////////////////////////
		  (D.f[dirZERO])[k] = ninetysix;
		  //////////////////////////////////////////////////////////////////////////
	  }
   }
}







////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitCompSP27( unsigned int* neighborX,
										   unsigned int* neighborY,
										   unsigned int* neighborZ,
										   unsigned int* geoD,
										   doubflo* rho,
										   doubflo* ux,
										   doubflo* uy,
										   doubflo* uz,
										   unsigned int size_Mat,
										   doubflo* DD,
										   bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID &&  BC != GEO_VOID)
      {
         Distributions27 D;
         if (EvenOrOdd==true)
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
            D.f[dirBSW ] = &DD[dirTNE *size_Mat];
            D.f[dirBNE ] = &DD[dirTSW *size_Mat];
            D.f[dirBNW ] = &DD[dirTSE *size_Mat];
            D.f[dirBSE ] = &DD[dirTNW *size_Mat];
            D.f[dirTSW ] = &DD[dirBNE *size_Mat];
            D.f[dirTNE ] = &DD[dirBSW *size_Mat];
            D.f[dirTNW ] = &DD[dirBSE *size_Mat];
            D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         doubflo drho = rho[k];//0.0f;//
         doubflo  vx1 = ux[k]; //0.0f;//
         doubflo  vx2 = uy[k]; //0.0f;//
         doubflo  vx3 = uz[k]; //0.0f;//
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
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
         doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq*(one+drho));
         (D.f[dirE   ])[ke   ] =   c2over27* (drho+ (one+drho) * (three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq));
         (D.f[dirW   ])[kw   ] =   c2over27* (drho+ (one+drho) * (three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq));
         (D.f[dirN   ])[kn   ] =   c2over27* (drho+ (one+drho) * (three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq));
         (D.f[dirS   ])[ks   ] =   c2over27* (drho+ (one+drho) * (three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq));
         (D.f[dirT   ])[kt   ] =   c2over27* (drho+ (one+drho) * (three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq));
         (D.f[dirB   ])[kb   ] =   c2over27* (drho+ (one+drho) * (three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq));
         (D.f[dirNE  ])[kne  ] =   c1over54* (drho+ (one+drho) * (three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
         (D.f[dirSW  ])[ksw  ] =   c1over54* (drho+ (one+drho) * (three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
         (D.f[dirSE  ])[kse  ] =   c1over54* (drho+ (one+drho) * (three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
         (D.f[dirNW  ])[knw  ] =   c1over54* (drho+ (one+drho) * (three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
         (D.f[dirTE  ])[kte  ] =   c1over54* (drho+ (one+drho) * (three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
         (D.f[dirBW  ])[kbw  ] =   c1over54* (drho+ (one+drho) * (three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
         (D.f[dirBE  ])[kbe  ] =   c1over54* (drho+ (one+drho) * (three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
         (D.f[dirTW  ])[ktw  ] =   c1over54* (drho+ (one+drho) * (three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
         (D.f[dirTN  ])[ktn  ] =   c1over54* (drho+ (one+drho) * (three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
         (D.f[dirBS  ])[kbs  ] =   c1over54* (drho+ (one+drho) * (three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
         (D.f[dirBN  ])[kbn  ] =   c1over54* (drho+ (one+drho) * (three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
         (D.f[dirTS  ])[kts  ] =   c1over54* (drho+ (one+drho) * (three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
         (D.f[dirTNE ])[ktne ] =   c1over216*(drho+ (one+drho) * (three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
         (D.f[dirBSW ])[kbsw ] =   c1over216*(drho+ (one+drho) * (three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
         (D.f[dirBNE ])[kbne ] =   c1over216*(drho+ (one+drho) * (three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
         (D.f[dirTSW ])[ktsw ] =   c1over216*(drho+ (one+drho) * (three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
         (D.f[dirTSE ])[ktse ] =   c1over216*(drho+ (one+drho) * (three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
         (D.f[dirBNW ])[kbnw ] =   c1over216*(drho+ (one+drho) * (three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
         (D.f[dirBSE ])[kbse ] =   c1over216*(drho+ (one+drho) * (three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
         (D.f[dirTNW ])[ktnw ] =   c1over216*(drho+ (one+drho) * (three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
      }
   }
}















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitThS7( unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int* geoD,
                                       doubflo* Conc,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       unsigned int size_Mat,
                                       doubflo* DD7,
                                       bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions7 D7;
         if (EvenOrOdd==true)
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
         //////////////////////////////////////////////////////////////////////////
         doubflo ConcD = Conc[k];
         doubflo   vx1 = ux[k];
         doubflo   vx2 = uy[k];
         doubflo   vx3 = uz[k];
         doubflo lambdaD     = -three + sqrt(three);
         doubflo Diffusivity = c1o20;
         doubflo Lam         = -(c1o2+one/lambdaD);
         doubflo nue_d       = Lam/three;
         doubflo ae          = Diffusivity/nue_d - one;
         doubflo ux_sq       = vx1 * vx1;
         doubflo uy_sq       = vx2 * vx2;
         doubflo uz_sq       = vx3 * vx3;
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
         unsigned int kzero= k;
         unsigned int ke   = k;
         unsigned int kw   = neighborX[k];
         unsigned int kn   = k;
         unsigned int ks   = neighborY[k];
         unsigned int kt   = k;
         unsigned int kb   = neighborZ[k];
         //////////////////////////////////////////////////////////////////////////

         (D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
         (D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
         (D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
         (D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
         (D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
         (D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
         (D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
      }
   }
}












////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitThS27(unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int* geoD,
                                       doubflo* Conc,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       unsigned int size_Mat,
                                       doubflo* DD27,
                                       bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions27 D27;
         if (EvenOrOdd==true)
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
         //////////////////////////////////////////////////////////////////////////
         doubflo ConcD = Conc[k];
         doubflo   vx1 = ux[k];
         doubflo   vx2 = uy[k];
         doubflo   vx3 = uz[k];
         //doubflo lambdaD     = -three + sqrt(three);
         //doubflo Diffusivity = c1o20;
         //doubflo Lam         = -(c1o2+one/lambdaD);
         //doubflo nue_d       = Lam/three;
         //doubflo ae          = Diffusivity/nue_d - one;
         //doubflo ux_sq       = vx1 * vx1;
         //doubflo uy_sq       = vx2 * vx2;
         //doubflo uz_sq       = vx3 * vx3;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //D3Q7
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //index
         //unsigned int kzero= k;
         //unsigned int ke   = k;
         //unsigned int kw   = neighborX[k];
         //unsigned int kn   = k;
         //unsigned int ks   = neighborY[k];
         //unsigned int kt   = k;
         //unsigned int kb   = neighborZ[k];
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //(D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
         //(D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
         //(D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
         //(D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
         //(D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
         //(D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
         //(D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //D3Q27
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D27.f[dirZERO])[kzero] =   c8over27* ConcD*(one-cu_sq);
         (D27.f[dirE   ])[ke   ] =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
         (D27.f[dirW   ])[kw   ] =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
         (D27.f[dirN   ])[kn   ] =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
         (D27.f[dirS   ])[ks   ] =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D27.f[dirT   ])[kt   ] =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
         (D27.f[dirB   ])[kb   ] =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
         (D27.f[dirNE  ])[kne  ] =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D27.f[dirSW  ])[ksw  ] =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D27.f[dirSE  ])[kse  ] =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D27.f[dirNW  ])[knw  ] =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D27.f[dirTE  ])[kte  ] =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D27.f[dirBW  ])[kbw  ] =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D27.f[dirBE  ])[kbe  ] =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D27.f[dirTW  ])[ktw  ] =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D27.f[dirTN  ])[ktn  ] =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D27.f[dirBS  ])[kbs  ] =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D27.f[dirBN  ])[kbn  ] =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D27.f[dirTS  ])[kts  ] =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D27.f[dirTNE ])[ktne ] =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D27.f[dirBSW ])[kbsw ] =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D27.f[dirBNE ])[kbne ] =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D27.f[dirTSW ])[ktsw ] =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D27.f[dirTSE ])[ktse ] =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D27.f[dirBNW ])[kbnw ] =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D27.f[dirBSE ])[kbse ] =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D27.f[dirTNW ])[ktnw ] =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
   }
}















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitIncompAD7(unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int* geoD,
											doubflo* Conc,
											doubflo* ux,
											doubflo* uy,
											doubflo* uz,
											unsigned int size_Mat,
											doubflo* DD7,
											bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions7 D7;
         if (EvenOrOdd==true)
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
         //////////////////////////////////////////////////////////////////////////
         doubflo ConcD = Conc[k];
         doubflo   vx1 = ux[k];
         doubflo   vx2 = uy[k];
         doubflo   vx3 = uz[k];
         doubflo lambdaD     = -three + sqrt(three);
         doubflo Diffusivity = c1o20;
         doubflo Lam         = -(c1o2+one/lambdaD);
         doubflo nue_d       = Lam/three;
         doubflo ae          = Diffusivity/nue_d - one;
         doubflo ux_sq       = vx1 * vx1;
         doubflo uy_sq       = vx2 * vx2;
         doubflo uz_sq       = vx3 * vx3;
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
         unsigned int kzero= k;
         unsigned int ke   = k;
         unsigned int kw   = neighborX[k];
         unsigned int kn   = k;
         unsigned int ks   = neighborY[k];
         unsigned int kt   = k;
         unsigned int kb   = neighborZ[k];
         //////////////////////////////////////////////////////////////////////////

         (D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
         (D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
         (D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
         (D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
         (D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
         (D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
         (D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
      }
   }
}












////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitIncompAD27(unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int* geoD,
											doubflo* Conc,
											doubflo* ux,
											doubflo* uy,
											doubflo* uz,
											unsigned int size_Mat,
											doubflo* DD27,
											bool EvenOrOdd)
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
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions27 D27;
         if (EvenOrOdd==true)
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
         //////////////////////////////////////////////////////////////////////////
         doubflo ConcD = Conc[k];
         doubflo   vx1 = ux[k];
         doubflo   vx2 = uy[k];
         doubflo   vx3 = uz[k];
         //doubflo lambdaD     = -three + sqrt(three);
         //doubflo Diffusivity = c1o20;
         //doubflo Lam         = -(c1o2+one/lambdaD);
         //doubflo nue_d       = Lam/three;
         //doubflo ae          = Diffusivity/nue_d - one;
         //doubflo ux_sq       = vx1 * vx1;
         //doubflo uy_sq       = vx2 * vx2;
         //doubflo uz_sq       = vx3 * vx3;
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //D3Q7
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //index
         //unsigned int kzero= k;
         //unsigned int ke   = k;
         //unsigned int kw   = neighborX[k];
         //unsigned int kn   = k;
         //unsigned int ks   = neighborY[k];
         //unsigned int kt   = k;
         //unsigned int kb   = neighborZ[k];
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //(D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
         //(D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
         //(D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
         //(D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
         //(D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
         //(D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
         //(D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //D3Q27
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D27.f[dirZERO])[kzero] =   c8over27* ConcD*(one-cu_sq);
         (D27.f[dirE   ])[ke   ] =   c2over27* ConcD*(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
         (D27.f[dirW   ])[kw   ] =   c2over27* ConcD*(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
         (D27.f[dirN   ])[kn   ] =   c2over27* ConcD*(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
         (D27.f[dirS   ])[ks   ] =   c2over27* ConcD*(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D27.f[dirT   ])[kt   ] =   c2over27* ConcD*(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
         (D27.f[dirB   ])[kb   ] =   c2over27* ConcD*(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
         (D27.f[dirNE  ])[kne  ] =   c1over54* ConcD*(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D27.f[dirSW  ])[ksw  ] =   c1over54* ConcD*(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D27.f[dirSE  ])[kse  ] =   c1over54* ConcD*(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D27.f[dirNW  ])[knw  ] =   c1over54* ConcD*(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D27.f[dirTE  ])[kte  ] =   c1over54* ConcD*(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D27.f[dirBW  ])[kbw  ] =   c1over54* ConcD*(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D27.f[dirBE  ])[kbe  ] =   c1over54* ConcD*(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D27.f[dirTW  ])[ktw  ] =   c1over54* ConcD*(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D27.f[dirTN  ])[ktn  ] =   c1over54* ConcD*(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D27.f[dirBS  ])[kbs  ] =   c1over54* ConcD*(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D27.f[dirBN  ])[kbn  ] =   c1over54* ConcD*(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D27.f[dirTS  ])[kts  ] =   c1over54* ConcD*(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D27.f[dirTNE ])[ktne ] =   c1over216*ConcD*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D27.f[dirBSW ])[kbsw ] =   c1over216*ConcD*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D27.f[dirBNE ])[kbne ] =   c1over216*ConcD*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D27.f[dirTSW ])[ktsw ] =   c1over216*ConcD*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D27.f[dirTSE ])[ktse ] =   c1over216*ConcD*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D27.f[dirBNW ])[kbnw ] =   c1over216*ConcD*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D27.f[dirBSE ])[kbse ] =   c1over216*ConcD*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D27.f[dirTNW ])[ktnw ] =   c1over216*ConcD*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
   }
}



//test