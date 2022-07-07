/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInit27( int myid,
                                     int numprocs,
                                     real u0,
                                     unsigned int* geoD,
                                     unsigned int* neighborX,
                                     unsigned int* neighborY,
                                     unsigned int* neighborZ,
                                     real* vParabel,
                                     unsigned int size_Mat,
                                     unsigned int grid_nx, 
                                     unsigned int grid_ny, 
                                     unsigned int grid_nz, 
                                     real* DD,
                                     int lev,
                                     int maxlev)
{
   Distributions27 D;
   D.f[E   ] = &DD[E   *size_Mat];
   D.f[W   ] = &DD[W   *size_Mat];
   D.f[N   ] = &DD[N   *size_Mat];
   D.f[S   ] = &DD[S   *size_Mat];
   D.f[T   ] = &DD[T   *size_Mat];
   D.f[B   ] = &DD[B   *size_Mat];
   D.f[NE  ] = &DD[NE  *size_Mat];
   D.f[SW  ] = &DD[SW  *size_Mat];
   D.f[SE  ] = &DD[SE  *size_Mat];
   D.f[NW  ] = &DD[NW  *size_Mat];
   D.f[TE  ] = &DD[TE  *size_Mat];
   D.f[BW  ] = &DD[BW  *size_Mat];
   D.f[BE  ] = &DD[BE  *size_Mat];
   D.f[TW  ] = &DD[TW  *size_Mat];
   D.f[TN  ] = &DD[TN  *size_Mat];
   D.f[BS  ] = &DD[BS  *size_Mat];
   D.f[BN  ] = &DD[BN  *size_Mat];
   D.f[TS  ] = &DD[TS  *size_Mat];
   D.f[REST] = &DD[REST*size_Mat];
   D.f[TNE ] = &DD[TNE *size_Mat];
   D.f[TSW ] = &DD[TSW *size_Mat];
   D.f[TSE ] = &DD[TSE *size_Mat];
   D.f[TNW ] = &DD[TNW *size_Mat];
   D.f[BNE ] = &DD[BNE *size_Mat];
   D.f[BSW ] = &DD[BSW *size_Mat];
   D.f[BSE ] = &DD[BSE *size_Mat];
   D.f[BNW ] = &DD[BNW *size_Mat];
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
      real        radius  = grid_ny / 2.56;

      unsigned int distSq = (centerX-tx)*(centerX-tx)+(centerY-by)*(centerY-by)+(centerZ-bz)*(centerZ-bz);
      real radiSq = radius*radius;

      if( distSq < radiSq)        geoD[k] = GEO_SOLID;
   }
   //////////////////////////////////////////////////////////////////////////
   real drho = c0o1;
   real  vx1 = c0o1;
   real  vx2 = c0o1;
   real  vx3 = u0;
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

   real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

   (D.f[REST])[kzero] =   c8o27* (drho-cu_sq);
   (D.f[E   ])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
   (D.f[W   ])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
   (D.f[N   ])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
   (D.f[S   ])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
   (D.f[T   ])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
   (D.f[B   ])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
   (D.f[NE  ])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
   (D.f[SW  ])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
   (D.f[SE  ])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
   (D.f[NW  ])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
   (D.f[TE  ])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
   (D.f[BW  ])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
   (D.f[BE  ])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
   (D.f[TW  ])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
   (D.f[TN  ])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
   (D.f[BS  ])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
   (D.f[BN  ])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
   (D.f[TS  ])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
   (D.f[TNE ])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
   (D.f[BSW ])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
   (D.f[BNE ])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
   (D.f[TSW ])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
   (D.f[TSE ])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
   (D.f[BNW ])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
   (D.f[BSE ])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
   (D.f[TNW ])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

}
////////////////////////////////////////////////////////////////////////////////










////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBInitNonEqPartSP27( unsigned int* neighborX,
                                                unsigned int* neighborY,
                                                unsigned int* neighborZ,
                                                unsigned int* neighborWSB,
                                                unsigned int* geoD,
                                                real* rho,
                                                real* ux,
                                                real* uy,
                                                real* uz,
                                                unsigned int size_Mat,
                                                real* DD,
                                                real omega,
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
        BC = geoD[k];

        if( BC != GEO_SOLID &&  BC != GEO_VOID)
        {
            Distributions27 D;
            if (EvenOrOdd==true)
            {
                D.f[E   ] = &DD[E   *size_Mat];
                D.f[W   ] = &DD[W   *size_Mat];
                D.f[N   ] = &DD[N   *size_Mat];
                D.f[S   ] = &DD[S   *size_Mat];
                D.f[T   ] = &DD[T   *size_Mat];
                D.f[B   ] = &DD[B   *size_Mat];
                D.f[NE  ] = &DD[NE  *size_Mat];
                D.f[SW  ] = &DD[SW  *size_Mat];
                D.f[SE  ] = &DD[SE  *size_Mat];
                D.f[NW  ] = &DD[NW  *size_Mat];
                D.f[TE  ] = &DD[TE  *size_Mat];
                D.f[BW  ] = &DD[BW  *size_Mat];
                D.f[BE  ] = &DD[BE  *size_Mat];
                D.f[TW  ] = &DD[TW  *size_Mat];
                D.f[TN  ] = &DD[TN  *size_Mat];
                D.f[BS  ] = &DD[BS  *size_Mat];
                D.f[BN  ] = &DD[BN  *size_Mat];
                D.f[TS  ] = &DD[TS  *size_Mat];
                D.f[REST] = &DD[REST*size_Mat];
                D.f[TNE ] = &DD[TNE *size_Mat];
                D.f[TSW ] = &DD[TSW *size_Mat];
                D.f[TSE ] = &DD[TSE *size_Mat];
                D.f[TNW ] = &DD[TNW *size_Mat];
                D.f[BNE ] = &DD[BNE *size_Mat];
                D.f[BSW ] = &DD[BSW *size_Mat];
                D.f[BSE ] = &DD[BSE *size_Mat];
                D.f[BNW ] = &DD[BNW *size_Mat];
            }
            else
            {
                D.f[W   ] = &DD[E   *size_Mat];
                D.f[E   ] = &DD[W   *size_Mat];
                D.f[S   ] = &DD[N   *size_Mat];
                D.f[N   ] = &DD[S   *size_Mat];
                D.f[B   ] = &DD[T   *size_Mat];
                D.f[T   ] = &DD[B   *size_Mat];
                D.f[SW  ] = &DD[NE  *size_Mat];
                D.f[NE  ] = &DD[SW  *size_Mat];
                D.f[NW  ] = &DD[SE  *size_Mat];
                D.f[SE  ] = &DD[NW  *size_Mat];
                D.f[BW  ] = &DD[TE  *size_Mat];
                D.f[TE  ] = &DD[BW  *size_Mat];
                D.f[TW  ] = &DD[BE  *size_Mat];
                D.f[BE  ] = &DD[TW  *size_Mat];
                D.f[BS  ] = &DD[TN  *size_Mat];
                D.f[TN  ] = &DD[BS  *size_Mat];
                D.f[TS  ] = &DD[BN  *size_Mat];
                D.f[BN  ] = &DD[TS  *size_Mat];
                D.f[REST] = &DD[REST*size_Mat];
                D.f[BSW ] = &DD[TNE *size_Mat];
                D.f[BNE ] = &DD[TSW *size_Mat];
                D.f[BNW ] = &DD[TSE *size_Mat];
                D.f[BSE ] = &DD[TNW *size_Mat];
                D.f[TSW ] = &DD[BNE *size_Mat];
                D.f[TNE ] = &DD[BSW *size_Mat];
                D.f[TNW ] = &DD[BSE *size_Mat];
                D.f[TSE ] = &DD[BNW *size_Mat];
            }
            //////////////////////////////////////////////////////////////////////////
            real drho = rho[k];//0.0f;//
            real  vx1 = ux[k]; //0.0f;//
            real  vx2 = uy[k]; //0.0f;//
            real  vx3 = uz[k]; //0.0f;//
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
	        //////////////////////////////////////////////////////////////////////////////
	        //neighbor index
	        uint kPx   = neighborX[k];
	        uint kPy   = neighborY[k];
	        uint kPz   = neighborZ[k];
	        uint kMxyz = neighborWSB[k];
	        uint kMx   = neighborZ[neighborY[kMxyz]];
	        uint kMy   = neighborZ[neighborX[kMxyz]];
	        uint kMz   = neighborY[neighborX[kMxyz]];
            //////////////////////////////////////////////////////////////////////////
	        //getVeloX//
	        real vx1NeighborPx = ux[kPx];
	        real vx1NeighborMx = ux[kMx];
	        real vx1NeighborPy = ux[kPy];
	        real vx1NeighborMy = ux[kMy];
	        real vx1NeighborPz = ux[kPz];
	        real vx1NeighborMz = ux[kMz];
	        //getVeloY//
	        real vx2NeighborPx = uy[kPx];
	        real vx2NeighborMx = uy[kMx];
	        real vx2NeighborPy = uy[kPy];
	        real vx2NeighborMy = uy[kMy];
	        real vx2NeighborPz = uy[kPz];
	        real vx2NeighborMz = uy[kMz];
	        //getVeloZ//
	        real vx3NeighborPx = uz[kPx];
	        real vx3NeighborMx = uz[kMx];
	        real vx3NeighborPy = uz[kPy];
	        real vx3NeighborMy = uz[kMy];
	        real vx3NeighborPz = uz[kPz];
	        real vx3NeighborMz = uz[kMz];
            //////////////////////////////////////////////////////////////////////////

	        real dvx1dx = (vx1NeighborPx - vx1NeighborMx) / c2o1;
	        real dvx1dy = (vx1NeighborPy - vx1NeighborMy) / c2o1;
	        real dvx1dz = (vx1NeighborPz - vx1NeighborMz) / c2o1;

	        real dvx2dx = (vx2NeighborPx - vx2NeighborMx) / c2o1;
	        real dvx2dy = (vx2NeighborPy - vx2NeighborMy) / c2o1;
	        real dvx2dz = (vx2NeighborPz - vx2NeighborMz) / c2o1;

	        real dvx3dx = (vx3NeighborPx - vx3NeighborMx) / c2o1;
	        real dvx3dy = (vx3NeighborPy - vx3NeighborMy) / c2o1;
	        real dvx3dz = (vx3NeighborPz - vx3NeighborMz) / c2o1;

            //////////////////////////////////////////////////////////////////////////

            // the following code is copy and pasted from VirtualFluidsCore/Visitors/InitDistributionsBlockVisitor.cpp
            // i.e. Konstantins code

            real ax = dvx1dx;
            real ay = dvx1dy;
            real az = dvx1dz;

            real bx = dvx2dx;
            real by = dvx2dy;
            real bz = dvx2dz;

            real cx = dvx3dx;
            real cy = dvx3dy;
            real cz = dvx3dz;

            real eps_new = c1o1;
            real op      = c1o1;
            real o       = omega;

            real f_E    = eps_new *((5.*ax*omega + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op)/(54.*o*op));

            real f_N    =    f_E   + eps_new *((2.*(ax - by))/(9.*o));
            real f_T    =    f_E   + eps_new *((2.*(ax - cz))/(9.*o));
            real f_NE   =            eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op))/(54.*o*op));
            real f_SE   =    f_NE  + eps_new *((  ay + bx )/(9.*o));
            real f_TE   =            eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op))/(54.*o*op));
            real f_BE   =    f_TE  + eps_new *((  az + cx )/(9.*o));
            real f_TN   =            eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op)/(54.*o*op));
            real f_BN   =    f_TN  + eps_new *((  bz + cy )/(9.*o));
            real f_ZERO =            eps_new *((5.*(ax + by + cz))/(9.*op));
            real f_TNE  =            eps_new *(-(ay + az + bx + bz + cx + cy)/(72.*o));
            real f_TSW  =  - f_TNE - eps_new *((ay + bx)/(36.*o));
            real f_TSE  =  - f_TNE - eps_new *((az + cx)/(36.*o));
            real f_TNW  =  - f_TNE - eps_new *((bz + cy)/(36.*o));

            //////////////////////////////////////////////////////////////////////////
            real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
            
            (D.f[REST])[kzero] =   c8o27* (drho-cu_sq);
            (D.f[E   ])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
            (D.f[W   ])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
            (D.f[N   ])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
            (D.f[S   ])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
            (D.f[T   ])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
            (D.f[B   ])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
            (D.f[NE  ])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
            (D.f[SW  ])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
            (D.f[SE  ])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
            (D.f[NW  ])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
            (D.f[TE  ])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
            (D.f[BW  ])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
            (D.f[BE  ])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
            (D.f[TW  ])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
            (D.f[TN  ])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
            (D.f[BS  ])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
            (D.f[BN  ])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
            (D.f[TS  ])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
            (D.f[TNE ])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
            (D.f[BSW ])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
            (D.f[BNE ])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
            (D.f[TSW ])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
            (D.f[TSE ])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
            (D.f[BNW ])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
            (D.f[BSE ])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
            (D.f[TNW ])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

            //////////////////////////////////////////////////////////////////////////

            (D.f[REST])[kzero] += f_ZERO;
            (D.f[E   ])[ke   ] += f_E   ;
            (D.f[W   ])[kw   ] += f_E   ;
            (D.f[N   ])[kn   ] += f_N   ;
            (D.f[S   ])[ks   ] += f_N   ;
            (D.f[T   ])[kt   ] += f_T   ;
            (D.f[B   ])[kb   ] += f_T   ;
            (D.f[NE  ])[kne  ] += f_NE  ;
            (D.f[SW  ])[ksw  ] += f_NE  ;
            (D.f[SE  ])[kse  ] += f_SE  ;
            (D.f[NW  ])[knw  ] += f_SE  ;
            (D.f[TE  ])[kte  ] += f_TE  ;
            (D.f[BW  ])[kbw  ] += f_TE  ;
            (D.f[BE  ])[kbe  ] += f_BE  ;
            (D.f[TW  ])[ktw  ] += f_BE  ;
            (D.f[TN  ])[ktn  ] += f_TN  ;
            (D.f[BS  ])[kbs  ] += f_TN  ;
            (D.f[BN  ])[kbn  ] += f_BN  ;
            (D.f[TS  ])[kts  ] += f_BN  ;
            (D.f[TNE ])[ktne ] += f_TNE ;
            (D.f[BSW ])[kbsw ] += f_TNE ;
            (D.f[BNE ])[kbne ] += f_TSW ;
            (D.f[TSW ])[ktsw ] += f_TSW ;
            (D.f[TSE ])[ktse ] += f_TSE ;
            (D.f[BNW ])[kbnw ] += f_TSE ;
            (D.f[BSE ])[kbse ] += f_TNW ;
            (D.f[TNW ])[ktnw ] += f_TNW ;

            //////////////////////////////////////////////////////////////////////////
        }
	    else
	    {
		    //////////////////////////////////////////////////////////////////////////
		    Distributions27 D;
		    D.f[REST] = &DD[REST*size_Mat];
		    //////////////////////////////////////////////////////////////////////////
		    (D.f[REST])[k] = c96o1;
		    //////////////////////////////////////////////////////////////////////////
	    }
   }
}











































