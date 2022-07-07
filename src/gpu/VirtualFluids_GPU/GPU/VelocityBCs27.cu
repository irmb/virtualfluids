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
#include "lbm/constants/NumericConstants.h"
#include "KernelUtilities.h"

using namespace vf::lbm::constant;

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompPlusSlip27(
													real* vx,
													real* vy,
													real* vz,
													real* DD, 
													int* k_Q, 
													real* QQ,
													unsigned int numberOfBCnodes, 
													real om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
													bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if (k < numberOfBCnodes)
   {
	   ////////////////////////////////////////////////////////////////////////////////
	   real VeloX = vx[k];
	   real VeloY = vy[k];
	   real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
	   ////////////////////////////////////////////////////////////////////////////////
	   real *q_dirE, *q_dirW, *q_dirN, *q_dirS, *q_dirT, *q_dirB,
		   *q_dirNE, *q_dirSW, *q_dirSE, *q_dirNW, *q_dirTE, *q_dirBW,
		   *q_dirBE, *q_dirTW, *q_dirTN, *q_dirBS, *q_dirBN, *q_dirTS,
		   *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
		   *q_dirBSE, *q_dirBNW;
	   q_dirE = &QQ[E   * numberOfBCnodes];
	   q_dirW = &QQ[W   * numberOfBCnodes];
	   q_dirN = &QQ[N   * numberOfBCnodes];
	   q_dirS = &QQ[S   * numberOfBCnodes];
	   q_dirT = &QQ[T   * numberOfBCnodes];
	   q_dirB = &QQ[B   * numberOfBCnodes];
	   q_dirNE = &QQ[NE  * numberOfBCnodes];
	   q_dirSW = &QQ[SW  * numberOfBCnodes];
	   q_dirSE = &QQ[SE  * numberOfBCnodes];
	   q_dirNW = &QQ[NW  * numberOfBCnodes];
	   q_dirTE = &QQ[TE  * numberOfBCnodes];
	   q_dirBW = &QQ[BW  * numberOfBCnodes];
	   q_dirBE = &QQ[BE  * numberOfBCnodes];
	   q_dirTW = &QQ[TW  * numberOfBCnodes];
	   q_dirTN = &QQ[TN  * numberOfBCnodes];
	   q_dirBS = &QQ[BS  * numberOfBCnodes];
	   q_dirBN = &QQ[BN  * numberOfBCnodes];
	   q_dirTS = &QQ[TS  * numberOfBCnodes];
	   q_dirTNE = &QQ[TNE * numberOfBCnodes];
	   q_dirTSW = &QQ[TSW * numberOfBCnodes];
	   q_dirTSE = &QQ[TSE * numberOfBCnodes];
	   q_dirTNW = &QQ[TNW * numberOfBCnodes];
	   q_dirBNE = &QQ[BNE * numberOfBCnodes];
	   q_dirBSW = &QQ[BSW * numberOfBCnodes];
	   q_dirBSE = &QQ[BSE * numberOfBCnodes];
	   q_dirBNW = &QQ[BNW * numberOfBCnodes];
	   ////////////////////////////////////////////////////////////////////////////////
	   //index
	   unsigned int KQK = k_Q[k];
	   unsigned int kzero = KQK;
	   unsigned int ke = KQK;
	   unsigned int kw = neighborX[KQK];
	   unsigned int kn = KQK;
	   unsigned int ks = neighborY[KQK];
	   unsigned int kt = KQK;
	   unsigned int kb = neighborZ[KQK];
	   unsigned int ksw = neighborY[kw];
	   unsigned int kne = KQK;
	   unsigned int kse = ks;
	   unsigned int knw = kw;
	   unsigned int kbw = neighborZ[kw];
	   unsigned int kte = KQK;
	   unsigned int kbe = kb;
	   unsigned int ktw = kw;
	   unsigned int kbs = neighborZ[ks];
	   unsigned int ktn = KQK;
	   unsigned int kbn = kb;
	   unsigned int kts = ks;
	   unsigned int ktse = ks;
	   unsigned int kbnw = kbw;
	   unsigned int ktnw = kw;
	   unsigned int kbse = kbs;
	   unsigned int ktsw = ksw;
	   unsigned int kbne = kb;
	   unsigned int ktne = KQK;
	   unsigned int kbsw = neighborZ[ksw];
	   ////////////////////////////////////////////////////////////////////////////////
	   real f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE,
		   f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

	   f_W = (D.f[E])[ke];
	   f_E = (D.f[W])[kw];
	   f_S = (D.f[N])[kn];
	   f_N = (D.f[S])[ks];
	   f_B = (D.f[T])[kt];
	   f_T = (D.f[B])[kb];
	   f_SW = (D.f[NE])[kne];
	   f_NE = (D.f[SW])[ksw];
	   f_NW = (D.f[SE])[kse];
	   f_SE = (D.f[NW])[knw];
	   f_BW = (D.f[TE])[kte];
	   f_TE = (D.f[BW])[kbw];
	   f_TW = (D.f[BE])[kbe];
	   f_BE = (D.f[TW])[ktw];
	   f_BS = (D.f[TN])[ktn];
	   f_TN = (D.f[BS])[kbs];
	   f_TS = (D.f[BN])[kbn];
	   f_BN = (D.f[TS])[kts];
	   f_BSW = (D.f[TNE])[ktne];
	   f_BNE = (D.f[TSW])[ktsw];
	   f_BNW = (D.f[TSE])[ktse];
	   f_BSE = (D.f[TNW])[ktnw];
	   f_TSW = (D.f[BNE])[kbne];
	   f_TNE = (D.f[BSW])[kbsw];
	   f_TNW = (D.f[BSE])[kbse];
	   f_TSE = (D.f[BNW])[kbnw];
	   ////////////////////////////////////////////////////////////////////////////////
	   real vx1, vx2, vx3, drho, feq, q;
	   drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
		   f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
		   f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]);

	   vx1 = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		   ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
		   (f_E - f_W)) / (c1o1 + drho);


	   vx2 = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		   ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
		   (f_N - f_S)) / (c1o1 + drho);

	   vx3 = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		   (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
		   (f_T - f_B)) / (c1o1 + drho);

	   real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3) * (c1o1 + drho);

	   //////////////////////////////////////////////////////////////////////////
	   if (isEvenTimestep == false)
	   {
		   D.f[E] = &DD[E   *size_Mat];
		   D.f[W] = &DD[W   *size_Mat];
		   D.f[N] = &DD[N   *size_Mat];
		   D.f[S] = &DD[S   *size_Mat];
		   D.f[T] = &DD[T   *size_Mat];
		   D.f[B] = &DD[B   *size_Mat];
		   D.f[NE] = &DD[NE  *size_Mat];
		   D.f[SW] = &DD[SW  *size_Mat];
		   D.f[SE] = &DD[SE  *size_Mat];
		   D.f[NW] = &DD[NW  *size_Mat];
		   D.f[TE] = &DD[TE  *size_Mat];
		   D.f[BW] = &DD[BW  *size_Mat];
		   D.f[BE] = &DD[BE  *size_Mat];
		   D.f[TW] = &DD[TW  *size_Mat];
		   D.f[TN] = &DD[TN  *size_Mat];
		   D.f[BS] = &DD[BS  *size_Mat];
		   D.f[BN] = &DD[BN  *size_Mat];
		   D.f[TS] = &DD[TS  *size_Mat];
		   D.f[dirREST] = &DD[dirREST*size_Mat];
		   D.f[TNE] = &DD[TNE *size_Mat];
		   D.f[TSW] = &DD[TSW *size_Mat];
		   D.f[TSE] = &DD[TSE *size_Mat];
		   D.f[TNW] = &DD[TNW *size_Mat];
		   D.f[BNE] = &DD[BNE *size_Mat];
		   D.f[BSW] = &DD[BSW *size_Mat];
		   D.f[BSE] = &DD[BSE *size_Mat];
		   D.f[BNW] = &DD[BNW *size_Mat];
	   }
	   else
	   {
		   D.f[W] = &DD[E   *size_Mat];
		   D.f[E] = &DD[W   *size_Mat];
		   D.f[S] = &DD[N   *size_Mat];
		   D.f[N] = &DD[S   *size_Mat];
		   D.f[B] = &DD[T   *size_Mat];
		   D.f[T] = &DD[B   *size_Mat];
		   D.f[SW] = &DD[NE  *size_Mat];
		   D.f[NE] = &DD[SW  *size_Mat];
		   D.f[NW] = &DD[SE  *size_Mat];
		   D.f[SE] = &DD[NW  *size_Mat];
		   D.f[BW] = &DD[TE  *size_Mat];
		   D.f[TE] = &DD[BW  *size_Mat];
		   D.f[TW] = &DD[BE  *size_Mat];
		   D.f[BE] = &DD[TW  *size_Mat];
		   D.f[BS] = &DD[TN  *size_Mat];
		   D.f[TN] = &DD[BS  *size_Mat];
		   D.f[TS] = &DD[BN  *size_Mat];
		   D.f[BN] = &DD[TS  *size_Mat];
		   D.f[dirREST] = &DD[dirREST*size_Mat];
		   D.f[TNE] = &DD[BSW *size_Mat];
		   D.f[TSW] = &DD[BNE *size_Mat];
		   D.f[TSE] = &DD[BNW *size_Mat];
		   D.f[TNW] = &DD[BSE *size_Mat];
		   D.f[BNE] = &DD[TSW *size_Mat];
		   D.f[BSW] = &DD[TNE *size_Mat];
		   D.f[BSE] = &DD[TNW *size_Mat];
		   D.f[BNW] = &DD[TSE *size_Mat];
	   }
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   //Test
	   //(D.f[dirREST])[k]=c1o10;
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   //ToDo anders Klammern

	   /////To Slip Or Not To Slip?
	   // We assume slip BC if |vec(V_BC)|=1. To avoid problems we take V_BC*V_BC>0.99 (c99o100)
	   if (VeloX*VeloX + VeloY*VeloY + VeloZ*VeloZ > c99o100)
		{
		   // vt=v-(n \dot v) *n
		   // n=(VeloX,VeloY,VeloZ) a misuse of the velocity variable!
		   real normalV = VeloX*vx1 + VeloY*vx2 + VeloZ*vx3;
		   vx1 = vx1 - normalV*VeloX;
		   vx2 = vx2 - normalV*VeloY;
		   vx3 = vx3 - normalV*VeloZ;
		}
	  ////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);// - c2over27 * drho;
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);// - c1over54 * drho;
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);// - c1over54 * drho;
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);// - c1over216 * drho;
         //(D.f[BSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVeloDeviceEQ27(real* VeloX,
										   real* VeloY,
										   real* VeloZ,
                                           real* DD, 
                                           int* k_Q, 
                                           int numberOfBCnodes, 
                                           real om1, 
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned int size_Mat, 
                                           bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==true)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[E   ])[ke   ];
			real mfabb = (D.f[W   ])[kw   ];
			real mfbcb = (D.f[N   ])[kn   ];
			real mfbab = (D.f[S   ])[ks   ];
			real mfbbc = (D.f[T   ])[kt   ];
			real mfbba = (D.f[B   ])[kb   ];
			real mfccb = (D.f[NE  ])[kne  ];
			real mfaab = (D.f[SW  ])[ksw  ];
			real mfcab = (D.f[SE  ])[kse  ];
			real mfacb = (D.f[NW  ])[knw  ];
			real mfcbc = (D.f[TE  ])[kte  ];
			real mfaba = (D.f[BW  ])[kbw  ];
			real mfcba = (D.f[BE  ])[kbe  ];
			real mfabc = (D.f[TW  ])[ktw  ];
			real mfbcc = (D.f[TN  ])[ktn  ];
			real mfbaa = (D.f[BS  ])[kbs  ];
			real mfbca = (D.f[BN  ])[kbn  ];
			real mfbac = (D.f[TS  ])[kts  ];
			real mfbbb = (D.f[dirREST])[kzero];
			real mfccc = (D.f[TNE ])[ktne ];
			real mfaac = (D.f[TSW ])[ktsw ];
			real mfcac = (D.f[TSE ])[ktse ];
			real mfacc = (D.f[TNW ])[ktnw ];
			real mfcca = (D.f[BNE ])[kbne ];
			real mfaaa = (D.f[BSW ])[kbsw ];
			real mfcaa = (D.f[BSE ])[kbse ];
			real mfaca = (D.f[BNW ])[kbnw ];
			////////////////////////////////////////////////////////////////////////////////////
			real rho   = (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
							 mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
							 mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + c1o1);//!!!!Achtung + one
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			if(VeloX[k]!=c0o1) vvx = VeloX[k];
			if(VeloY[k]!=c0o1) vvy = VeloY[k];
			if(VeloZ[k]!=c0o1) vvz = VeloZ[k];
			////////////////////////////////////////////////////////////////////////////////////
			real vx2    = vvx * vvx;
			real vy2    = vvy * vvy;
			real vz2    = vvz * vvz;
			////////////////////////////////////////////////////////////////////////////////////
            real XXb    = -c2o3 + vx2;
            real XXc    = -c1o2 * (XXb + c1o1 + vvx);
            real XXa    = XXc + vvx;
            real YYb    = -c2o3 + vy2;
            real YYc    = -c1o2 * (YYb + c1o1 + vvy);
            real YYa    = YYc + vvy;
            real ZZb    = -c2o3 + vz2;
            real ZZc    = -c1o2 * (ZZb + c1o1 + vvz);
            real ZZa    = ZZc + vvz;
			////////////////////////////////////////////////////////////////////////////////////
            mfcbb = -rho * XXc * YYb * ZZb - c2o27 ; 
			mfabb = -rho * XXa * YYb * ZZb - c2o27 ;
			mfbcb = -rho * XXb * YYc * ZZb - c2o27 ;
			mfbab = -rho * XXb * YYa * ZZb - c2o27 ;
			mfbbc = -rho * XXb * YYb * ZZc - c2o27 ;
			mfbba = -rho * XXb * YYb * ZZa - c2o27 ;
			mfccb = -rho * XXc * YYc * ZZb - c1o54 ;
			mfaab = -rho * XXa * YYa * ZZb - c1o54 ;
			mfcab = -rho * XXc * YYa * ZZb - c1o54 ;
			mfacb = -rho * XXa * YYc * ZZb - c1o54 ;
			mfcbc = -rho * XXc * YYb * ZZc - c1o54 ;
			mfaba = -rho * XXa * YYb * ZZa - c1o54 ;
			mfcba = -rho * XXc * YYb * ZZa - c1o54 ;
			mfabc = -rho * XXa * YYb * ZZc - c1o54 ;
			mfbcc = -rho * XXb * YYc * ZZc - c1o54 ;
			mfbaa = -rho * XXb * YYa * ZZa - c1o54 ;
			mfbca = -rho * XXb * YYc * ZZa - c1o54 ;
			mfbac = -rho * XXb * YYa * ZZc - c1o54 ;
			mfbbb = -rho * XXb * YYb * ZZb - c8o27 ;
			mfccc = -rho * XXc * YYc * ZZc - c1o216;
			mfaac = -rho * XXa * YYa * ZZc - c1o216;
			mfcac = -rho * XXc * YYa * ZZc - c1o216;
			mfacc = -rho * XXa * YYc * ZZc - c1o216;
			mfcca = -rho * XXc * YYc * ZZa - c1o216;
			mfaaa = -rho * XXa * YYa * ZZa - c1o216;
			mfcaa = -rho * XXc * YYa * ZZa - c1o216;
			mfaca = -rho * XXa * YYc * ZZa - c1o216;
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			(D.f[E   ])[ke   ] = mfabb;//mfcbb;
			(D.f[W   ])[kw   ] = mfcbb;//mfabb;
			(D.f[N   ])[kn   ] = mfbab;//mfbcb;
			(D.f[S   ])[ks   ] = mfbcb;//mfbab;
			(D.f[T   ])[kt   ] = mfbba;//mfbbc;
			(D.f[B   ])[kb   ] = mfbbc;//mfbba;
			(D.f[NE  ])[kne  ] = mfaab;//mfccb;
			(D.f[SW  ])[ksw  ] = mfccb;//mfaab;
			(D.f[SE  ])[kse  ] = mfacb;//mfcab;
			(D.f[NW  ])[knw  ] = mfcab;//mfacb;
			(D.f[TE  ])[kte  ] = mfaba;//mfcbc;
			(D.f[BW  ])[kbw  ] = mfcbc;//mfaba;
			(D.f[BE  ])[kbe  ] = mfabc;//mfcba;
			(D.f[TW  ])[ktw  ] = mfcba;//mfabc;
			(D.f[TN  ])[ktn  ] = mfbaa;//mfbcc;
			(D.f[BS  ])[kbs  ] = mfbcc;//mfbaa;
			(D.f[BN  ])[kbn  ] = mfbac;//mfbca;
			(D.f[TS  ])[kts  ] = mfbca;//mfbac;
			(D.f[dirREST])[kzero] = mfbbb;//mfbbb;
			(D.f[TNE ])[ktne ] = mfaaa;//mfccc;
			(D.f[TSW ])[ktsw ] = mfcca;//mfaac;
			(D.f[TSE ])[ktse ] = mfaca;//mfcac;
			(D.f[TNW ])[ktnw ] = mfcaa;//mfacc;
			(D.f[BNE ])[kbne ] = mfaac;//mfcca;
			(D.f[BSW ])[kbsw ] = mfccc;//mfaaa;
			(D.f[BSE ])[kbse ] = mfacc;//mfcaa;
			(D.f[BNW ])[kbnw ] = mfcac;//mfaca;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVeloStreetDeviceEQ27(
	real* veloXfraction,
	real* veloYfraction,
	int*  naschVelo,
	real* DD,
	int*  naschIndex,
	int   numberOfStreetNodes,
	real  velocityRatio,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint  size_Mat,
	bool  isEvenTimestep)
{
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if (k < numberOfStreetNodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK   = naschIndex[k];
		unsigned int kzero = KQK;
		unsigned int ke    = KQK;
		unsigned int kw    = neighborX[KQK];
		unsigned int kn    = KQK;
		unsigned int ks    = neighborY[KQK];
		unsigned int kt    = KQK;
		unsigned int kb    = neighborZ[KQK];
		unsigned int ksw   = neighborY[kw];
		unsigned int kne   = KQK;
		unsigned int kse   = ks;
		unsigned int knw   = kw;
		unsigned int kbw   = neighborZ[kw];
		unsigned int kte   = KQK;
		unsigned int kbe   = kb;
		unsigned int ktw   = kw;
		unsigned int kbs   = neighborZ[ks];
		unsigned int ktn   = KQK;
		unsigned int kbn   = kb;
		unsigned int kts   = ks;
		unsigned int ktse  = ks;
		unsigned int kbnw  = kbw;
		unsigned int ktnw  = kw;
		unsigned int kbse  = kbs;
		unsigned int ktsw  = ksw;
		unsigned int kbne  = kb;
		unsigned int ktne  = KQK;
		unsigned int kbsw  = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////
		Distributions27 D;
		if (isEvenTimestep == true)
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
			D.f[TNE ] = &DD[BSW *size_Mat];
			D.f[TSW ] = &DD[BNE *size_Mat];
			D.f[TSE ] = &DD[BNW *size_Mat];
			D.f[TNW ] = &DD[BSE *size_Mat];
			D.f[BNE ] = &DD[TSW *size_Mat];
			D.f[BSW ] = &DD[TNE *size_Mat];
			D.f[BSE ] = &DD[TNW *size_Mat];
			D.f[BNW ] = &DD[TSE *size_Mat];
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// based on BGK Plus Comp
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		real mfcbb = (D.f[E   ])[ke   ];
		real mfabb = (D.f[W   ])[kw   ];
		real mfbcb = (D.f[N   ])[kn   ];
		real mfbab = (D.f[S   ])[ks   ];
		real mfbbc = (D.f[T   ])[kt   ];
		real mfbba = (D.f[B   ])[kb   ];
		real mfccb = (D.f[NE  ])[kne  ];
		real mfaab = (D.f[SW  ])[ksw  ];
		real mfcab = (D.f[SE  ])[kse  ];
		real mfacb = (D.f[NW  ])[knw  ];
		real mfcbc = (D.f[TE  ])[kte  ];
		real mfaba = (D.f[BW  ])[kbw  ];
		real mfcba = (D.f[BE  ])[kbe  ];
		real mfabc = (D.f[TW  ])[ktw  ];
		real mfbcc = (D.f[TN  ])[ktn  ];
		real mfbaa = (D.f[BS  ])[kbs  ];
		real mfbca = (D.f[BN  ])[kbn  ];
		real mfbac = (D.f[TS  ])[kts  ];
		real mfbbb = (D.f[dirREST])[kzero];
		real mfccc = (D.f[TNE ])[ktne ];
		real mfaac = (D.f[TSW ])[ktsw ];
		real mfcac = (D.f[TSE ])[ktse ];
		real mfacc = (D.f[TNW ])[ktnw ];
		real mfcca = (D.f[BNE ])[kbne ];
		real mfaaa = (D.f[BSW ])[kbsw ];
		real mfcaa = (D.f[BSE ])[kbse ];
		real mfaca = (D.f[BNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////////
		real rho = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
			        mfbac + mfbca + mfbaa + mfbcc + mfabc + mfcba + mfaba + mfcbc + mfacb + mfcab + mfaab + mfccb +
			        mfabb + mfcbb + mfbab + mfbcb + mfbba + mfbbc + mfbbb + c1o1);
		//!!!!Achtung + one
		////////////////////////////////////////////////////////////////////////////////////
		real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
			        (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
			          (mfcbb - mfabb)) / rho;
		real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
			        (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
			          (mfbcb - mfbab)) / rho;
		real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
			        (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
			          (mfbbc - mfbba)) / rho;
		////////////////////////////////////////////////////////////////////////////////////
		if (naschVelo[k] >= 0)
		{
			real VeloX = naschVelo[k] * veloXfraction[k] / velocityRatio;
			real VeloY = naschVelo[k] * veloYfraction[k] / velocityRatio;
			vvx = VeloX;
			vvy = VeloY;
		}
		////////////////////////////////////////////////////////////////////////////////////
		real vx2 = vvx * vvx;
		real vy2 = vvy * vvy;
		real vz2 = vvz * vvz;
		////////////////////////////////////////////////////////////////////////////////////
		real XXb = -c2o3 + vx2;
		real XXc = -c1o2 * (XXb + c1o1 + vvx);
		real XXa = XXc + vvx;
		real YYb = -c2o3 + vy2;
		real YYc = -c1o2 * (YYb + c1o1 + vvy);
		real YYa = YYc + vvy;
		real ZZb = -c2o3 + vz2;
		real ZZc = -c1o2 * (ZZb + c1o1 + vvz);
		real ZZa = ZZc + vvz;
		////////////////////////////////////////////////////////////////////////////////////
		mfcbb = -rho * XXc * YYb * ZZb - c2o27;
		mfabb = -rho * XXa * YYb * ZZb - c2o27;
		mfbcb = -rho * XXb * YYc * ZZb - c2o27;
		mfbab = -rho * XXb * YYa * ZZb - c2o27;
		mfbbc = -rho * XXb * YYb * ZZc - c2o27;
		mfbba = -rho * XXb * YYb * ZZa - c2o27;
		mfccb = -rho * XXc * YYc * ZZb - c1o54;
		mfaab = -rho * XXa * YYa * ZZb - c1o54;
		mfcab = -rho * XXc * YYa * ZZb - c1o54;
		mfacb = -rho * XXa * YYc * ZZb - c1o54;
		mfcbc = -rho * XXc * YYb * ZZc - c1o54;
		mfaba = -rho * XXa * YYb * ZZa - c1o54;
		mfcba = -rho * XXc * YYb * ZZa - c1o54;
		mfabc = -rho * XXa * YYb * ZZc - c1o54;
		mfbcc = -rho * XXb * YYc * ZZc - c1o54;
		mfbaa = -rho * XXb * YYa * ZZa - c1o54;
		mfbca = -rho * XXb * YYc * ZZa - c1o54;
		mfbac = -rho * XXb * YYa * ZZc - c1o54;
		mfbbb = -rho * XXb * YYb * ZZb - c8o27;
		mfccc = -rho * XXc * YYc * ZZc - c1o216;
		mfaac = -rho * XXa * YYa * ZZc - c1o216;
		mfcac = -rho * XXc * YYa * ZZc - c1o216;
		mfacc = -rho * XXa * YYc * ZZc - c1o216;
		mfcca = -rho * XXc * YYc * ZZa - c1o216;
		mfaaa = -rho * XXa * YYa * ZZa - c1o216;
		mfcaa = -rho * XXc * YYa * ZZa - c1o216;
		mfaca = -rho * XXa * YYc * ZZa - c1o216;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		(D.f[E   ])[ke   ] = mfabb;//mfcbb;
		(D.f[W   ])[kw   ] = mfcbb;//mfabb;
		(D.f[N   ])[kn   ] = mfbab;//mfbcb;
		(D.f[S   ])[ks   ] = mfbcb;//mfbab;
		(D.f[T   ])[kt   ] = mfbba;//mfbbc;
		(D.f[B   ])[kb   ] = mfbbc;//mfbba;
		(D.f[NE  ])[kne  ] = mfaab;//mfccb;
		(D.f[SW  ])[ksw  ] = mfccb;//mfaab;
		(D.f[SE  ])[kse  ] = mfacb;//mfcab;
		(D.f[NW  ])[knw  ] = mfcab;//mfacb;
		(D.f[TE  ])[kte  ] = mfaba;//mfcbc;
		(D.f[BW  ])[kbw  ] = mfcbc;//mfaba;
		(D.f[BE  ])[kbe  ] = mfabc;//mfcba;
		(D.f[TW  ])[ktw  ] = mfcba;//mfabc;
		(D.f[TN  ])[ktn  ] = mfbaa;//mfbcc;
		(D.f[BS  ])[kbs  ] = mfbcc;//mfbaa;
		(D.f[BN  ])[kbn  ] = mfbac;//mfbca;
		(D.f[TS  ])[kts  ] = mfbca;//mfbac;
		(D.f[dirREST])[kzero] = mfbbb;//mfbbb;
		(D.f[TNE ])[ktne ] = mfaaa;//mfccc;
		(D.f[TSW ])[ktsw ] = mfcca;//mfaac;
		(D.f[TSE ])[ktse ] = mfaca;//mfcac;
		(D.f[TNW ])[ktnw ] = mfcaa;//mfacc;
		(D.f[BNE ])[kbne ] = mfaac;//mfcca;
		(D.f[BSW ])[kbsw ] = mfccc;//mfaaa;
		(D.f[BSE ])[kbse ] = mfacc;//mfcaa;
		(D.f[BNW ])[kbnw ] = mfcac;//mfaca;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceIncompHighNu27(int inx,
													int iny,
													real* vx,
													real* vy,
													real* vz,
													real* DD, 
													int* k_Q, 
													real* QQ,
													unsigned int numberOfBCnodes, 
													real om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
													bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real VeloX = vx[k];
      real VeloY = vy[k];
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[E   ])[ke   ];
      f_W   = (D.f[W   ])[kw   ];
      f_N   = (D.f[N   ])[kn   ];
      f_S   = (D.f[S   ])[ks   ];
      f_T   = (D.f[T   ])[kt   ];
      f_B   = (D.f[B   ])[kb   ];
      f_NE  = (D.f[NE  ])[kne  ];
      f_SW  = (D.f[SW  ])[ksw  ];
      f_SE  = (D.f[SE  ])[kse  ];
      f_NW  = (D.f[NW  ])[knw  ];
      f_TE  = (D.f[TE  ])[kte  ];
      f_BW  = (D.f[BW  ])[kbw  ];
      f_BE  = (D.f[BE  ])[kbe  ];
      f_TW  = (D.f[TW  ])[ktw  ];
      f_TN  = (D.f[TN  ])[ktn  ];
      f_BS  = (D.f[BS  ])[kbs  ];
      f_BN  = (D.f[BN  ])[kbn  ];
      f_TS  = (D.f[TS  ])[kts  ];
      f_TNE = (D.f[TNE ])[ktne ];
      f_TSW = (D.f[TSW ])[ktsw ];
      f_TSE = (D.f[TSE ])[ktse ];
      f_TNW = (D.f[TNW ])[ktnw ];
      f_BNE = (D.f[BNE ])[kbne ];
      f_BSW = (D.f[BSW ])[kbsw ];
      f_BSE = (D.f[BSE ])[kbse ];
      f_BNW = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W));// / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S));// / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B));// / (one + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);// * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[W])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[E])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[S])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[N])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[B])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[T])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[SW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[NE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[NW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[SE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[TNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[BSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompHighNu27(
													real* vx,
													real* vy,
													real* vz,
													real* DD,
													int* k_Q,
													real* QQ,
													unsigned int numberOfBCnodes, 
													real om1,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
													bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real VeloX = vx[k];
      real VeloY = vy[k];
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[E   ])[ke   ];
      f_W   = (D.f[W   ])[kw   ];
      f_N   = (D.f[N   ])[kn   ];
      f_S   = (D.f[S   ])[ks   ];
      f_T   = (D.f[T   ])[kt   ];
      f_B   = (D.f[B   ])[kb   ];
      f_NE  = (D.f[NE  ])[kne  ];
      f_SW  = (D.f[SW  ])[ksw  ];
      f_SE  = (D.f[SE  ])[kse  ];
      f_NW  = (D.f[NW  ])[knw  ];
      f_TE  = (D.f[TE  ])[kte  ];
      f_BW  = (D.f[BW  ])[kbw  ];
      f_BE  = (D.f[BE  ])[kbe  ];
      f_TW  = (D.f[TW  ])[ktw  ];
      f_TN  = (D.f[TN  ])[ktn  ];
      f_BS  = (D.f[BS  ])[kbs  ];
      f_BN  = (D.f[BN  ])[kbn  ];
      f_TS  = (D.f[TS  ])[kts  ];
      f_TNE = (D.f[TNE ])[ktne ];
      f_TSW = (D.f[TSW ])[ktsw ];
      f_TSE = (D.f[TSE ])[ktse ];
      f_TNW = (D.f[TNW ])[ktnw ];
      f_BNE = (D.f[BNE ])[kbne ];
      f_BSW = (D.f[BSW ])[kbsw ];
      f_BSE = (D.f[BSE ])[kbse ];
      f_BNW = (D.f[BNW ])[kbnw ];
      //f_W    = (D.f[E   ])[ke   ];
      //f_E    = (D.f[W   ])[kw   ];
      //f_S    = (D.f[N   ])[kn   ];
      //f_N    = (D.f[S   ])[ks   ];
      //f_B    = (D.f[T   ])[kt   ];
      //f_T    = (D.f[B   ])[kb   ];
      //f_SW   = (D.f[NE  ])[kne  ];
      //f_NE   = (D.f[SW  ])[ksw  ];
      //f_NW   = (D.f[SE  ])[kse  ];
      //f_SE   = (D.f[NW  ])[knw  ];
      //f_BW   = (D.f[TE  ])[kte  ];
      //f_TE   = (D.f[BW  ])[kbw  ];
      //f_TW   = (D.f[BE  ])[kbe  ];
      //f_BE   = (D.f[TW  ])[ktw  ];
      //f_BS   = (D.f[TN  ])[ktn  ];
      //f_TN   = (D.f[BS  ])[kbs  ];
      //f_TS   = (D.f[BN  ])[kbn  ];
      //f_BN   = (D.f[TS  ])[kts  ];
      //f_BSW  = (D.f[TNE ])[ktne ];
      //f_BNE  = (D.f[TSW ])[ktsw ];
      //f_BNW  = (D.f[TSE ])[ktse ];
      //f_BSE  = (D.f[TNW ])[ktnw ];
      //f_TSW  = (D.f[BNE ])[kbne ];
      //f_TNE  = (D.f[BSW ])[kbsw ];
      //f_TNW  = (D.f[BSE ])[kbse ];
      //f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]=((c1o1 - q) * f_E + q * ((f_E + f_W) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloX     )) / (q + c1o1) ;
         //(D.f[W])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]=((c1o1 - q) * f_W + q * ((f_W + f_E) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloX     )) / (q + c1o1) ;
         //(D.f[E])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]=((c1o1 - q) * f_N + q * ((f_N + f_S) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloY     )) / (q + c1o1) ;
         //(D.f[S])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]=((c1o1 - q) * f_S + q * ((f_S + f_N) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloY     )) / (q + c1o1) ;
         //(D.f[N])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]=((c1o1 - q) * f_T + q * ((f_T + f_B) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*( VeloZ     )) / (q + c1o1) ;
         //(D.f[B])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]=((c1o1 - q) * f_B + q * ((f_B + f_T) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c2o27*(-VeloZ     )) / (q + c1o1) ;
         //(D.f[T])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]=((c1o1 - q) * f_NE + q * ((f_NE + f_SW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[SW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]=((c1o1 - q) * f_SW + q * ((f_SW + f_NE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[NE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]=((c1o1 - q) * f_SE + q * ((f_SE + f_NW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloY)) / (q + c1o1) ;
         //(D.f[NW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]=((c1o1 - q) * f_NW + q * ((f_NW + f_SE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloY)) / (q + c1o1) ;
         //(D.f[SE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]=((c1o1 - q) * f_TE + q * ((f_TE + f_BW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[BW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]=((c1o1 - q) * f_BW + q * ((f_BW + f_TE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[TE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]=((c1o1 - q) * f_BE + q * ((f_BE + f_TW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloX-VeloZ)) / (q + c1o1) ;
         //(D.f[TW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]=((c1o1 - q) * f_TW + q * ((f_TW + f_BE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloX+VeloZ)) / (q + c1o1) ;
         //(D.f[BE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]=((c1o1 - q) * f_TN + q * ((f_TN + f_BS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]=((c1o1 - q) * f_BS + q * ((f_BS + f_TN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]=((c1o1 - q) * f_BN + q * ((f_BN + f_TS) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*( VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]=((c1o1 - q) * f_TS + q * ((f_TS + f_BN) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o54*(-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]=((c1o1 - q) * f_TNE + q * ((f_TNE + f_BSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]=((c1o1 - q) * f_BSW + q * ((f_BSW + f_TNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]=((c1o1 - q) * f_BNE + q * ((f_BNE + f_TSW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]=((c1o1 - q) * f_TSW + q * ((f_TSW + f_BNE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]=((c1o1 - q) * f_TSE + q * ((f_TSE + f_BNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]=((c1o1 - q) * f_BNW + q * ((f_BNW + f_TSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]=((c1o1 - q) * f_BSE + q * ((f_BSE + f_TNW) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*( VeloX-VeloY-VeloZ)) / (q + c1o1) ;
         //(D.f[TNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]=((c1o1 - q) * f_TNW + q * ((f_TNW + f_BSE) * (c1o1 - om1) + om1 * c2o1 * feq) - c6o1*c1o216*(-VeloX+VeloY+VeloZ)) / (q + c1o1) ;
         //(D.f[BSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[BSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompZeroPress27(
														real* velocityX,
														real* velocityY,
														real* velocityZ,
														real* distribution, 
														int* subgridDistanceIndices, 
														real* subgridDistances,
														unsigned int numberOfBCnodes, 
														real omega, 
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														unsigned int numberOfLBnodes, 
														bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
	//! The velocity boundary condition is executed in the following steps
	//!
	////////////////////////////////////////////////////////////////////////////////
	//! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
	//!
   const unsigned  x = threadIdx.x;  // global x-index 
   const unsigned  y = blockIdx.x;   // global y-index 
   const unsigned  z = blockIdx.y;   // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(k < numberOfBCnodes)
   {

      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distribution, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[k];
      real VeloY = velocityY[k];
      real VeloZ = velocityZ[k];


      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
     
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int KQK  = subgridDistanceIndices[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                     f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                     f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[dirREST])[kzero]); 

      real vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                      ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                      (f_E - f_W)) / (c1o1 + drho); 
         

      real vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                       ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                       (f_N - f_S)) / (c1o1 + drho); 

      real vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                       (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                       (f_T - f_B)) / (c1o1 + drho); 
    
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distribution, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      real feq, q, velocityLB, velocityBC;
      q = (subgridD.q[E])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[W])[kw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_E, f_W, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[W])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[E])[ke] = getInterpolatedDistributionForVeloWithPressureBC(q, f_W, f_E, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[N])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[S])[ks] = getInterpolatedDistributionForVeloWithPressureBC(q, f_N, f_S, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[S])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[N])[kn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_S, f_N, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[T])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[B])[kb] = getInterpolatedDistributionForVeloWithPressureBC(q, f_T, f_B, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[B])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[T])[kt] = getInterpolatedDistributionForVeloWithPressureBC(q, f_B, f_T, feq, omega, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[NE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[SW])[ksw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NE, f_SW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[SW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[NE])[kne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SW, f_NE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[SE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[NW])[knw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SE, f_NW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[NW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[SE])[kse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NW, f_SE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[TE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[BW])[kbw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TE, f_BW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[BW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[TE])[kte] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BW, f_TE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[BE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[TW])[ktw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BE, f_TW, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[TW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[BE])[kbe] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TW, f_BE, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[TN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[BS])[kbs] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TN, f_BS, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[BS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[TN])[ktn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BS, f_TN, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[BN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[TS])[kts] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BN, f_TS, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[TS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[BN])[kbn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TS, f_BN, feq, omega, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[TNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[BSW])[kbsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNE, f_BSW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[BSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[TNE])[ktne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSW, f_TNE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[BNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[TSW])[ktsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNE, f_TSW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[TSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[BNE])[kbne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSW, f_BNE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[TSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[BNW])[kbnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSE, f_BNW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[BNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[TSE])[ktse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNW, f_TSE, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[BSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[TNW])[ktnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSE, f_TNW, feq, omega, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[TNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[BSE])[kbse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNW, f_BSE, feq, omega, drho, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompZeroPress1h27( int inx,
														int iny,
														real* vx,
														real* vy,
														real* vz,
														real* DD, 
														int* k_Q, 
														real* QQ,
														unsigned int numberOfBCnodes,
														real om1, 
														real Phi,
														real angularVelocity,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* coordX,
														real* coordY,
														real* coordZ,
														unsigned int size_Mat, 
														bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //real VeloX = vx[k];
      //real VeloY = vy[k];
      //real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////
		real VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		real VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		real coord0X = 281.125f;//7.5f;
		real coord0Y = 388.125f;//7.5f;
		real ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		real uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		real VeloXpur=VeloX;
		real VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real vx1, vx2, vx3, drho, feq, q, cu_sq;
	  ///////// equilibrium BC
	  cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
	  VeloXpur*=-c1o1;
	  VeloYpur*=-c1o1;
	  vx1=VeloX;
	  vx2=VeloY;
	  vx3=c0o1;
	  drho=c0o1;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*( VeloXpur        )+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[W])[kw]= feq - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(-VeloXpur        )+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[E])[ke]= feq - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(    VeloYpur     )+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[S])[ks]= feq - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(   -VeloYpur     )+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[N])[kn]= feq - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[B])[kb]= feq - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[T])[kt]= feq - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur+VeloYpur    )+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SW])[ksw]= feq - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur-VeloYpur    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NE])[kne]= feq - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur-VeloYpur    )+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[NW])[knw]= feq - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur+VeloYpur    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[SE])[kse]= feq - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BW])[kbw]= feq - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TE])[kte]= feq - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*( VeloXpur    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TW])[ktw]= feq - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(-VeloXpur    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BE])[kbe]= feq - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(     VeloYpur+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BS])[kbs]= feq - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(    -VeloYpur-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TN])[ktn]= feq - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(     VeloYpur-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TS])[kts]= feq - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c3o1*(    -VeloYpur+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BN])[kbn]= feq - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSW])[kbsw]= feq - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNE])[ktne]= feq - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSW])[ktsw]= feq - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNE])[kbne]= feq - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BNW])[kbnw]= feq - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TSE])[ktse]= feq - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[TNW])[ktnw]= feq - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[BSE])[kbse]= feq - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_BC_Vel_West_27( int nx, 
                                              int ny, 
                                              int nz, 
                                              int itz, 
                                              unsigned int* bcMatD, 
                                              unsigned int* neighborX,
                                              unsigned int* neighborY,
                                              unsigned int* neighborZ,
                                              real* DD, 
                                              unsigned int size_Mat, 
                                              bool isEvenTimestep, 
                                              real u0x, 
                                              unsigned int grid_nx, 
                                              unsigned int grid_ny, 
                                              real om) 
{
   //thread-index
   unsigned int ity = blockIdx.x;
   unsigned int itx = threadIdx.x;

   unsigned int  k, nxny;                   // Zugriff auf arrays im device

   unsigned int  x = itx + STARTOFFX;  // Globaler x-Index 
   unsigned int  y = ity + STARTOFFY;  // Globaler y-Index 
   unsigned int  z = itz + STARTOFFZ;  // Globaler z-Index 

   k = nx*(ny*z + y) + x;
   nxny = nx*ny;
   unsigned int k1 = k+nxny;

   if( bcMatD[k] == GEO_VELO )
   {
      Distributions27 D;
      if (isEvenTimestep==true)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }

      ////////////////////////////////////////////////////////////////////////////////
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
      //unsigned int kzero= k;
      //unsigned int ke   = k;
      //unsigned int kw   = k + 1;
      //unsigned int kn   = k;
      //unsigned int ks   = k + nx;
      //unsigned int kt   = k;
      //unsigned int kb   = k + nxny;
      //unsigned int ksw  = k + nx + 1;
      //unsigned int kne  = k;
      //unsigned int kse  = k + nx;
      //unsigned int knw  = k + 1;
      //unsigned int kbw  = k + nxny + 1;
      //unsigned int kte  = k;
      //unsigned int kbe  = k + nxny;
      //unsigned int ktw  = k + 1;
      //unsigned int kbs  = k + nxny + nx;
      //unsigned int ktn  = k;
      //unsigned int kbn  = k + nxny;
      //unsigned int kts  = k + nx;
      //unsigned int ktse = k + nx;
      //unsigned int kbnw = k + nxny + 1;
      //unsigned int ktnw = k + 1;
      //unsigned int kbse = k + nxny + nx;
      //unsigned int ktsw = k + nx + 1;
      //unsigned int kbne = k + nxny;
      //unsigned int ktne = k;
      //unsigned int kbsw = k + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int k1zero= k1;
      unsigned int k1e   = k1;
      unsigned int k1w   = neighborX[k1];
      unsigned int k1n   = k1;
      unsigned int k1s   = neighborY[k1];
      unsigned int k1t   = k1;
      unsigned int k1b   = neighborZ[k1];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = k1;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = k1;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = k1;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = k1;
      unsigned int k1bsw = neighborZ[k1sw];
      //unsigned int k1zero= k1;
      //unsigned int k1e   = k1;
      //unsigned int k1w   = k1 + 1;
      //unsigned int k1n   = k1;
      //unsigned int k1s   = k1 + nx;
      //unsigned int k1t   = k1;
      //unsigned int k1b   = k1 + nxny;
      //unsigned int k1sw  = k1 + nx + 1;
      //unsigned int k1ne  = k1;
      //unsigned int k1se  = k1 + nx;
      //unsigned int k1nw  = k1 + 1;
      //unsigned int k1bw  = k1 + nxny + 1;
      //unsigned int k1te  = k1;
      //unsigned int k1be  = k1 + nxny;
      //unsigned int k1tw  = k1 + 1;
      //unsigned int k1bs  = k1 + nxny + nx;
      //unsigned int k1tn  = k1;
      //unsigned int k1bn  = k1 + nxny;
      //unsigned int k1ts  = k1 + nx;
      //unsigned int k1tse = k1 + nx;
      //unsigned int k1bnw = k1 + nxny + 1;
      //unsigned int k1tnw = k1 + 1;
      //unsigned int k1bse = k1 + nxny + nx;
      //unsigned int k1tsw = k1 + nx + 1;
      //unsigned int k1bne = k1 + nxny;
      //unsigned int k1tne = k1;
      //unsigned int k1bsw = k1 + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[E   ])[k1e   ];
      f1_E    = (D.f[W   ])[k1w   ];
      f1_S    = (D.f[N   ])[k1n   ];
      f1_N    = (D.f[S   ])[k1s   ];
      f1_B    = (D.f[T   ])[k1t   ];
      f1_T    = (D.f[B   ])[k1b   ];
      f1_SW   = (D.f[NE  ])[k1ne  ];
      f1_NE   = (D.f[SW  ])[k1sw  ];
      f1_NW   = (D.f[SE  ])[k1se  ];
      f1_SE   = (D.f[NW  ])[k1nw  ];
      f1_BW   = (D.f[TE  ])[k1te  ];
      f1_TE   = (D.f[BW  ])[k1bw  ];
      f1_TW   = (D.f[BE  ])[k1be  ];
      f1_BE   = (D.f[TW  ])[k1tw  ];
      f1_BS   = (D.f[TN  ])[k1tn  ];
      f1_TN   = (D.f[BS  ])[k1bs  ];
      f1_TS   = (D.f[BN  ])[k1bn  ];
      f1_BN   = (D.f[TS  ])[k1ts  ];
      f1_ZERO = (D.f[dirREST])[k1zero];
      f1_BSW  = (D.f[TNE ])[k1tne ];
      f1_BNE  = (D.f[TSW ])[k1tsw ];
      f1_BNW  = (D.f[TSE ])[k1tse ];
      f1_BSE  = (D.f[TNW ])[k1tnw ];
      f1_TSW  = (D.f[BNE ])[k1bne ];
      f1_TNE  = (D.f[BSW ])[k1bsw ];
      f1_TNW  = (D.f[BSE ])[k1bse ];
      f1_TSE  = (D.f[BNW ])[k1bnw ];

      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real drho = drho1;
      real  vx1 = c0o1;
      real  vx2 = c0o1;
      real  vx3 = u0x;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D.f[dirREST])[kzero] =   c8o27* (drho-cu_sq);
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
   __syncthreads();
}          
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevPlainBB27(
   real* velocityX,
   real* velocityY,
   real* velocityZ,
   real* distributions,
   int* subgridDistanceIndices,
   real* subgridDistances,
   uint numberOfBCnodes,
   uint* neighborX,
   uint* neighborY,
   uint* neighborZ,
   uint numberOfLBnodes,
   bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;   // global x-index
   const unsigned  y = blockIdx.x;    // global y-index
   const unsigned  z = blockIdx.y;    // global z-index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   // run for all indices in size of boundary condition (numberOfBCnodes)
   if(k < numberOfBCnodes)
   {
       //////////////////////////////////////////////////////////////////////////
       //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
       //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
       //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[k];
      real VeloY = velocityY[k];
      real VeloZ = velocityZ[k];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      uint indexOfBCnode = subgridDistanceIndices[k];
      uint ke   = indexOfBCnode;
      uint kw   = neighborX[indexOfBCnode];
      uint kn   = indexOfBCnode;
      uint ks   = neighborY[indexOfBCnode];
      uint kt   = indexOfBCnode;
      uint kb   = neighborZ[indexOfBCnode];
      uint ksw  = neighborY[kw];
      uint kne  = indexOfBCnode;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = indexOfBCnode;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = indexOfBCnode;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = indexOfBCnode;
      uint kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - rewrite distributions if there is a sub-grid distance (q) in same direction
      real q;
      q = (subgridD.q[E  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[W  ])[kw  ]=f_E   + c4o9  * (-VeloX);
      q = (subgridD.q[W  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[E  ])[ke  ]=f_W   + c4o9  * ( VeloX);
      q = (subgridD.q[N  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[S  ])[ks  ]=f_N   + c4o9  * (-VeloY);
      q = (subgridD.q[S  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[N  ])[kn  ]=f_S   + c4o9  * ( VeloY);
      q = (subgridD.q[T  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[B  ])[kb  ]=f_T   + c4o9  * (-VeloZ);
      q = (subgridD.q[B  ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[T  ])[kt  ]=f_B   + c4o9  * ( VeloZ);
      q = (subgridD.q[NE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[SW ])[ksw ]=f_NE  + c1o9  * (-VeloX - VeloY);
      q = (subgridD.q[SW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[NE ])[kne ]=f_SW  + c1o9  * ( VeloX + VeloY);
      q = (subgridD.q[SE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[NW ])[knw ]=f_SE  + c1o9  * (-VeloX + VeloY);
      q = (subgridD.q[NW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[SE ])[kse ]=f_NW  + c1o9  * ( VeloX - VeloY);
      q = (subgridD.q[TE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BW ])[kbw ]=f_TE  + c1o9  * (-VeloX - VeloZ);
      q = (subgridD.q[BW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TE ])[kte ]=f_BW  + c1o9  * ( VeloX + VeloZ);
      q = (subgridD.q[BE ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TW ])[ktw ]=f_BE  + c1o9  * (-VeloX + VeloZ);
      q = (subgridD.q[TW ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BE ])[kbe ]=f_TW  + c1o9  * ( VeloX - VeloZ);
      q = (subgridD.q[TN ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BS ])[kbs ]=f_TN  + c1o9  * (-VeloY - VeloZ);
      q = (subgridD.q[BS ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TN ])[ktn ]=f_BS  + c1o9  * ( VeloY + VeloZ);
      q = (subgridD.q[BN ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TS ])[kts ]=f_BN  + c1o9  * (-VeloY + VeloZ);
      q = (subgridD.q[TS ])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BN ])[kbn ]=f_TS  + c1o9  * ( VeloY - VeloZ);
      q = (subgridD.q[TNE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BSW])[kbsw]=f_TNE + c1o36 * (-VeloX - VeloY - VeloZ);
      q = (subgridD.q[BSW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TNE])[ktne]=f_BSW + c1o36 * ( VeloX + VeloY + VeloZ);
      q = (subgridD.q[BNE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TSW])[ktsw]=f_BNE + c1o36 * (-VeloX - VeloY + VeloZ);
      q = (subgridD.q[TSW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BNE])[kbne]=f_TSW + c1o36 * ( VeloX + VeloY - VeloZ);
      q = (subgridD.q[TSE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BNW])[kbnw]=f_TSE + c1o36 * (-VeloX + VeloY - VeloZ);
      q = (subgridD.q[BNW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TSE])[ktse]=f_BNW + c1o36 * ( VeloX - VeloY + VeloZ);
      q = (subgridD.q[BSE])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[TNW])[ktnw]=f_BSE + c1o36 * (-VeloX + VeloY + VeloZ);
      q = (subgridD.q[TNW])[k];   if (q>=c0o1 && q<=c1o1)    (dist.f[BSE])[kbse]=f_TNW + c1o36 * ( VeloX - VeloY - VeloZ);
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevCouette27(real* vx,
											real* vy,
	 										real* vz,
											real* DD,
											int* k_Q, 
											real* QQ,
											unsigned int numberOfBCnodes, 
											real om1, 
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
											bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
	  ////////////////////////////////////////////////////////////////////////////////
	  real VeloX = vx[k];
	  real VeloY = vy[k];
	  real VeloZ = vz[k];
      ////////////////////////////////////////////////////////////////////////////////
      real*q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			 *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			 *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			 *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			 *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
     
      ////////////////////////////////////////////////////////////////////////////////
      real f_W    = (D.f[E   ])[ke   ];
      real f_E    = (D.f[W   ])[kw   ];
      real f_S    = (D.f[N   ])[kn   ];
      real f_N    = (D.f[S   ])[ks   ];
      real f_B    = (D.f[T   ])[kt   ];
      real f_T    = (D.f[B   ])[kb   ];
      real f_SW   = (D.f[NE  ])[kne  ];
      real f_NE   = (D.f[SW  ])[ksw  ];
      real f_NW   = (D.f[SE  ])[kse  ];
      real f_SE   = (D.f[NW  ])[knw  ];
      real f_BW   = (D.f[TE  ])[kte  ];
      real f_TE   = (D.f[BW  ])[kbw  ];
      real f_TW   = (D.f[BE  ])[kbe  ];
      real f_BE   = (D.f[TW  ])[ktw  ];
      real f_BS   = (D.f[TN  ])[ktn  ];
      real f_TN   = (D.f[BS  ])[kbs  ];
      real f_TS   = (D.f[BN  ])[kbn  ];
      real f_BN   = (D.f[TS  ])[kts  ];
      real f_BSW  = (D.f[TNE ])[ktne ];
      real f_BNE  = (D.f[TSW ])[ktsw ];
      real f_BNW  = (D.f[TSE ])[ktse ];
      real f_BSE  = (D.f[TNW ])[ktnw ];
      real f_TSW  = (D.f[BNE ])[kbne ];
      real f_TNE  = (D.f[BSW ])[kbsw ];
      real f_TNW  = (D.f[BSE ])[kbse ];
      real f_TSE  = (D.f[BNW ])[kbnw ];
	  ////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////               FlowDirection Y !!!!!!!!!!                                                           ///////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //calculate velocity
	  //real vx1 = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W);
	  real vx2 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S);
	  //real vx3 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B);
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //constant
	  real on=c0o1;//c1o2;//one;
	  real ms=-c6o1;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //2nd order moment
	  real kxxMyyFromfcNEQ = c0o1;//-c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));		//all E+W minus all N+S (no combinations of xy left)

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set distributions
      real q;
      q = q_dirE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[W  ])[kw  ]=f_E   + ms*c2o27  * VeloX;	
      q = q_dirW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[E  ])[ke  ]=f_W   - ms*c2o27  * VeloX;	
      q = q_dirN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[S  ])[ks  ]=f_N   + ms*c2o27  * VeloY;	
      q = q_dirS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[N  ])[kn  ]=f_S   - ms*c2o27  * VeloY;	
	  q = q_dirT[k];	if (q>=c0o1 && q<=c1o1)	(D.f[B  ])[kb  ]=f_T   + ms*c2o27  * VeloZ - c3o2*c2o27*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirB[k];	if (q>=c0o1 && q<=c1o1)	(D.f[T  ])[kt  ]=f_B   - ms*c2o27  * VeloZ;
      q = q_dirNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[SW ])[ksw ]=f_NE  + ms*c1o54  * VeloX + ms*c1o54  * VeloY;
	  q = q_dirSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[NE ])[kne ]=f_SW  - ms*c1o54  * VeloX - ms*c1o54  * VeloY;
	  q = q_dirSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[NW ])[knw ]=f_SE  + ms*c1o54  * VeloX - ms*c1o54  * VeloY;
	  q = q_dirNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[SE ])[kse ]=f_NW  - ms*c1o54  * VeloX + ms*c1o54  * VeloY;
	  q = q_dirTE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BW ])[kbw ]=f_TE  + ms*c1o54  * VeloX + ms*c1o54  * VeloZ - c3o2*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TE ])[kte ]=f_BW  - ms*c1o54  * VeloX - ms*c1o54  * VeloZ;
	  q = q_dirBE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TW ])[ktw ]=f_BE  + ms*c1o54  * VeloX - ms*c1o54  * VeloZ;
	  q = q_dirTW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BE ])[kbe ]=f_TW  - ms*c1o54  * VeloX + ms*c1o54  * VeloZ - c3o2*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirTN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BS ])[kbs ]=f_TN  + ms*c1o54  * VeloY + ms*c1o54  * VeloZ + c3o1*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TN ])[ktn ]=f_BS  - ms*c1o54  * VeloY - ms*c1o54  * VeloZ;
	  q = q_dirBN[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TS ])[kts ]=f_BN  + ms*c1o54  * VeloY - ms*c1o54  * VeloZ;
	  q = q_dirTS[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BN ])[kbn ]=f_TS  - ms*c1o54  * VeloY + ms*c1o54  * VeloZ + c3o1*c1o54*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
      q = q_dirTNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BSW])[kbsw]=f_TNE + ms*c1o216 * VeloX + ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TNE])[ktne]=f_BSW - ms*c1o216 * VeloX - ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirBNE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TSW])[ktsw]=f_BNE + ms*c1o216 * VeloX + ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirTSW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BNE])[kbne]=f_TSW - ms*c1o216 * VeloX - ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirTSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BNW])[kbnw]=f_TSE + ms*c1o216 * VeloX - ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TSE])[ktse]=f_BNW - ms*c1o216 * VeloX + ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirBSE[k];	if (q>=c0o1 && q<=c1o1)	(D.f[TNW])[ktnw]=f_BSE + ms*c1o216 * VeloX - ms*c1o216 * VeloY - ms*c1o216 * VeloZ;
      q = q_dirTNW[k];	if (q>=c0o1 && q<=c1o1)	(D.f[BSE])[kbse]=f_TNW - ms*c1o216 * VeloX + ms*c1o216 * VeloY + ms*c1o216 * VeloZ + c3o1*c1o216*((c2o1*VeloY-vx2)*(c2o1*VeloY-vx2)-vx2*vx2)*on;
      //q = q_dirE[k];	if (q>=zero && q<=one)	(D.f[W  ])[kw  ]=f_E   + ms*c2over27  * VeloX;	
   //   q = q_dirW[k];	if (q>=zero && q<=one)	(D.f[E  ])[ke  ]=f_W   - ms*c2over27  * VeloX;	
   //   q = q_dirN[k];	if (q>=zero && q<=one)	(D.f[S  ])[ks  ]=f_N   + ms*c2over27  * VeloY;	
   //   q = q_dirS[k];	if (q>=zero && q<=one)	(D.f[N  ])[kn  ]=f_S   - ms*c2over27  * VeloY;	
	  //q = q_dirT[k];	if (q>=zero && q<=one)	(D.f[B  ])[kb  ]=f_T   + ms*c2over27  * VeloZ - c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirB[k];	if (q>=zero && q<=one)	(D.f[T  ])[kt  ]=f_B   - ms*c2over27  * VeloZ;
   //   q = q_dirNE[k];	if (q>=zero && q<=one)	(D.f[SW ])[ksw ]=f_NE  + ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirSW[k];	if (q>=zero && q<=one)	(D.f[NE ])[kne ]=f_SW  - ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirSE[k];	if (q>=zero && q<=one)	(D.f[NW ])[knw ]=f_SE  + ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirNW[k];	if (q>=zero && q<=one)	(D.f[SE ])[kse ]=f_NW  - ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirTE[k];	if (q>=zero && q<=one)	(D.f[BW ])[kbw ]=f_TE  + ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBW[k];	if (q>=zero && q<=one)	(D.f[TE ])[kte ]=f_BW  - ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirBE[k];	if (q>=zero && q<=one)	(D.f[TW ])[ktw ]=f_BE  + ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirTW[k];	if (q>=zero && q<=one)	(D.f[BE ])[kbe ]=f_TW  - ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirTN[k];	if (q>=zero && q<=one)	(D.f[BS ])[kbs ]=f_TN  + ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBS[k];	if (q>=zero && q<=one)	(D.f[TN ])[ktn ]=f_BS  - ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirBN[k];	if (q>=zero && q<=one)	(D.f[TS ])[kts ]=f_BN  + ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirTS[k];	if (q>=zero && q<=one)	(D.f[BN ])[kbn ]=f_TS  - ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTNE[k];	if (q>=zero && q<=one)	(D.f[BSW])[kbsw]=f_TNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBSW[k];	if (q>=zero && q<=one)	(D.f[TNE])[ktne]=f_BSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBNE[k];	if (q>=zero && q<=one)	(D.f[TSW])[ktsw]=f_BNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTSW[k];	if (q>=zero && q<=one)	(D.f[BNE])[kbne]=f_TSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTSE[k];	if (q>=zero && q<=one)	(D.f[BNW])[kbnw]=f_TSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBNW[k];	if (q>=zero && q<=one)	(D.f[TSE])[ktse]=f_BNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBSE[k];	if (q>=zero && q<=one)	(D.f[TNW])[ktnw]=f_BSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTNW[k];	if (q>=zero && q<=one)	(D.f[BSE])[kbse]=f_TNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDev1h27( int inx,
										int iny,
										real* vx,
										real* vy,
										real* vz,
										real* DD, 
										int* k_Q, 
										real* QQ,
										unsigned int numberOfBCnodes, 
										real om1,
										real Phi,
										real angularVelocity,
										unsigned int* neighborX,
										unsigned int* neighborY,
										unsigned int* neighborZ,
										real* coordX,
										real* coordY,
										real* coordZ,
										unsigned int size_Mat, 
										bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
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
		D.f[dirREST] = &DD[dirREST*size_Mat];
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
		D.f[dirREST] = &DD[dirREST*size_Mat];
		D.f[TNE ] = &DD[BSW *size_Mat];
		D.f[TSW ] = &DD[BNE *size_Mat];
		D.f[TSE ] = &DD[BNW *size_Mat];
		D.f[TNW ] = &DD[BSE *size_Mat];
		D.f[BNE ] = &DD[TSW *size_Mat];
		D.f[BSW ] = &DD[TNE *size_Mat];
		D.f[BSE ] = &DD[TNW *size_Mat];
		D.f[BNW ] = &DD[TSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		real VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		real VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		real coord0X = 281.125f;//7.5f;
		real coord0Y = 388.125f;//7.5f;
		real ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		real uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		real VeloXpur=VeloX;
		real VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[E   * numberOfBCnodes];
		q_dirW   = &QQ[W   * numberOfBCnodes];
		q_dirN   = &QQ[N   * numberOfBCnodes];
		q_dirS   = &QQ[S   * numberOfBCnodes];
		q_dirT   = &QQ[T   * numberOfBCnodes];
		q_dirB   = &QQ[B   * numberOfBCnodes];
		q_dirNE  = &QQ[NE  * numberOfBCnodes];
		q_dirSW  = &QQ[SW  * numberOfBCnodes];
		q_dirSE  = &QQ[SE  * numberOfBCnodes];
		q_dirNW  = &QQ[NW  * numberOfBCnodes];
		q_dirTE  = &QQ[TE  * numberOfBCnodes];
		q_dirBW  = &QQ[BW  * numberOfBCnodes];
		q_dirBE  = &QQ[BE  * numberOfBCnodes];
		q_dirTW  = &QQ[TW  * numberOfBCnodes];
		q_dirTN  = &QQ[TN  * numberOfBCnodes];
		q_dirBS  = &QQ[BS  * numberOfBCnodes];
		q_dirBN  = &QQ[BN  * numberOfBCnodes];
		q_dirTS  = &QQ[TS  * numberOfBCnodes];
		q_dirTNE = &QQ[TNE * numberOfBCnodes];
		q_dirTSW = &QQ[TSW * numberOfBCnodes];
		q_dirTSE = &QQ[TSE * numberOfBCnodes];
		q_dirTNW = &QQ[TNW * numberOfBCnodes];
		q_dirBNE = &QQ[BNE * numberOfBCnodes];
		q_dirBSW = &QQ[BSW * numberOfBCnodes];
		q_dirBSE = &QQ[BSE * numberOfBCnodes];
		q_dirBNW = &QQ[BNW * numberOfBCnodes];
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = k_Q[k];
		//unsigned int kzero= KQK;
		unsigned int ke   = KQK;
		unsigned int kw   = neighborX[KQK];
		unsigned int kn   = KQK;
		unsigned int ks   = neighborY[KQK];
		unsigned int kt   = KQK;
		unsigned int kb   = neighborZ[KQK];
		unsigned int ksw  = neighborY[kw];
		unsigned int kne  = KQK;
		unsigned int kse  = ks;
		unsigned int knw  = kw;
		unsigned int kbw  = neighborZ[kw];
		unsigned int kte  = KQK;
		unsigned int kbe  = kb;
		unsigned int ktw  = kw;
		unsigned int kbs  = neighborZ[ks];
		unsigned int ktn  = KQK;
		unsigned int kbn  = kb;
		unsigned int kts  = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
		//unsigned int nxny = nx*ny;
		//unsigned int kzero= KQK;
		//unsigned int ke   = KQK;
		//unsigned int kw   = KQK + 1;
		//unsigned int kn   = KQK;
		//unsigned int ks   = KQK + nx;
		//unsigned int kt   = KQK;
		//unsigned int kb   = KQK + nxny;
		//unsigned int ksw  = KQK + nx + 1;
		//unsigned int kne  = KQK;
		//unsigned int kse  = KQK + nx;
		//unsigned int knw  = KQK + 1;
		//unsigned int kbw  = KQK + nxny + 1;
		//unsigned int kte  = KQK;
		//unsigned int kbe  = KQK + nxny;
		//unsigned int ktw  = KQK + 1;
		//unsigned int kbs  = KQK + nxny + nx;
		//unsigned int ktn  = KQK;
		//unsigned int kbn  = KQK + nxny;
		//unsigned int kts  = KQK + nx;
		//unsigned int ktse = KQK + nx;
		//unsigned int kbnw = KQK + nxny + 1;
		//unsigned int ktnw = KQK + 1;
		//unsigned int kbse = KQK + nxny + nx;
		//unsigned int ktsw = KQK + nx + 1;
		//unsigned int kbne = KQK + nxny;
		//unsigned int ktne = KQK;
		//unsigned int kbsw = KQK + nxny + nx + 1;
		////////////////////////////////////////////////////////////////////////////////
		//real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		//	f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		//f_W    = (D.f[E   ])[ke   ];
		//f_E    = (D.f[W   ])[kw   ];
		//f_S    = (D.f[N   ])[kn   ];
		//f_N    = (D.f[S   ])[ks   ];
		//f_B    = (D.f[T   ])[kt   ];
		//f_T    = (D.f[B   ])[kb   ];
		//f_SW   = (D.f[NE  ])[kne  ];
		//f_NE   = (D.f[SW  ])[ksw  ];
		//f_NW   = (D.f[SE  ])[kse  ];
		//f_SE   = (D.f[NW  ])[knw  ];
		//f_BW   = (D.f[TE  ])[kte  ];
		//f_TE   = (D.f[BW  ])[kbw  ];
		//f_TW   = (D.f[BE  ])[kbe  ];
		//f_BE   = (D.f[TW  ])[ktw  ];
		//f_BS   = (D.f[TN  ])[ktn  ];
		//f_TN   = (D.f[BS  ])[kbs  ];
		//f_TS   = (D.f[BN  ])[kbn  ];
		//f_BN   = (D.f[TS  ])[kts  ];
		//f_BSW  = (D.f[TNE ])[ktne ];
		//f_BNE  = (D.f[TSW ])[ktsw ];
		//f_BNW  = (D.f[TSE ])[ktse ];
		//f_BSE  = (D.f[TNW ])[ktnw ];
		//f_TSW  = (D.f[BNE ])[kbne ];
		//f_TNE  = (D.f[BSW ])[kbsw ];
		//f_TNW  = (D.f[BSE ])[kbse ];
		//f_TSE  = (D.f[BNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////
		real /*vx1, vx2,*/ vx3, drho, feq, q, cu_sq;
		//drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
		//	f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
		//	f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

		//vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		//	((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
		//	(f_E - f_W); 


		//vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		//	((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
		//	(f_N - f_S); 

		//vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		//	(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
		//	(f_T - f_B); 

		//cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

		//////////////////////////////////////////////////////////////////////////
		if (isEvenTimestep==false)
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
			D.f[TNE ] = &DD[BSW *size_Mat];
			D.f[TSW ] = &DD[BNE *size_Mat];
			D.f[TSE ] = &DD[BNW *size_Mat];
			D.f[TNW ] = &DD[BSE *size_Mat];
			D.f[BNE ] = &DD[TSW *size_Mat];
			D.f[BSW ] = &DD[TNE *size_Mat];
			D.f[BSE ] = &DD[TNW *size_Mat];
			D.f[BNW ] = &DD[TSE *size_Mat];
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Test
		//(D.f[dirREST])[k]=c1o10;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//ToDo anders Klammern

		//q = q_dirE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        )-cu_sq); 
		//	(D.f[W])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
		//	//(D.f[W])[kw]=zero;
		//}

		//q = q_dirW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
		//	(D.f[E])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
		//	//(D.f[E])[ke]=zero;
		//}

		//q = q_dirN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
		//	(D.f[S])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
		//	//(D.f[S])[ks]=zero;
		//}

		//q = q_dirS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
		//	(D.f[N])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
		//	//(D.f[N])[kn]=zero;
		//}

		//q = q_dirT[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3)-cu_sq); 
		//	(D.f[B])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
		//	//(D.f[B])[kb]=one;
		//}

		//q = q_dirB[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
		//	(D.f[T])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
		//	//(D.f[T])[kt]=zero;
		//}

		//q = q_dirNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
		//	(D.f[SW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
		//	//(D.f[SW])[ksw]=zero;
		//}

		//q = q_dirSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
		//	(D.f[NE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
		//	//(D.f[NE])[kne]=zero;
		//}

		//q = q_dirSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
		//	(D.f[NW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
		//	//(D.f[NW])[knw]=zero;
		//}

		//q = q_dirNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
		//	(D.f[SE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
		//	//(D.f[SE])[kse]=zero;
		//}

		//q = q_dirTE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
		//	(D.f[BW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
		//	//(D.f[BW])[kbw]=zero;
		//}

		//q = q_dirBW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
		//	(D.f[TE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
		//	//(D.f[TE])[kte]=zero;
		//}

		//q = q_dirBE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
		//	(D.f[TW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
		//	//(D.f[TW])[ktw]=zero;
		//}

		//q = q_dirTW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
		//	(D.f[BE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
		//	//(D.f[BE])[kbe]=zero;
		//}

		//q = q_dirTN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
		//	(D.f[BS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
		//	//(D.f[BS])[kbs]=zero;
		//}

		//q = q_dirBS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
		//	(D.f[TN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
		//	//(D.f[TN])[ktn]=zero;
		//}

		//q = q_dirBN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
		//	(D.f[TS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
		//	//(D.f[TS])[kts]=zero;
		//}

		//q = q_dirTS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
		//	(D.f[BN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
		//	//(D.f[BN])[kbn]=zero;
		//}

		//q = q_dirTNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
		//	(D.f[BSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[BSW])[kbsw]=zero;
		//}

		//q = q_dirBSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
		//	(D.f[TNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[TNE])[ktne]=zero;
		//}

		//q = q_dirBNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
		//	(D.f[TSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[TSW])[ktsw]=zero;
		//}

		//q = q_dirTSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
		//	(D.f[BNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[BNE])[kbne]=zero;
		//}

		//q = q_dirTSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
		//	(D.f[BNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[BNW])[kbnw]=zero;
		//}

		//q = q_dirBNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
		//	(D.f[TSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[TSE])[ktse]=zero;
		//}

		//q = q_dirBSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
		//	(D.f[TNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[TNW])[ktnw]=zero;
		//}

		//q = q_dirTNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
		//	(D.f[BSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[BSE])[kbse]=zero;
		//}

		///////// equilibrium BC
		cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
		VeloXpur*=-c1o1;
		VeloYpur*=-c1o1;
		vx3=c0o1;
		drho=c0o1;
		q = q_dirE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*( VeloXpur        )+c9o2*( VeloX        )*( VeloX        )-cu_sq); 
			(D.f[W])[kw]=feq;
			//(D.f[W])[kw]=zero;
		}

		q = q_dirW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(-VeloXpur        )+c9o2*(-VeloX        )*(-VeloX        )-cu_sq); 
			(D.f[E])[ke]=feq;
			//(D.f[E])[ke]=zero;
		}

		q = q_dirN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(    VeloYpur     )+c9o2*(     VeloY    )*(     VeloY    )-cu_sq); 
			(D.f[S])[ks]=feq;
			//(D.f[S])[ks]=zero;
		}

		q = q_dirS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(   -VeloYpur     )+c9o2*(    -VeloY    )*(    -VeloY    )-cu_sq); 
			(D.f[N])[kn]=feq;
			//(D.f[N])[kn]=zero;
		}

		q = q_dirT[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); 
			(D.f[B])[kb]=feq;
			//(D.f[B])[kb]=one;
		}

		q = q_dirB[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
			(D.f[T])[kt]=feq;
			//(D.f[T])[kt]=zero;
		}

		q = q_dirNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur+VeloYpur    )+c9o2*( VeloX+VeloY    )*( VeloX+VeloY    )-cu_sq); 
			(D.f[SW])[ksw]=feq;
			//(D.f[SW])[ksw]=zero;
		}

		q = q_dirSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur-VeloYpur    )+c9o2*(-VeloX-VeloY    )*(-VeloX-VeloY    )-cu_sq); 
			(D.f[NE])[kne]=feq;
			//(D.f[NE])[kne]=zero;
		}

		q = q_dirSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur-VeloYpur    )+c9o2*( VeloX-VeloY    )*( VeloX-VeloY    )-cu_sq); 
			(D.f[NW])[knw]=feq;
			//(D.f[NW])[knw]=zero;
		}

		q = q_dirNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur+VeloYpur    )+c9o2*(-VeloX+VeloY    )*(-VeloX+VeloY    )-cu_sq); 
			(D.f[SE])[kse]=feq;
			//(D.f[SE])[kse]=zero;
		}

		q = q_dirTE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur    +vx3)+c9o2*( VeloX    +vx3)*( VeloX    +vx3)-cu_sq); 
			(D.f[BW])[kbw]=feq;
			//(D.f[BW])[kbw]=zero;
		}

		q = q_dirBW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur    -vx3)+c9o2*(-VeloX    -vx3)*(-VeloX    -vx3)-cu_sq); 
			(D.f[TE])[kte]=feq;
			//(D.f[TE])[kte]=zero;
		}

		q = q_dirBE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*( VeloXpur    -vx3)+c9o2*( VeloX    -vx3)*( VeloX    -vx3)-cu_sq); 
			(D.f[TW])[ktw]=feq;
			//(D.f[TW])[ktw]=zero;
		}

		q = q_dirTW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(-VeloXpur    +vx3)+c9o2*(-VeloX    +vx3)*(-VeloX    +vx3)-cu_sq); 
			(D.f[BE])[kbe]=feq;
			//(D.f[BE])[kbe]=zero;
		}

		q = q_dirTN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(     VeloYpur+vx3)+c9o2*(     VeloY+vx3)*(     VeloY+vx3)-cu_sq); 
			(D.f[BS])[kbs]=feq;
			//(D.f[BS])[kbs]=zero;
		}

		q = q_dirBS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(    -VeloYpur-vx3)+c9o2*(    -VeloY-vx3)*(    -VeloY-vx3)-cu_sq); 
			(D.f[TN])[ktn]=feq;
			//(D.f[TN])[ktn]=zero;
		}

		q = q_dirBN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(     VeloYpur-vx3)+c9o2*(     VeloY-vx3)*(     VeloY-vx3)-cu_sq); 
			(D.f[TS])[kts]=feq;
			//(D.f[TS])[kts]=zero;
		}

		q = q_dirTS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c3o1*(    -VeloYpur+vx3)+c9o2*(    -VeloY+vx3)*(    -VeloY+vx3)-cu_sq); 
			(D.f[BN])[kbn]=feq;
			//(D.f[BN])[kbn]=zero;
		}

		q = q_dirTNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur+vx3)+c9o2*( VeloX+VeloY+vx3)*( VeloX+VeloY+vx3)-cu_sq); 
			(D.f[BSW])[kbsw]=feq;
			//(D.f[BSW])[kbsw]=zero;
		}

		q = q_dirBSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur-vx3)+c9o2*(-VeloX-VeloY-vx3)*(-VeloX-VeloY-vx3)-cu_sq); 
			(D.f[TNE])[ktne]=feq;
			//(D.f[TNE])[ktne]=zero;
		}

		q = q_dirBNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur+VeloYpur-vx3)+c9o2*( VeloX+VeloY-vx3)*( VeloX+VeloY-vx3)-cu_sq); 
			(D.f[TSW])[ktsw]=feq;
			//(D.f[TSW])[ktsw]=zero;
		}

		q = q_dirTSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur-VeloYpur+vx3)+c9o2*(-VeloX-VeloY+vx3)*(-VeloX-VeloY+vx3)-cu_sq); 
			(D.f[BNE])[kbne]=feq;
			//(D.f[BNE])[kbne]=zero;
		}

		q = q_dirTSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur+vx3)+c9o2*( VeloX-VeloY+vx3)*( VeloX-VeloY+vx3)-cu_sq); 
			(D.f[BNW])[kbnw]=feq;
			//(D.f[BNW])[kbnw]=zero;
		}

		q = q_dirBNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur-vx3)+c9o2*(-VeloX+VeloY-vx3)*(-VeloX+VeloY-vx3)-cu_sq); 
			(D.f[TSE])[ktse]=feq;
			//(D.f[TSE])[ktse]=zero;
		}

		q = q_dirBSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*( VeloXpur-VeloYpur-vx3)+c9o2*( VeloX-VeloY-vx3)*( VeloX-VeloY-vx3)-cu_sq); 
			(D.f[TNW])[ktnw]=feq;
			//(D.f[TNW])[ktnw]=zero;
		}

		q = q_dirTNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c3o1*(-VeloXpur+VeloYpur+vx3)+c9o2*(-VeloX+VeloY+vx3)*(-VeloX+VeloY+vx3)-cu_sq); 
			(D.f[BSE])[kbse]=feq;
			//(D.f[BSE])[kbse]=zero;
		}
	
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceComp27(
											real* velocityX,
											real* velocityY,
											real* velocityZ,
											real* distributions,
											int* subgridDistanceIndices,
											real* subgridDistances,
											unsigned int numberOfBCnodes,
											real omega,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int numberOfLBnodes,
											bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // global x-index 
   const unsigned  y = blockIdx.x;   // global y-index 
   const unsigned  z = blockIdx.y;   // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(k < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local velocities
      //!
      real VeloX = velocityX[k];
      real VeloY = velocityY[k];
      real VeloZ = velocityZ[k];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[k];
      unsigned int kzero= indexOfBCnode;
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[E   ])[ke   ];
      real f_E    = (dist.f[W   ])[kw   ];
      real f_S    = (dist.f[N   ])[kn   ];
      real f_N    = (dist.f[S   ])[ks   ];
      real f_B    = (dist.f[T   ])[kt   ];
      real f_T    = (dist.f[B   ])[kb   ];
      real f_SW   = (dist.f[NE  ])[kne  ];
      real f_NE   = (dist.f[SW  ])[ksw  ];
      real f_NW   = (dist.f[SE  ])[kse  ];
      real f_SE   = (dist.f[NW  ])[knw  ];
      real f_BW   = (dist.f[TE  ])[kte  ];
      real f_TE   = (dist.f[BW  ])[kbw  ];
      real f_TW   = (dist.f[BE  ])[kbe  ];
      real f_BE   = (dist.f[TW  ])[ktw  ];
      real f_BS   = (dist.f[TN  ])[ktn  ];
      real f_TN   = (dist.f[BS  ])[kbs  ];
      real f_TS   = (dist.f[BN  ])[kbn  ];
      real f_BN   = (dist.f[TS  ])[kts  ];
      real f_BSW  = (dist.f[TNE ])[ktne ];
      real f_BNE  = (dist.f[TSW ])[ktsw ];
      real f_BNW  = (dist.f[TSE ])[ktse ];
      real f_BSE  = (dist.f[TNW ])[ktnw ];
      real f_TSW  = (dist.f[BNE ])[kbne ];
      real f_TNE  = (dist.f[BSW ])[kbsw ];
      real f_TNW  = (dist.f[BSE ])[kbse ];
      real f_TSE  = (dist.f[BNW ])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[dirREST])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);          

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho); 

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;
      q = (subgridD.q[E])[k];
      if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[W])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[W])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[E])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[N])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[S])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[S])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[N])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[T])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[B])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[B])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[T])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[NE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[SW])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[SW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[NE])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[SE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[NW])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[NW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[SE])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[BW])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[TE])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[TW])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[BE])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[BS])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[TN])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[BN])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[TS])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TS])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[BN])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[TNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[BSW])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[TNE])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BNE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[TSW])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TSW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[BNE])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[BNW])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[TSE])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[BSE])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[TNW])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[TNW])[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[BSE])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, omega, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevice27(int inx,
                                        int iny,
                                        real* vx,
                                        real* vy,
                                        real* vz,
                                        real* DD, 
                                        int* k_Q, 
                                        real* QQ,
                                        unsigned int numberOfBCnodes, 
                                        real om1, 
                                        unsigned int* neighborX,
                                        unsigned int* neighborY,
                                        unsigned int* neighborZ,
                                        unsigned int size_Mat, 
                                        bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
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
      D.f[dirREST] = &DD[dirREST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real VeloX = vx[k];
      real VeloY = vy[k];
      real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      //unsigned int nxny = nx*ny;
      //unsigned int kzero= KQK;
      //unsigned int ke   = KQK;
      //unsigned int kw   = KQK + 1;
      //unsigned int kn   = KQK;
      //unsigned int ks   = KQK + nx;
      //unsigned int kt   = KQK;
      //unsigned int kb   = KQK + nxny;
      //unsigned int ksw  = KQK + nx + 1;
      //unsigned int kne  = KQK;
      //unsigned int kse  = KQK + nx;
      //unsigned int knw  = KQK + 1;
      //unsigned int kbw  = KQK + nxny + 1;
      //unsigned int kte  = KQK;
      //unsigned int kbe  = KQK + nxny;
      //unsigned int ktw  = KQK + 1;
      //unsigned int kbs  = KQK + nxny + nx;
      //unsigned int ktn  = KQK;
      //unsigned int kbn  = KQK + nxny;
      //unsigned int kts  = KQK + nx;
      //unsigned int ktse = KQK + nx;
      //unsigned int kbnw = KQK + nxny + 1;
      //unsigned int ktnw = KQK + 1;
      //unsigned int kbse = KQK + nxny + nx;
      //unsigned int ktsw = KQK + nx + 1;
      //unsigned int kbne = KQK + nxny;
      //unsigned int ktne = KQK;
      //unsigned int kbsw = KQK + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[E   ])[ke   ];
      f_E    = (D.f[W   ])[kw   ];
      f_S    = (D.f[N   ])[kn   ];
      f_N    = (D.f[S   ])[ks   ];
      f_B    = (D.f[T   ])[kt   ];
      f_T    = (D.f[B   ])[kb   ];
      f_SW   = (D.f[NE  ])[kne  ];
      f_NE   = (D.f[SW  ])[ksw  ];
      f_NW   = (D.f[SE  ])[kse  ];
      f_SE   = (D.f[NW  ])[knw  ];
      f_BW   = (D.f[TE  ])[kte  ];
      f_TE   = (D.f[BW  ])[kbw  ];
      f_TW   = (D.f[BE  ])[kbe  ];
      f_BE   = (D.f[TW  ])[ktw  ];
      f_BS   = (D.f[TN  ])[ktn  ];
      f_TN   = (D.f[BS  ])[kbs  ];
      f_TS   = (D.f[BN  ])[kbn  ];
      f_BN   = (D.f[TS  ])[kts  ];
      f_BSW  = (D.f[TNE ])[ktne ];
      f_BNE  = (D.f[TSW ])[ktsw ];
      f_BNW  = (D.f[TSE ])[ktse ];
      f_BSE  = (D.f[TNW ])[ktnw ];
      f_TSW  = (D.f[BNE ])[kbne ];
      f_TNE  = (D.f[BSW ])[kbsw ];
      f_TNW  = (D.f[BSE ])[kbse ];
      f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirREST])[kzero]); 

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 
         

      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
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
         D.f[dirREST] = &DD[dirREST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirREST])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);
         //(D.f[W])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         //(D.f[E])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);
         //(D.f[S])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         //(D.f[N])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         //(D.f[B])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         //(D.f[T])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         //(D.f[SW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         //(D.f[NE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         //(D.f[NW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         //(D.f[SE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         //(D.f[BW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         //(D.f[TE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         //(D.f[TW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         //(D.f[BE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         //(D.f[BS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         //(D.f[TN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         //(D.f[TS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         //(D.f[BN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[BSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[TNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[TSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[BNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[BNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[TSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[TNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[BSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void PropellerBC(unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       real* rho,
                                       real* ux,
                                       real* uy,
                                       real* uz,
                                       int* k_Q, 
									   unsigned int size_Prop,
                                       unsigned int size_Mat,
                                       unsigned int* bcMatD,
                                       real* DD,
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

   if(k<size_Prop)
   {
    ////////////////////////////////////////////////////////////////////////////////
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
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
			D.f[dirREST] = &DD[dirREST*size_Mat];
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
		unsigned int KQK = k_Q[k];
		unsigned int BC  = bcMatD[KQK];
		if( (BC != GEO_SOLID) && (BC != GEO_VOID))
		{		
		//////////////////////////////////////////////////////////////////////////
        real  vx1 = ux[k];
        real  vx2 = uy[k];
        real  vx3 = uz[k];
        //real  vx1 = -c1o100;
        //real  vx2 = zero;
        //real  vx3 = zero;
        //////////////////////////////////////////////////////////////////////////
        //index
        //////////////////////////////////////////////////////////////////////////
		unsigned int kzero= KQK;
		unsigned int ke   = KQK;
		unsigned int kw   = neighborX[KQK];
		unsigned int kn   = KQK;
		unsigned int ks   = neighborY[KQK];
		unsigned int kt   = KQK;
		unsigned int kb   = neighborZ[KQK];
		unsigned int ksw  = neighborY[kw];
		unsigned int kne  = KQK;
		unsigned int kse  = ks;
		unsigned int knw  = kw;
		unsigned int kbw  = neighborZ[kw];
		unsigned int kte  = KQK;
		unsigned int kbe  = kb;
		unsigned int ktw  = kw;
		unsigned int kbs  = neighborZ[ks];
		unsigned int ktn  = KQK;
		unsigned int kbn  = kb;
		unsigned int kts  = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
		real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

		f_ZERO= (D.f[dirREST])[kzero];
		f_E   = (D.f[E   ])[ke   ];
		f_W   = (D.f[W   ])[kw   ];
		f_N   = (D.f[N   ])[kn   ];
		f_S   = (D.f[S   ])[ks   ];
		f_T   = (D.f[T   ])[kt   ];
		f_B   = (D.f[B   ])[kb   ];
		f_NE  = (D.f[NE  ])[kne  ];
		f_SW  = (D.f[SW  ])[ksw  ];
		f_SE  = (D.f[SE  ])[kse  ];
		f_NW  = (D.f[NW  ])[knw  ];
		f_TE  = (D.f[TE  ])[kte  ];
		f_BW  = (D.f[BW  ])[kbw  ];
		f_BE  = (D.f[BE  ])[kbe  ];
		f_TW  = (D.f[TW  ])[ktw  ];
		f_TN  = (D.f[TN  ])[ktn  ];
		f_BS  = (D.f[BS  ])[kbs  ];
		f_BN  = (D.f[BN  ])[kbn  ];
		f_TS  = (D.f[TS  ])[kts  ];
		f_TNE = (D.f[TNE ])[ktne ];
		f_BSW = (D.f[BSW ])[kbsw ];
		f_BNE = (D.f[BNE ])[kbne ];
		f_TSW = (D.f[TSW ])[ktsw ];
		f_TSE = (D.f[TSE ])[ktse ];
		f_BNW = (D.f[BNW ])[kbnw ];
		f_BSE = (D.f[BSE ])[kbse ];
		f_TNW = (D.f[TNW ])[ktnw ];
		//f_W    = (D.f[E   ])[ke   ];
		//f_E    = (D.f[W   ])[kw   ];
		//f_S    = (D.f[N   ])[kn   ];
		//f_N    = (D.f[S   ])[ks   ];
		//f_B    = (D.f[T   ])[kt   ];
		//f_T    = (D.f[B   ])[kb   ];
		//f_SW   = (D.f[NE  ])[kne  ];
		//f_NE   = (D.f[SW  ])[ksw  ];
		//f_NW   = (D.f[SE  ])[kse  ];
		//f_SE   = (D.f[NW  ])[knw  ];
		//f_BW   = (D.f[TE  ])[kte  ];
		//f_TE   = (D.f[BW  ])[kbw  ];
		//f_TW   = (D.f[BE  ])[kbe  ];
		//f_BE   = (D.f[TW  ])[ktw  ];
		//f_BS   = (D.f[TN  ])[ktn  ];
		//f_TN   = (D.f[BS  ])[kbs  ];
		//f_TS   = (D.f[BN  ])[kbn  ];
		//f_BN   = (D.f[TS  ])[kts  ];
		//f_BSW  = (D.f[TNE ])[ktne ];
		//f_TNE  = (D.f[BSW ])[kbsw ];
		//f_TSW  = (D.f[BNE ])[kbne ];
		//f_BNE  = (D.f[TSW ])[ktsw ];
		//f_BNW  = (D.f[TSE ])[ktse ];
		//f_TSE  = (D.f[BNW ])[kbnw ];
		//f_TNW  = (D.f[BSE ])[kbse ];
		//f_BSE  = (D.f[TNW ])[ktnw ];
		//////////////////////////////////////////////////////////////////////////////////
		real vxo1, vxo2, vxo3, drho;
		drho   =  /*zero;*/f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				  f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

		vxo1   =   (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
					(f_E - f_W) )/ (c1o1 + drho); 
        

		vxo2   =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
					(f_N - f_S) )/ (c1o1 + drho); 

		vxo3   =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		 			(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
					(f_T - f_B) )/ (c1o1 + drho); 

		real cusq=c3o2*(vxo1*vxo1+vxo2*vxo2+vxo3*vxo3);
		//vx1 = vx1 * two - vxo1;
		//vx2 = vx2 * two - vxo2;
		//vx3 = vx3 * two - vxo3;
		real cusq2=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         //f_ZERO = ((one+drho) * (   c8over27 *(one+(-cusq2)))) - c8over27;
         //f_E    = ((one+drho) * (   c2over27 *(one+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq2))) - c2over27 ;
         //f_W    = ((one+drho) * (   c2over27 *(one+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq2))) - c2over27 ;
         //f_N    = ((one+drho) * (   c2over27 *(one+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq2))) - c2over27 ;
         //f_S    = ((one+drho) * (   c2over27 *(one+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq2))) - c2over27 ;
         //f_T    = ((one+drho) * (   c2over27 *(one+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq2))) - c2over27 ;
         //f_B    = ((one+drho) * (   c2over27 *(one+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq2))) - c2over27 ;
         //f_NE   = ((one+drho) * (   c1over54 *(one+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq2))) - c1over54 ;
         //f_SW   = ((one+drho) * (   c1over54 *(one+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq2))) - c1over54 ;
         //f_SE   = ((one+drho) * (   c1over54 *(one+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq2))) - c1over54 ;
         //f_NW   = ((one+drho) * (   c1over54 *(one+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq2))) - c1over54 ;
         //f_TE   = ((one+drho) * (   c1over54 *(one+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq2))) - c1over54 ;
         //f_BW   = ((one+drho) * (   c1over54 *(one+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq2))) - c1over54 ;
         //f_BE   = ((one+drho) * (   c1over54 *(one+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq2))) - c1over54 ;
         //f_TW   = ((one+drho) * (   c1over54 *(one+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq2))) - c1over54 ;
         //f_TN   = ((one+drho) * (   c1over54 *(one+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq2))) - c1over54 ;
         //f_BS   = ((one+drho) * (   c1over54 *(one+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq2))) - c1over54 ;
         //f_BN   = ((one+drho) * (   c1over54 *(one+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq2))) - c1over54 ;
         //f_TS   = ((one+drho) * (   c1over54 *(one+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq2))) - c1over54 ;
         //f_TNE  = ((one+drho) * (   c1over216*(one+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq2))) - c1over216;
         //f_BSW  = ((one+drho) * (   c1over216*(one+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq2))) - c1over216;
         //f_BNE  = ((one+drho) * (   c1over216*(one+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq2))) - c1over216;
         //f_TSW  = ((one+drho) * (   c1over216*(one+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq2))) - c1over216;
         //f_TSE  = ((one+drho) * (   c1over216*(one+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq2))) - c1over216;
         //f_BNW  = ((one+drho) * (   c1over216*(one+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq2))) - c1over216;
         //f_BSE  = ((one+drho) * (   c1over216*(one+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq2))) - c1over216;
         //f_TNW  = ((one+drho) * (   c1over216*(one+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq2))) - c1over216;
         f_ZERO = f_ZERO + ((c1o1+drho) * (-  c8o27* (-cusq)																   +   c8o27* (-cusq2)));
         f_E    = f_E    + ((c1o1+drho) * (-  c2o27* (c3o1*( vxo1          )+c9o2*( vxo1          )*( vxo1          )-cusq) +   c2o27* (c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq2)));
         f_W    = f_W    + ((c1o1+drho) * (-  c2o27* (c3o1*(-vxo1          )+c9o2*(-vxo1          )*(-vxo1          )-cusq) +   c2o27* (c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq2)));
         f_N    = f_N    + ((c1o1+drho) * (-  c2o27* (c3o1*(      vxo2     )+c9o2*(      vxo2     )*(      vxo2     )-cusq) +   c2o27* (c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq2)));
         f_S    = f_S    + ((c1o1+drho) * (-  c2o27* (c3o1*(     -vxo2     )+c9o2*(     -vxo2     )*(     -vxo2     )-cusq) +   c2o27* (c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq2)));
         f_T    = f_T    + ((c1o1+drho) * (-  c2o27* (c3o1*(           vxo3)+c9o2*(           vxo3)*(           vxo3)-cusq) +   c2o27* (c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq2)));
         f_B    = f_B    + ((c1o1+drho) * (-  c2o27* (c3o1*(          -vxo3)+c9o2*(          -vxo3)*(          -vxo3)-cusq) +   c2o27* (c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq2)));
         f_NE   = f_NE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1+vxo2     )+c9o2*( vxo1+vxo2     )*( vxo1+vxo2     )-cusq) +   c1o54* (c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq2)));
         f_SW   = f_SW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1-vxo2     )+c9o2*(-vxo1-vxo2     )*(-vxo1-vxo2     )-cusq) +   c1o54* (c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq2)));
         f_SE   = f_SE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1-vxo2     )+c9o2*( vxo1-vxo2     )*( vxo1-vxo2     )-cusq) +   c1o54* (c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq2)));
         f_NW   = f_NW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1+vxo2     )+c9o2*(-vxo1+vxo2     )*(-vxo1+vxo2     )-cusq) +   c1o54* (c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq2)));
         f_TE   = f_TE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1     +vxo3)+c9o2*( vxo1     +vxo3)*( vxo1     +vxo3)-cusq) +   c1o54* (c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq2)));
         f_BW   = f_BW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1     -vxo3)+c9o2*(-vxo1     -vxo3)*(-vxo1     -vxo3)-cusq) +   c1o54* (c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq2)));
         f_BE   = f_BE   + ((c1o1+drho) * (-  c1o54* (c3o1*( vxo1     -vxo3)+c9o2*( vxo1     -vxo3)*( vxo1     -vxo3)-cusq) +   c1o54* (c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq2)));
         f_TW   = f_TW   + ((c1o1+drho) * (-  c1o54* (c3o1*(-vxo1     +vxo3)+c9o2*(-vxo1     +vxo3)*(-vxo1     +vxo3)-cusq) +   c1o54* (c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq2)));
         f_TN   = f_TN   + ((c1o1+drho) * (-  c1o54* (c3o1*(      vxo2+vxo3)+c9o2*(      vxo2+vxo3)*(      vxo2+vxo3)-cusq) +   c1o54* (c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq2)));
         f_BS   = f_BS   + ((c1o1+drho) * (-  c1o54* (c3o1*(     -vxo2-vxo3)+c9o2*(     -vxo2-vxo3)*(     -vxo2-vxo3)-cusq) +   c1o54* (c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq2)));
         f_BN   = f_BN   + ((c1o1+drho) * (-  c1o54* (c3o1*(      vxo2-vxo3)+c9o2*(      vxo2-vxo3)*(      vxo2-vxo3)-cusq) +   c1o54* (c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq2)));
         f_TS   = f_TS   + ((c1o1+drho) * (-  c1o54* (c3o1*(     -vxo2+vxo3)+c9o2*(     -vxo2+vxo3)*(     -vxo2+vxo3)-cusq) +   c1o54* (c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq2)));
         f_TNE  = f_TNE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1+vxo2+vxo3)+c9o2*( vxo1+vxo2+vxo3)*( vxo1+vxo2+vxo3)-cusq) +   c1o216*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq2)));
         f_BSW  = f_BSW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1-vxo2-vxo3)+c9o2*(-vxo1-vxo2-vxo3)*(-vxo1-vxo2-vxo3)-cusq) +   c1o216*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq2)));
         f_BNE  = f_BNE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1+vxo2-vxo3)+c9o2*( vxo1+vxo2-vxo3)*( vxo1+vxo2-vxo3)-cusq) +   c1o216*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq2)));
         f_TSW  = f_TSW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1-vxo2+vxo3)+c9o2*(-vxo1-vxo2+vxo3)*(-vxo1-vxo2+vxo3)-cusq) +   c1o216*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq2)));
         f_TSE  = f_TSE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1-vxo2+vxo3)+c9o2*( vxo1-vxo2+vxo3)*( vxo1-vxo2+vxo3)-cusq) +   c1o216*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq2)));
         f_BNW  = f_BNW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1+vxo2-vxo3)+c9o2*(-vxo1+vxo2-vxo3)*(-vxo1+vxo2-vxo3)-cusq) +   c1o216*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq2)));
         f_BSE  = f_BSE  + ((c1o1+drho) * (-  c1o216*(c3o1*( vxo1-vxo2-vxo3)+c9o2*( vxo1-vxo2-vxo3)*( vxo1-vxo2-vxo3)-cusq) +   c1o216*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq2)));
         f_TNW  = f_TNW  + ((c1o1+drho) * (-  c1o216*(c3o1*(-vxo1+vxo2+vxo3)+c9o2*(-vxo1+vxo2+vxo3)*(-vxo1+vxo2+vxo3)-cusq) +   c1o216*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq2)));

		(D.f[dirREST])[kzero] =  f_ZERO;
        (D.f[E   ])[ke   ] =  f_E   ;	// f_W   ;//    	
        (D.f[W   ])[kw   ] =  f_W   ;	// f_E   ;//    	
        (D.f[N   ])[kn   ] =  f_N   ;	// f_S   ;//    	
        (D.f[S   ])[ks   ] =  f_S   ;	// f_N   ;//    	
        (D.f[T   ])[kt   ] =  f_T   ;	// f_B   ;//    	
        (D.f[B   ])[kb   ] =  f_B   ;	// f_T   ;//    	
        (D.f[NE  ])[kne  ] =  f_NE  ;	// f_SW  ;//    	
        (D.f[SW  ])[ksw  ] =  f_SW  ;	// f_NE  ;//    	
        (D.f[SE  ])[kse  ] =  f_SE  ;	// f_NW  ;//    	
        (D.f[NW  ])[knw  ] =  f_NW  ;	// f_SE  ;//    	
        (D.f[TE  ])[kte  ] =  f_TE  ;	// f_BW  ;//    	
        (D.f[BW  ])[kbw  ] =  f_BW  ;	// f_TE  ;//    	
        (D.f[BE  ])[kbe  ] =  f_BE  ;	// f_TW  ;//    	
        (D.f[TW  ])[ktw  ] =  f_TW  ;	// f_BE  ;//    	
        (D.f[TN  ])[ktn  ] =  f_TN  ;	// f_BS  ;//    	
        (D.f[BS  ])[kbs  ] =  f_BS  ;	// f_TN  ;//    	
        (D.f[BN  ])[kbn  ] =  f_BN  ;	// f_TS  ;//    	
        (D.f[TS  ])[kts  ] =  f_TS  ;	// f_BN  ;//    	
        (D.f[TNE ])[ktne ] =  f_TNE ;	// f_BSW ;//    	
        (D.f[BSW ])[kbsw ] =  f_BSW ;	// f_BNE ;//    	
        (D.f[BNE ])[kbne ] =  f_BNE ;	// f_BNW ;//    	
        (D.f[TSW ])[ktsw ] =  f_TSW ;	// f_BSE ;//    	
        (D.f[TSE ])[ktse ] =  f_TSE ;	// f_TSW ;//    	
        (D.f[BNW ])[kbnw ] =  f_BNW ;	// f_TNE ;//    	
        (D.f[BSE ])[kbse ] =  f_BSE ;	// f_TNW ;//    	
        (D.f[TNW ])[ktnw ] =  f_TNW ;	// f_TSE ;//    	

		//////////////////////////////////////////////////////////////////////////
        ////(D.f[dirREST])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[E   ])[ke   ] =   three*c2over27* ( vx1        );		//six
        //(D.f[W   ])[kw   ] =   three*c2over27* (-vx1        );		//six
        //(D.f[N   ])[kn   ] =   three*c2over27* (     vx2    );		//six
        //(D.f[S   ])[ks   ] =   three*c2over27* (    -vx2    );		//six
        //(D.f[T   ])[kt   ] =   three*c2over27* (         vx3);		//six
        //(D.f[B   ])[kb   ] =   three*c2over27* (        -vx3);		//six
        //(D.f[NE  ])[kne  ] =   three*c1over54* ( vx1+vx2    );		//six
        //(D.f[SW  ])[ksw  ] =   three*c1over54* (-vx1-vx2    );		//six
        //(D.f[SE  ])[kse  ] =   three*c1over54* ( vx1-vx2    );		//six
        //(D.f[NW  ])[knw  ] =   three*c1over54* (-vx1+vx2    );		//six
        //(D.f[TE  ])[kte  ] =   three*c1over54* ( vx1    +vx3);		//six
        //(D.f[BW  ])[kbw  ] =   three*c1over54* (-vx1    -vx3);		//six
        //(D.f[BE  ])[kbe  ] =   three*c1over54* ( vx1    -vx3);		//six
        //(D.f[TW  ])[ktw  ] =   three*c1over54* (-vx1    +vx3);		//six
        //(D.f[TN  ])[ktn  ] =   three*c1over54* (     vx2+vx3);		//six
        //(D.f[BS  ])[kbs  ] =   three*c1over54* (    -vx2-vx3);		//six
        //(D.f[BN  ])[kbn  ] =   three*c1over54* (     vx2-vx3);		//six
        //(D.f[TS  ])[kts  ] =   three*c1over54* (    -vx2+vx3);		//six
        //(D.f[TNE ])[ktne ] =   three*c1over216*( vx1+vx2+vx3);		//six
        //(D.f[BSW ])[kbsw ] =   three*c1over216*(-vx1-vx2-vx3);		//six
        //(D.f[BNE ])[kbne ] =   three*c1over216*( vx1+vx2-vx3);		//six
        //(D.f[TSW ])[ktsw ] =   three*c1over216*(-vx1-vx2+vx3);		//six
        //(D.f[TSE ])[ktse ] =   three*c1over216*( vx1-vx2+vx3);		//six
        //(D.f[BNW ])[kbnw ] =   three*c1over216*(-vx1+vx2-vx3);		//six
        //(D.f[BSE ])[kbse ] =   three*c1over216*( vx1-vx2-vx3);		//six
        //(D.f[TNW ])[ktnw ] =   three*c1over216*(-vx1+vx2+vx3);		//six
        //(D.f[dirREST])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[E   ])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
        //(D.f[W   ])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
        //(D.f[N   ])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
        //(D.f[S   ])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
        //(D.f[T   ])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
        //(D.f[B   ])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
        //(D.f[NE  ])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
        //(D.f[SW  ])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
        //(D.f[SE  ])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
        //(D.f[NW  ])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
        //(D.f[TE  ])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
        //(D.f[BW  ])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
        //(D.f[BE  ])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
        //(D.f[TW  ])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
        //(D.f[TN  ])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
        //(D.f[BS  ])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
        //(D.f[BN  ])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
        //(D.f[TS  ])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
        //(D.f[TNE ])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
        //(D.f[BSW ])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
        //(D.f[BNE ])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
        //(D.f[TSW ])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
        //(D.f[TSE ])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
        //(D.f[BNW ])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
        //(D.f[BSE ])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
        //(D.f[TNW ])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
		}
    }
}
//////////////////////////////////////////////////////////////////////////


