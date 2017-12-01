/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompPlusSlip27(int inx,
													int iny,
													doubflo* vx,
													doubflo* vy,
													doubflo* vz,
													doubflo* DD, 
													int* k_Q, 
													doubflo* QQ,
													unsigned int sizeQ,
													int kQ, 
													doubflo om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
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

   if (k < kQ)
   {
	   ////////////////////////////////////////////////////////////////////////////////
	   doubflo VeloX = vx[k];
	   doubflo VeloY = vy[k];
	   doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
	   ////////////////////////////////////////////////////////////////////////////////
	   doubflo *q_dirE, *q_dirW, *q_dirN, *q_dirS, *q_dirT, *q_dirB,
		   *q_dirNE, *q_dirSW, *q_dirSE, *q_dirNW, *q_dirTE, *q_dirBW,
		   *q_dirBE, *q_dirTW, *q_dirTN, *q_dirBS, *q_dirBN, *q_dirTS,
		   *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
		   *q_dirBSE, *q_dirBNW;
	   q_dirE = &QQ[dirE   *sizeQ];
	   q_dirW = &QQ[dirW   *sizeQ];
	   q_dirN = &QQ[dirN   *sizeQ];
	   q_dirS = &QQ[dirS   *sizeQ];
	   q_dirT = &QQ[dirT   *sizeQ];
	   q_dirB = &QQ[dirB   *sizeQ];
	   q_dirNE = &QQ[dirNE  *sizeQ];
	   q_dirSW = &QQ[dirSW  *sizeQ];
	   q_dirSE = &QQ[dirSE  *sizeQ];
	   q_dirNW = &QQ[dirNW  *sizeQ];
	   q_dirTE = &QQ[dirTE  *sizeQ];
	   q_dirBW = &QQ[dirBW  *sizeQ];
	   q_dirBE = &QQ[dirBE  *sizeQ];
	   q_dirTW = &QQ[dirTW  *sizeQ];
	   q_dirTN = &QQ[dirTN  *sizeQ];
	   q_dirBS = &QQ[dirBS  *sizeQ];
	   q_dirBN = &QQ[dirBN  *sizeQ];
	   q_dirTS = &QQ[dirTS  *sizeQ];
	   q_dirTNE = &QQ[dirTNE *sizeQ];
	   q_dirTSW = &QQ[dirTSW *sizeQ];
	   q_dirTSE = &QQ[dirTSE *sizeQ];
	   q_dirTNW = &QQ[dirTNW *sizeQ];
	   q_dirBNE = &QQ[dirBNE *sizeQ];
	   q_dirBSW = &QQ[dirBSW *sizeQ];
	   q_dirBSE = &QQ[dirBSE *sizeQ];
	   q_dirBNW = &QQ[dirBNW *sizeQ];
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
	   doubflo f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE,
		   f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

	   f_W = (D.f[dirE])[ke];
	   f_E = (D.f[dirW])[kw];
	   f_S = (D.f[dirN])[kn];
	   f_N = (D.f[dirS])[ks];
	   f_B = (D.f[dirT])[kt];
	   f_T = (D.f[dirB])[kb];
	   f_SW = (D.f[dirNE])[kne];
	   f_NE = (D.f[dirSW])[ksw];
	   f_NW = (D.f[dirSE])[kse];
	   f_SE = (D.f[dirNW])[knw];
	   f_BW = (D.f[dirTE])[kte];
	   f_TE = (D.f[dirBW])[kbw];
	   f_TW = (D.f[dirBE])[kbe];
	   f_BE = (D.f[dirTW])[ktw];
	   f_BS = (D.f[dirTN])[ktn];
	   f_TN = (D.f[dirBS])[kbs];
	   f_TS = (D.f[dirBN])[kbn];
	   f_BN = (D.f[dirTS])[kts];
	   f_BSW = (D.f[dirTNE])[ktne];
	   f_BNE = (D.f[dirTSW])[ktsw];
	   f_BNW = (D.f[dirTSE])[ktse];
	   f_BSE = (D.f[dirTNW])[ktnw];
	   f_TSW = (D.f[dirBNE])[kbne];
	   f_TNE = (D.f[dirBSW])[kbsw];
	   f_TNW = (D.f[dirBSE])[kbse];
	   f_TSE = (D.f[dirBNW])[kbnw];
	   ////////////////////////////////////////////////////////////////////////////////
	   doubflo vx1, vx2, vx3, drho, feq, q;
	   drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
		   f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
		   f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]);

	   vx1 = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		   ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
		   (f_E - f_W)) / (one + drho);


	   vx2 = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
		   ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
		   (f_N - f_S)) / (one + drho);

	   vx3 = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		   (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
		   (f_T - f_B)) / (one + drho);

	   doubflo cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3) * (one + drho);

	   //////////////////////////////////////////////////////////////////////////
	   if (evenOrOdd == false)
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
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   //Test
	   //(D.f[dirZERO])[k]=c1o10;
	   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   //ToDo anders Klammern

	   /////To Slip Or Not To Slip?
	   // We assume slip BC if |vec(V_BC)|=1. To avoid problems we take V_BC*V_BC>0.99 (c99o100)
	   if (VeloX*VeloX + VeloY*VeloY + VeloZ*VeloZ > c99o100)
		{
		   // vt=v-(n \dot v) *n
		   // n=(VeloX,VeloY,VeloZ) a misuse of the velocity variable!
		   doubflo normalV = VeloX*vx1 + VeloY*vx2 + VeloZ*vx3;
		   vx1 = vx1 - normalV*VeloX;
		   vx2 = vx2 - normalV*VeloY;
		   vx3 = vx3 - normalV*VeloZ;
		}
	  ////////////////

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);// - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);// - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);// - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);// - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);// - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);// - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);// - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVeloDeviceEQ27(doubflo* VeloX,
										   doubflo* VeloY,
										   doubflo* VeloZ,
                                           doubflo* DD, 
                                           int* k_Q, 
                                           int kQ, 
                                           doubflo om1, 
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
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

   if(k<kQ)
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

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			doubflo mfcbb = (D.f[dirE   ])[ke   ];
			doubflo mfabb = (D.f[dirW   ])[kw   ];
			doubflo mfbcb = (D.f[dirN   ])[kn   ];
			doubflo mfbab = (D.f[dirS   ])[ks   ];
			doubflo mfbbc = (D.f[dirT   ])[kt   ];
			doubflo mfbba = (D.f[dirB   ])[kb   ];
			doubflo mfccb = (D.f[dirNE  ])[kne  ];
			doubflo mfaab = (D.f[dirSW  ])[ksw  ];
			doubflo mfcab = (D.f[dirSE  ])[kse  ];
			doubflo mfacb = (D.f[dirNW  ])[knw  ];
			doubflo mfcbc = (D.f[dirTE  ])[kte  ];
			doubflo mfaba = (D.f[dirBW  ])[kbw  ];
			doubflo mfcba = (D.f[dirBE  ])[kbe  ];
			doubflo mfabc = (D.f[dirTW  ])[ktw  ];
			doubflo mfbcc = (D.f[dirTN  ])[ktn  ];
			doubflo mfbaa = (D.f[dirBS  ])[kbs  ];
			doubflo mfbca = (D.f[dirBN  ])[kbn  ];
			doubflo mfbac = (D.f[dirTS  ])[kts  ];
			doubflo mfbbb = (D.f[dirZERO])[kzero];
			doubflo mfccc = (D.f[dirTNE ])[ktne ];
			doubflo mfaac = (D.f[dirTSW ])[ktsw ];
			doubflo mfcac = (D.f[dirTSE ])[ktse ];
			doubflo mfacc = (D.f[dirTNW ])[ktnw ];
			doubflo mfcca = (D.f[dirBNE ])[kbne ];
			doubflo mfaaa = (D.f[dirBSW ])[kbsw ];
			doubflo mfcaa = (D.f[dirBSE ])[kbse ];
			doubflo mfaca = (D.f[dirBNW ])[kbnw ];
			////////////////////////////////////////////////////////////////////////////////////
			doubflo rho   = (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
							 mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
							 mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + one);//!!!!Achtung + one
			////////////////////////////////////////////////////////////////////////////////////
			doubflo vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			doubflo vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			doubflo vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			if(VeloX[k]!=zero) vvx = VeloX[k];
			if(VeloY[k]!=zero) vvy = VeloY[k];
			if(VeloZ[k]!=zero) vvz = VeloZ[k];
			////////////////////////////////////////////////////////////////////////////////////
			doubflo vx2    = vvx * vvx;
			doubflo vy2    = vvy * vvy;
			doubflo vz2    = vvz * vvz;
			////////////////////////////////////////////////////////////////////////////////////
            doubflo XXb    = -c2o3 + vx2;
            doubflo XXc    = -c1o2 * (XXb + one + vvx);
            doubflo XXa    = XXc + vvx;
            doubflo YYb    = -c2o3 + vy2;
            doubflo YYc    = -c1o2 * (YYb + one + vvy);
            doubflo YYa    = YYc + vvy;
            doubflo ZZb    = -c2o3 + vz2;
            doubflo ZZc    = -c1o2 * (ZZb + one + vvz);
            doubflo ZZa    = ZZc + vvz;
			////////////////////////////////////////////////////////////////////////////////////
            mfcbb = -rho * XXc * YYb * ZZb - c2over27 ; 
			mfabb = -rho * XXa * YYb * ZZb - c2over27 ;
			mfbcb = -rho * XXb * YYc * ZZb - c2over27 ;
			mfbab = -rho * XXb * YYa * ZZb - c2over27 ;
			mfbbc = -rho * XXb * YYb * ZZc - c2over27 ;
			mfbba = -rho * XXb * YYb * ZZa - c2over27 ;
			mfccb = -rho * XXc * YYc * ZZb - c1over54 ;
			mfaab = -rho * XXa * YYa * ZZb - c1over54 ;
			mfcab = -rho * XXc * YYa * ZZb - c1over54 ;
			mfacb = -rho * XXa * YYc * ZZb - c1over54 ;
			mfcbc = -rho * XXc * YYb * ZZc - c1over54 ;
			mfaba = -rho * XXa * YYb * ZZa - c1over54 ;
			mfcba = -rho * XXc * YYb * ZZa - c1over54 ;
			mfabc = -rho * XXa * YYb * ZZc - c1over54 ;
			mfbcc = -rho * XXb * YYc * ZZc - c1over54 ;
			mfbaa = -rho * XXb * YYa * ZZa - c1over54 ;
			mfbca = -rho * XXb * YYc * ZZa - c1over54 ;
			mfbac = -rho * XXb * YYa * ZZc - c1over54 ;
			mfbbb = -rho * XXb * YYb * ZZb - c8over27 ;
			mfccc = -rho * XXc * YYc * ZZc - c1over216;
			mfaac = -rho * XXa * YYa * ZZc - c1over216;
			mfcac = -rho * XXc * YYa * ZZc - c1over216;
			mfacc = -rho * XXa * YYc * ZZc - c1over216;
			mfcca = -rho * XXc * YYc * ZZa - c1over216;
			mfaaa = -rho * XXa * YYa * ZZa - c1over216;
			mfcaa = -rho * XXc * YYa * ZZa - c1over216;
			mfaca = -rho * XXa * YYc * ZZa - c1over216;
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			(D.f[dirE   ])[ke   ] = mfabb;//mfcbb;
			(D.f[dirW   ])[kw   ] = mfcbb;//mfabb;
			(D.f[dirN   ])[kn   ] = mfbab;//mfbcb;
			(D.f[dirS   ])[ks   ] = mfbcb;//mfbab;
			(D.f[dirT   ])[kt   ] = mfbba;//mfbbc;
			(D.f[dirB   ])[kb   ] = mfbbc;//mfbba;
			(D.f[dirNE  ])[kne  ] = mfaab;//mfccb;
			(D.f[dirSW  ])[ksw  ] = mfccb;//mfaab;
			(D.f[dirSE  ])[kse  ] = mfacb;//mfcab;
			(D.f[dirNW  ])[knw  ] = mfcab;//mfacb;
			(D.f[dirTE  ])[kte  ] = mfaba;//mfcbc;
			(D.f[dirBW  ])[kbw  ] = mfcbc;//mfaba;
			(D.f[dirBE  ])[kbe  ] = mfabc;//mfcba;
			(D.f[dirTW  ])[ktw  ] = mfcba;//mfabc;
			(D.f[dirTN  ])[ktn  ] = mfbaa;//mfbcc;
			(D.f[dirBS  ])[kbs  ] = mfbcc;//mfbaa;
			(D.f[dirBN  ])[kbn  ] = mfbac;//mfbca;
			(D.f[dirTS  ])[kts  ] = mfbca;//mfbac;
			(D.f[dirZERO])[kzero] = mfbbb;//mfbbb;
			(D.f[dirTNE ])[ktne ] = mfaaa;//mfccc;
			(D.f[dirTSW ])[ktsw ] = mfcca;//mfaac;
			(D.f[dirTSE ])[ktse ] = mfaca;//mfcac;
			(D.f[dirTNW ])[ktnw ] = mfcaa;//mfacc;
			(D.f[dirBNE ])[kbne ] = mfaac;//mfcca;
			(D.f[dirBSW ])[kbsw ] = mfccc;//mfaaa;
			(D.f[dirBSE ])[kbse ] = mfacc;//mfcaa;
			(D.f[dirBNW ])[kbnw ] = mfcac;//mfaca;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceIncompHighNu27(int inx,
													int iny,
													doubflo* vx,
													doubflo* vy,
													doubflo* vz,
													doubflo* DD, 
													int* k_Q, 
													doubflo* QQ,
													unsigned int sizeQ,
													int kQ, 
													doubflo om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
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

   if(k<kQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      doubflo VeloX = vx[k];
      doubflo VeloY = vy[k];
      doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W));// / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S));// / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B));// / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);// * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirW])[kw]=((one - q) * f_E + q * ((f_E + f_W) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloX     )) / (q + one) ;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) /** (one + drho)*/-cu_sq); 
         (D.f[dirE])[ke]=((one - q) * f_W + q * ((f_W + f_E) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloX     )) / (q + one) ;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirS])[ks]=((one - q) * f_N + q * ((f_N + f_S) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloY     )) / (q + one) ;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirN])[kn]=((one - q) * f_S + q * ((f_S + f_N) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloY     )) / (q + one) ;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirB])[kb]=((one - q) * f_T + q * ((f_T + f_B) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloZ     )) / (q + one) ;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirT])[kt]=((one - q) * f_B + q * ((f_B + f_T) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloZ     )) / (q + one) ;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSW])[ksw]=((one - q) * f_NE + q * ((f_NE + f_SW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloY)) / (q + one) ;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNE])[kne]=((one - q) * f_SW + q * ((f_SW + f_NE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloY)) / (q + one) ;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirNW])[knw]=((one - q) * f_SE + q * ((f_SE + f_NW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloY)) / (q + one) ;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) /** (one + drho)*/-cu_sq); 
         (D.f[dirSE])[kse]=((one - q) * f_NW + q * ((f_NW + f_SE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloY)) / (q + one) ;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBW])[kbw]=((one - q) * f_TE + q * ((f_TE + f_BW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloZ)) / (q + one) ;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTE])[kte]=((one - q) * f_BW + q * ((f_BW + f_TE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloZ)) / (q + one) ;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTW])[ktw]=((one - q) * f_BE + q * ((f_BE + f_TW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloZ)) / (q + one) ;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBE])[kbe]=((one - q) * f_TW + q * ((f_TW + f_BE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloZ)) / (q + one) ;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBS])[kbs]=((one - q) * f_TN + q * ((f_TN + f_BS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTN])[ktn]=((one - q) * f_BS + q * ((f_BS + f_TN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTS])[kts]=((one - q) * f_BN + q * ((f_BN + f_TS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBN])[kbn]=((one - q) * f_TS + q * ((f_TS + f_BN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSW])[kbsw]=((one - q) * f_TNE + q * ((f_TNE + f_BSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNE])[ktne]=((one - q) * f_BSW + q * ((f_BSW + f_TNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSW])[ktsw]=((one - q) * f_BNE + q * ((f_BNE + f_TSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNE])[kbne]=((one - q) * f_TSW + q * ((f_TSW + f_BNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBNW])[kbnw]=((one - q) * f_TSE + q * ((f_TSE + f_BNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY+VeloZ)) / (q + one) ;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTSE])[ktse]=((one - q) * f_BNW + q * ((f_BNW + f_TSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirTNW])[ktnw]=((one - q) * f_BSE + q * ((f_BSE + f_TNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY-VeloZ)) / (q + one) ;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) /** (one + drho)*/-cu_sq); 
         (D.f[dirBSE])[kbse]=((one - q) * f_TNW + q * ((f_TNW + f_BSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY+VeloZ)) / (q + one) ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompHighNu27(  int inx,
													int iny,
													doubflo* vx,
													doubflo* vy,
													doubflo* vz,
													doubflo* DD, 
													int* k_Q, 
													doubflo* QQ,
													unsigned int sizeQ,
													int kQ, 
													doubflo om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
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

   if(k<kQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      doubflo VeloX = vx[k];
      doubflo VeloY = vy[k];
      doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_E   = (D.f[dirE   ])[ke   ];
      f_W   = (D.f[dirW   ])[kw   ];
      f_N   = (D.f[dirN   ])[kn   ];
      f_S   = (D.f[dirS   ])[ks   ];
      f_T   = (D.f[dirT   ])[kt   ];
      f_B   = (D.f[dirB   ])[kb   ];
      f_NE  = (D.f[dirNE  ])[kne  ];
      f_SW  = (D.f[dirSW  ])[ksw  ];
      f_SE  = (D.f[dirSE  ])[kse  ];
      f_NW  = (D.f[dirNW  ])[knw  ];
      f_TE  = (D.f[dirTE  ])[kte  ];
      f_BW  = (D.f[dirBW  ])[kbw  ];
      f_BE  = (D.f[dirBE  ])[kbe  ];
      f_TW  = (D.f[dirTW  ])[ktw  ];
      f_TN  = (D.f[dirTN  ])[ktn  ];
      f_BS  = (D.f[dirBS  ])[kbs  ];
      f_BN  = (D.f[dirBN  ])[kbn  ];
      f_TS  = (D.f[dirTS  ])[kts  ];
      f_TNE = (D.f[dirTNE ])[ktne ];
      f_TSW = (D.f[dirTSW ])[ktsw ];
      f_TSE = (D.f[dirTSE ])[ktse ];
      f_TNW = (D.f[dirTNW ])[ktnw ];
      f_BNE = (D.f[dirBNE ])[kbne ];
      f_BSW = (D.f[dirBSW ])[kbsw ];
      f_BSE = (D.f[dirBSE ])[kbse ];
      f_BNW = (D.f[dirBNW ])[kbnw ];
      //f_W    = (D.f[dirE   ])[ke   ];
      //f_E    = (D.f[dirW   ])[kw   ];
      //f_S    = (D.f[dirN   ])[kn   ];
      //f_N    = (D.f[dirS   ])[ks   ];
      //f_B    = (D.f[dirT   ])[kt   ];
      //f_T    = (D.f[dirB   ])[kb   ];
      //f_SW   = (D.f[dirNE  ])[kne  ];
      //f_NE   = (D.f[dirSW  ])[ksw  ];
      //f_NW   = (D.f[dirSE  ])[kse  ];
      //f_SE   = (D.f[dirNW  ])[knw  ];
      //f_BW   = (D.f[dirTE  ])[kte  ];
      //f_TE   = (D.f[dirBW  ])[kbw  ];
      //f_TW   = (D.f[dirBE  ])[kbe  ];
      //f_BE   = (D.f[dirTW  ])[ktw  ];
      //f_BS   = (D.f[dirTN  ])[ktn  ];
      //f_TN   = (D.f[dirBS  ])[kbs  ];
      //f_TS   = (D.f[dirBN  ])[kbn  ];
      //f_BN   = (D.f[dirTS  ])[kts  ];
      //f_BSW  = (D.f[dirTNE ])[ktne ];
      //f_BNE  = (D.f[dirTSW ])[ktsw ];
      //f_BNW  = (D.f[dirTSE ])[ktse ];
      //f_BSE  = (D.f[dirTNW ])[ktnw ];
      //f_TSW  = (D.f[dirBNE ])[kbne ];
      //f_TNE  = (D.f[dirBSW ])[kbsw ];
      //f_TNW  = (D.f[dirBSE ])[kbse ];
      //f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=((one - q) * f_E + q * ((f_E + f_W) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloX     )) / (q + one) ;
         //(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=((one - q) * f_W + q * ((f_W + f_E) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloX     )) / (q + one) ;
         //(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=((one - q) * f_N + q * ((f_N + f_S) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloY     )) / (q + one) ;
         //(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=((one - q) * f_S + q * ((f_S + f_N) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloY     )) / (q + one) ;
         //(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=((one - q) * f_T + q * ((f_T + f_B) * (one - om1) + om1 * two * feq) - six*c2over27*( VeloZ     )) / (q + one) ;
         //(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=((one - q) * f_B + q * ((f_B + f_T) * (one - om1) + om1 * two * feq) - six*c2over27*(-VeloZ     )) / (q + one) ;
         //(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=((one - q) * f_NE + q * ((f_NE + f_SW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloY)) / (q + one) ;
         //(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=((one - q) * f_SW + q * ((f_SW + f_NE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloY)) / (q + one) ;
         //(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=((one - q) * f_SE + q * ((f_SE + f_NW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloY)) / (q + one) ;
         //(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=((one - q) * f_NW + q * ((f_NW + f_SE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloY)) / (q + one) ;
         //(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=((one - q) * f_TE + q * ((f_TE + f_BW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX+VeloZ)) / (q + one) ;
         //(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=((one - q) * f_BW + q * ((f_BW + f_TE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX-VeloZ)) / (q + one) ;
         //(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=((one - q) * f_BE + q * ((f_BE + f_TW) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloX-VeloZ)) / (q + one) ;
         //(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=((one - q) * f_TW + q * ((f_TW + f_BE) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloX+VeloZ)) / (q + one) ;
         //(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=((one - q) * f_TN + q * ((f_TN + f_BS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=((one - q) * f_BS + q * ((f_BS + f_TN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=((one - q) * f_BN + q * ((f_BN + f_TS) * (one - om1) + om1 * two * feq) - six*c1over54*( VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=((one - q) * f_TS + q * ((f_TS + f_BN) * (one - om1) + om1 * two * feq) - six*c1over54*(-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=((one - q) * f_TNE + q * ((f_TNE + f_BSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=((one - q) * f_BSW + q * ((f_BSW + f_TNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=((one - q) * f_BNE + q * ((f_BNE + f_TSW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX+VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=((one - q) * f_TSW + q * ((f_TSW + f_BNE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=((one - q) * f_TSE + q * ((f_TSE + f_BNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=((one - q) * f_BNW + q * ((f_BNW + f_TSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=((one - q) * f_BSE + q * ((f_BSE + f_TNW) * (one - om1) + om1 * two * feq) - six*c1over216*( VeloX-VeloY-VeloZ)) / (q + one) ;
         //(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=((one - q) * f_TNW + q * ((f_TNW + f_BSE) * (one - om1) + om1 * two * feq) - six*c1over216*(-VeloX+VeloY+VeloZ)) / (q + one) ;
         //(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompZeroPress27(   int inx,
														int iny,
														doubflo* vx,
														doubflo* vy,
														doubflo* vz,
														doubflo* DD, 
														int* k_Q, 
														doubflo* QQ,
														unsigned int sizeQ,
														//int kQ, 
														doubflo om1, 
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														unsigned int size_Mat, 
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

   if(k<sizeQ/*kQ*/)
   {
      ////////////////////////////////////////////////////////////////////////////////
      doubflo VeloX = vx[k];
      doubflo VeloY = vy[k];
      doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q) - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q) - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q) - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q) - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q) - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q) - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceCompZeroPress1h27( int inx,
														int iny,
														doubflo* vx,
														doubflo* vy,
														doubflo* vz,
														doubflo* DD, 
														int* k_Q, 
														doubflo* QQ,
														unsigned int sizeQ,
														int kQ, 
														doubflo om1, 
														doubflo Phi,
														doubflo angularVelocity,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														doubflo* coordX,
														doubflo* coordY,
														doubflo* coordZ,
														unsigned int size_Mat, 
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

   if(k<kQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //doubflo VeloX = vx[k];
      //doubflo VeloY = vy[k];
      //doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////
		doubflo VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		doubflo VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		doubflo coord0X = 281.125f;//7.5f;
		doubflo coord0Y = 388.125f;//7.5f;
		doubflo ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		doubflo uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		doubflo VeloXpur=VeloX;
		doubflo VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  doubflo vx1, vx2, vx3, drho, feq, q, cu_sq;
	  ///////// equilibrium BC
	  cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
	  VeloXpur*=-one;
	  VeloYpur*=-one;
	  vx1=VeloX;
	  vx2=VeloY;
	  vx3=zero;
	  drho=zero;

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*( VeloXpur        )+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]= feq - c2over27 * drho;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*(-VeloXpur        )+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]= feq - c2over27 * drho;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*(    VeloYpur     )+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]= feq - c2over27 * drho;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*(   -VeloYpur     )+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]= feq - c2over27 * drho;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]= feq - c2over27 * drho;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]= feq - c2over27 * drho;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*( VeloXpur+VeloYpur    )+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]= feq - c1over54 * drho;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(-VeloXpur-VeloYpur    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]= feq - c1over54 * drho;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*( VeloXpur-VeloYpur    )+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]= feq - c1over54 * drho;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(-VeloXpur+VeloYpur    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]= feq - c1over54 * drho;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*( VeloXpur    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]= feq - c1over54 * drho;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(-VeloXpur    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]= feq - c1over54 * drho;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*( VeloXpur    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]= feq - c1over54 * drho;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(-VeloXpur    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]= feq - c1over54 * drho;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(     VeloYpur+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]= feq - c1over54 * drho;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(    -VeloYpur-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]= feq - c1over54 * drho;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(     VeloYpur-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]= feq - c1over54 * drho;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho+three*(    -VeloYpur+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]= feq - c1over54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*( VeloXpur+VeloYpur+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]= feq - c1over216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*(-VeloXpur-VeloYpur-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]= feq - c1over216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*( VeloXpur+VeloYpur-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]= feq - c1over216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*(-VeloXpur-VeloYpur+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]= feq - c1over216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*( VeloXpur-VeloYpur+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]= feq - c1over216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*(-VeloXpur+VeloYpur-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]= feq - c1over216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*( VeloXpur-VeloYpur-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]= feq - c1over216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho+three*(-VeloXpur+VeloYpur+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]= feq - c1over216 * drho;
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
                                              doubflo* DD, 
                                              unsigned int size_Mat, 
                                              bool evenOrOdd, 
                                              doubflo u0x, 
                                              unsigned int grid_nx, 
                                              unsigned int grid_ny, 
                                              doubflo om) 
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
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo drho = drho1;
      doubflo  vx1 = zero;
      doubflo  vx2 = zero;
      doubflo  vx3 = u0x;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
   __syncthreads();
}          
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevPlainBB27(doubflo* vx,
											doubflo* vy,
	 										doubflo* vz,
											doubflo* DD,
											int* k_Q, 
											doubflo* QQ,
											unsigned int sizeQ,
											int kQ, 
											doubflo om1, 
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
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

   if(k<kQ)
   {
	  ////////////////////////////////////////////////////////////////////////////////
	  doubflo VeloX = vx[k];
	  doubflo VeloY = vy[k];
	  doubflo VeloZ = vz[k];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo*q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			 *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			 *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			 *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			 *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
	  ////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      doubflo q;
      q = q_dirE[k];	if (q>=zero && q<=one)	(D.f[dirW  ])[kw  ]=f_E   + c2over27  * VeloX;	
      q = q_dirW[k];	if (q>=zero && q<=one)	(D.f[dirE  ])[ke  ]=f_W   - c2over27  * VeloX;	
      q = q_dirN[k];	if (q>=zero && q<=one)	(D.f[dirS  ])[ks  ]=f_N   + c2over27  * VeloY;	
      q = q_dirS[k];	if (q>=zero && q<=one)	(D.f[dirN  ])[kn  ]=f_S   - c2over27  * VeloY;	
      q = q_dirT[k];	if (q>=zero && q<=one)	(D.f[dirB  ])[kb  ]=f_T   + c2over27  * VeloZ;
      q = q_dirB[k];	if (q>=zero && q<=one)	(D.f[dirT  ])[kt  ]=f_B   - c2over27  * VeloZ;
      q = q_dirNE[k];	if (q>=zero && q<=one)	(D.f[dirSW ])[ksw ]=f_NE  + c1over54  * VeloX + c1over54  * VeloY;
	  q = q_dirSW[k];	if (q>=zero && q<=one)	(D.f[dirNE ])[kne ]=f_SW  - c1over54  * VeloX - c1over54  * VeloY;
	  q = q_dirSE[k];	if (q>=zero && q<=one)	(D.f[dirNW ])[knw ]=f_SE  + c1over54  * VeloX - c1over54  * VeloY;
	  q = q_dirNW[k];	if (q>=zero && q<=one)	(D.f[dirSE ])[kse ]=f_NW  - c1over54  * VeloX + c1over54  * VeloY;
	  q = q_dirTE[k];	if (q>=zero && q<=one)	(D.f[dirBW ])[kbw ]=f_TE  + c1over54  * VeloX + c1over54  * VeloZ;
	  q = q_dirBW[k];	if (q>=zero && q<=one)	(D.f[dirTE ])[kte ]=f_BW  - c1over54  * VeloX - c1over54  * VeloZ;
	  q = q_dirBE[k];	if (q>=zero && q<=one)	(D.f[dirTW ])[ktw ]=f_BE  + c1over54  * VeloX - c1over54  * VeloZ;
	  q = q_dirTW[k];	if (q>=zero && q<=one)	(D.f[dirBE ])[kbe ]=f_TW  - c1over54  * VeloX + c1over54  * VeloZ;
	  q = q_dirTN[k];	if (q>=zero && q<=one)	(D.f[dirBS ])[kbs ]=f_TN  + c1over54  * VeloY + c1over54  * VeloZ;
	  q = q_dirBS[k];	if (q>=zero && q<=one)	(D.f[dirTN ])[ktn ]=f_BS  - c1over54  * VeloY - c1over54  * VeloZ;
	  q = q_dirBN[k];	if (q>=zero && q<=one)	(D.f[dirTS ])[kts ]=f_BN  + c1over54  * VeloY - c1over54  * VeloZ;
	  q = q_dirTS[k];	if (q>=zero && q<=one)	(D.f[dirBN ])[kbn ]=f_TS  - c1over54  * VeloY + c1over54  * VeloZ;
      q = q_dirTNE[k];	if (q>=zero && q<=one)	(D.f[dirBSW])[kbsw]=f_TNE + c1over216 * VeloX + c1over216 * VeloY + c1over216 * VeloZ;
      q = q_dirBSW[k];	if (q>=zero && q<=one)	(D.f[dirTNE])[ktne]=f_BSW - c1over216 * VeloX - c1over216 * VeloY - c1over216 * VeloZ;
      q = q_dirBNE[k];	if (q>=zero && q<=one)	(D.f[dirTSW])[ktsw]=f_BNE + c1over216 * VeloX + c1over216 * VeloY - c1over216 * VeloZ;
      q = q_dirTSW[k];	if (q>=zero && q<=one)	(D.f[dirBNE])[kbne]=f_TSW - c1over216 * VeloX - c1over216 * VeloY + c1over216 * VeloZ;
      q = q_dirTSE[k];	if (q>=zero && q<=one)	(D.f[dirBNW])[kbnw]=f_TSE + c1over216 * VeloX - c1over216 * VeloY + c1over216 * VeloZ;
      q = q_dirBNW[k];	if (q>=zero && q<=one)	(D.f[dirTSE])[ktse]=f_BNW - c1over216 * VeloX + c1over216 * VeloY - c1over216 * VeloZ;
      q = q_dirBSE[k];	if (q>=zero && q<=one)	(D.f[dirTNW])[ktnw]=f_BSE + c1over216 * VeloX - c1over216 * VeloY - c1over216 * VeloZ;
      q = q_dirTNW[k];	if (q>=zero && q<=one)	(D.f[dirBSE])[kbse]=f_TNW - c1over216 * VeloX + c1over216 * VeloY + c1over216 * VeloZ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevCouhette27(doubflo* vx,
											doubflo* vy,
	 										doubflo* vz,
											doubflo* DD,
											int* k_Q, 
											doubflo* QQ,
											unsigned int sizeQ,
											int kQ, 
											doubflo om1, 
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
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

   if(k<kQ)
   {
	  ////////////////////////////////////////////////////////////////////////////////
	  doubflo VeloX = vx[k];
	  doubflo VeloY = vy[k];
	  doubflo VeloZ = vz[k];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo*q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			 *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			 *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			 *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			 *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_W    = (D.f[dirE   ])[ke   ];
      doubflo f_E    = (D.f[dirW   ])[kw   ];
      doubflo f_S    = (D.f[dirN   ])[kn   ];
      doubflo f_N    = (D.f[dirS   ])[ks   ];
      doubflo f_B    = (D.f[dirT   ])[kt   ];
      doubflo f_T    = (D.f[dirB   ])[kb   ];
      doubflo f_SW   = (D.f[dirNE  ])[kne  ];
      doubflo f_NE   = (D.f[dirSW  ])[ksw  ];
      doubflo f_NW   = (D.f[dirSE  ])[kse  ];
      doubflo f_SE   = (D.f[dirNW  ])[knw  ];
      doubflo f_BW   = (D.f[dirTE  ])[kte  ];
      doubflo f_TE   = (D.f[dirBW  ])[kbw  ];
      doubflo f_TW   = (D.f[dirBE  ])[kbe  ];
      doubflo f_BE   = (D.f[dirTW  ])[ktw  ];
      doubflo f_BS   = (D.f[dirTN  ])[ktn  ];
      doubflo f_TN   = (D.f[dirBS  ])[kbs  ];
      doubflo f_TS   = (D.f[dirBN  ])[kbn  ];
      doubflo f_BN   = (D.f[dirTS  ])[kts  ];
      doubflo f_BSW  = (D.f[dirTNE ])[ktne ];
      doubflo f_BNE  = (D.f[dirTSW ])[ktsw ];
      doubflo f_BNW  = (D.f[dirTSE ])[ktse ];
      doubflo f_BSE  = (D.f[dirTNW ])[ktnw ];
      doubflo f_TSW  = (D.f[dirBNE ])[kbne ];
      doubflo f_TNE  = (D.f[dirBSW ])[kbsw ];
      doubflo f_TNW  = (D.f[dirBSE ])[kbse ];
      doubflo f_TSE  = (D.f[dirBNW ])[kbnw ];
	  ////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////               FlowDirection Y !!!!!!!!!!                                                           ///////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //calculate velocity
	  //doubflo vx1 = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W);
	  doubflo vx2 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S);
	  //doubflo vx3 = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B);
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //constant
	  doubflo on=zero;//c1o2;//one;
	  doubflo ms=-six;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //2nd order moment
	  doubflo kxxMyyFromfcNEQ = zero;//-c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));		//all E+W minus all N+S (no combinations of xy left)

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set distributions
      doubflo q;
      q = q_dirE[k];	if (q>=zero && q<=one)	(D.f[dirW  ])[kw  ]=f_E   + ms*c2over27  * VeloX;	
      q = q_dirW[k];	if (q>=zero && q<=one)	(D.f[dirE  ])[ke  ]=f_W   - ms*c2over27  * VeloX;	
      q = q_dirN[k];	if (q>=zero && q<=one)	(D.f[dirS  ])[ks  ]=f_N   + ms*c2over27  * VeloY;	
      q = q_dirS[k];	if (q>=zero && q<=one)	(D.f[dirN  ])[kn  ]=f_S   - ms*c2over27  * VeloY;	
	  q = q_dirT[k];	if (q>=zero && q<=one)	(D.f[dirB  ])[kb  ]=f_T   + ms*c2over27  * VeloZ - c3o2*c2over27*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirB[k];	if (q>=zero && q<=one)	(D.f[dirT  ])[kt  ]=f_B   - ms*c2over27  * VeloZ;
      q = q_dirNE[k];	if (q>=zero && q<=one)	(D.f[dirSW ])[ksw ]=f_NE  + ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  q = q_dirSW[k];	if (q>=zero && q<=one)	(D.f[dirNE ])[kne ]=f_SW  - ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  q = q_dirSE[k];	if (q>=zero && q<=one)	(D.f[dirNW ])[knw ]=f_SE  + ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  q = q_dirNW[k];	if (q>=zero && q<=one)	(D.f[dirSE ])[kse ]=f_NW  - ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  q = q_dirTE[k];	if (q>=zero && q<=one)	(D.f[dirBW ])[kbw ]=f_TE  + ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c3o2*c1over54*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBW[k];	if (q>=zero && q<=one)	(D.f[dirTE ])[kte ]=f_BW  - ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  q = q_dirBE[k];	if (q>=zero && q<=one)	(D.f[dirTW ])[ktw ]=f_BE  + ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  q = q_dirTW[k];	if (q>=zero && q<=one)	(D.f[dirBE ])[kbe ]=f_TW  - ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c3o2*c1over54*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on-c1o12*kxxMyyFromfcNEQ;
	  q = q_dirTN[k];	if (q>=zero && q<=one)	(D.f[dirBS ])[kbs ]=f_TN  + ms*c1over54  * VeloY + ms*c1over54  * VeloZ + three*c1over54*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
	  q = q_dirBS[k];	if (q>=zero && q<=one)	(D.f[dirTN ])[ktn ]=f_BS  - ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  q = q_dirBN[k];	if (q>=zero && q<=one)	(D.f[dirTS ])[kts ]=f_BN  + ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  q = q_dirTS[k];	if (q>=zero && q<=one)	(D.f[dirBN ])[kbn ]=f_TS  - ms*c1over54  * VeloY + ms*c1over54  * VeloZ + three*c1over54*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on+c1o12*kxxMyyFromfcNEQ;
      q = q_dirTNE[k];	if (q>=zero && q<=one)	(D.f[dirBSW])[kbsw]=f_TNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + three*c1over216*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBSW[k];	if (q>=zero && q<=one)	(D.f[dirTNE])[ktne]=f_BSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
      q = q_dirBNE[k];	if (q>=zero && q<=one)	(D.f[dirTSW])[ktsw]=f_BNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
      q = q_dirTSW[k];	if (q>=zero && q<=one)	(D.f[dirBNE])[kbne]=f_TSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + three*c1over216*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirTSE[k];	if (q>=zero && q<=one)	(D.f[dirBNW])[kbnw]=f_TSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + three*c1over216*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
      q = q_dirBNW[k];	if (q>=zero && q<=one)	(D.f[dirTSE])[ktse]=f_BNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
      q = q_dirBSE[k];	if (q>=zero && q<=one)	(D.f[dirTNW])[ktnw]=f_BSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
      q = q_dirTNW[k];	if (q>=zero && q<=one)	(D.f[dirBSE])[kbse]=f_TNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + three*c1over216*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
      //q = q_dirE[k];	if (q>=zero && q<=one)	(D.f[dirW  ])[kw  ]=f_E   + ms*c2over27  * VeloX;	
   //   q = q_dirW[k];	if (q>=zero && q<=one)	(D.f[dirE  ])[ke  ]=f_W   - ms*c2over27  * VeloX;	
   //   q = q_dirN[k];	if (q>=zero && q<=one)	(D.f[dirS  ])[ks  ]=f_N   + ms*c2over27  * VeloY;	
   //   q = q_dirS[k];	if (q>=zero && q<=one)	(D.f[dirN  ])[kn  ]=f_S   - ms*c2over27  * VeloY;	
	  //q = q_dirT[k];	if (q>=zero && q<=one)	(D.f[dirB  ])[kb  ]=f_T   + ms*c2over27  * VeloZ - c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirB[k];	if (q>=zero && q<=one)	(D.f[dirT  ])[kt  ]=f_B   - ms*c2over27  * VeloZ;
   //   q = q_dirNE[k];	if (q>=zero && q<=one)	(D.f[dirSW ])[ksw ]=f_NE  + ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirSW[k];	if (q>=zero && q<=one)	(D.f[dirNE ])[kne ]=f_SW  - ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirSE[k];	if (q>=zero && q<=one)	(D.f[dirNW ])[knw ]=f_SE  + ms*c1over54  * VeloX - ms*c1over54  * VeloY;
	  //q = q_dirNW[k];	if (q>=zero && q<=one)	(D.f[dirSE ])[kse ]=f_NW  - ms*c1over54  * VeloX + ms*c1over54  * VeloY;
	  //q = q_dirTE[k];	if (q>=zero && q<=one)	(D.f[dirBW ])[kbw ]=f_TE  + ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBW[k];	if (q>=zero && q<=one)	(D.f[dirTE ])[kte ]=f_BW  - ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirBE[k];	if (q>=zero && q<=one)	(D.f[dirTW ])[ktw ]=f_BE  + ms*c1over54  * VeloX - ms*c1over54  * VeloZ;
	  //q = q_dirTW[k];	if (q>=zero && q<=one)	(D.f[dirBE ])[kbe ]=f_TW  - ms*c1over54  * VeloX + ms*c1over54  * VeloZ - c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirTN[k];	if (q>=zero && q<=one)	(D.f[dirBS ])[kbs ]=f_TN  + ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //q = q_dirBS[k];	if (q>=zero && q<=one)	(D.f[dirTN ])[ktn ]=f_BS  - ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirBN[k];	if (q>=zero && q<=one)	(D.f[dirTS ])[kts ]=f_BN  + ms*c1over54  * VeloY - ms*c1over54  * VeloZ;
	  //q = q_dirTS[k];	if (q>=zero && q<=one)	(D.f[dirBN ])[kbn ]=f_TS  - ms*c1over54  * VeloY + ms*c1over54  * VeloZ + c1o2*c1o9*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTNE[k];	if (q>=zero && q<=one)	(D.f[dirBSW])[kbsw]=f_TNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBSW[k];	if (q>=zero && q<=one)	(D.f[dirTNE])[ktne]=f_BSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBNE[k];	if (q>=zero && q<=one)	(D.f[dirTSW])[ktsw]=f_BNE + ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTSW[k];	if (q>=zero && q<=one)	(D.f[dirBNE])[kbne]=f_TSW - ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirTSE[k];	if (q>=zero && q<=one)	(D.f[dirBNW])[kbnw]=f_TSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
   //   q = q_dirBNW[k];	if (q>=zero && q<=one)	(D.f[dirTSE])[ktse]=f_BNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirBSE[k];	if (q>=zero && q<=one)	(D.f[dirTNW])[ktnw]=f_BSE + ms*c1over216 * VeloX - ms*c1over216 * VeloY - ms*c1over216 * VeloZ;
   //   q = q_dirTNW[k];	if (q>=zero && q<=one)	(D.f[dirBSE])[kbse]=f_TNW - ms*c1over216 * VeloX + ms*c1over216 * VeloY + ms*c1over216 * VeloZ + c1o2*c1o36*((two*VeloY-vx2)*(two*VeloY-vx2)-vx2*vx2)*on;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDev1h27( int inx,
										int iny,
										doubflo* vx,
										doubflo* vy,
										doubflo* vz,
										doubflo* DD, 
										int* k_Q, 
										doubflo* QQ,
										unsigned int sizeQ,
										int kQ, 
										doubflo om1,
										doubflo Phi,
										doubflo angularVelocity,
										unsigned int* neighborX,
										unsigned int* neighborY,
										unsigned int* neighborZ,
										doubflo* coordX,
										doubflo* coordY,
										doubflo* coordZ,
										unsigned int size_Mat, 
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

	if(k<kQ)
	{
		////////////////////////////////////////////////////////////////////////////////
		doubflo VeloX = cosf(Phi)*vx[k] - sinf(Phi)*vy[k];
		doubflo VeloY = sinf(Phi)*vx[k] + cosf(Phi)*vy[k];
		//doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////////
		//Ship
		doubflo coord0X = 281.125f;//7.5f;
		doubflo coord0Y = 388.125f;//7.5f;
		doubflo ux = - angularVelocity * (coordY[k_Q[k]] - coord0Y);
		doubflo uy =   angularVelocity * (coordX[k_Q[k]] - coord0X);
		doubflo VeloXpur=VeloX;
		doubflo VeloYpur=VeloY;
		VeloX-=ux;
		VeloY-=uy;
		////////////////////////////////////////////////////////////////////////////////
		doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[dirE   *sizeQ];
		q_dirW   = &QQ[dirW   *sizeQ];
		q_dirN   = &QQ[dirN   *sizeQ];
		q_dirS   = &QQ[dirS   *sizeQ];
		q_dirT   = &QQ[dirT   *sizeQ];
		q_dirB   = &QQ[dirB   *sizeQ];
		q_dirNE  = &QQ[dirNE  *sizeQ];
		q_dirSW  = &QQ[dirSW  *sizeQ];
		q_dirSE  = &QQ[dirSE  *sizeQ];
		q_dirNW  = &QQ[dirNW  *sizeQ];
		q_dirTE  = &QQ[dirTE  *sizeQ];
		q_dirBW  = &QQ[dirBW  *sizeQ];
		q_dirBE  = &QQ[dirBE  *sizeQ];
		q_dirTW  = &QQ[dirTW  *sizeQ];
		q_dirTN  = &QQ[dirTN  *sizeQ];
		q_dirBS  = &QQ[dirBS  *sizeQ];
		q_dirBN  = &QQ[dirBN  *sizeQ];
		q_dirTS  = &QQ[dirTS  *sizeQ];
		q_dirTNE = &QQ[dirTNE *sizeQ];
		q_dirTSW = &QQ[dirTSW *sizeQ];
		q_dirTSE = &QQ[dirTSE *sizeQ];
		q_dirTNW = &QQ[dirTNW *sizeQ];
		q_dirBNE = &QQ[dirBNE *sizeQ];
		q_dirBSW = &QQ[dirBSW *sizeQ];
		q_dirBSE = &QQ[dirBSE *sizeQ];
		q_dirBNW = &QQ[dirBNW *sizeQ];
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
		//doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		//	f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		//f_W    = (D.f[dirE   ])[ke   ];
		//f_E    = (D.f[dirW   ])[kw   ];
		//f_S    = (D.f[dirN   ])[kn   ];
		//f_N    = (D.f[dirS   ])[ks   ];
		//f_B    = (D.f[dirT   ])[kt   ];
		//f_T    = (D.f[dirB   ])[kb   ];
		//f_SW   = (D.f[dirNE  ])[kne  ];
		//f_NE   = (D.f[dirSW  ])[ksw  ];
		//f_NW   = (D.f[dirSE  ])[kse  ];
		//f_SE   = (D.f[dirNW  ])[knw  ];
		//f_BW   = (D.f[dirTE  ])[kte  ];
		//f_TE   = (D.f[dirBW  ])[kbw  ];
		//f_TW   = (D.f[dirBE  ])[kbe  ];
		//f_BE   = (D.f[dirTW  ])[ktw  ];
		//f_BS   = (D.f[dirTN  ])[ktn  ];
		//f_TN   = (D.f[dirBS  ])[kbs  ];
		//f_TS   = (D.f[dirBN  ])[kbn  ];
		//f_BN   = (D.f[dirTS  ])[kts  ];
		//f_BSW  = (D.f[dirTNE ])[ktne ];
		//f_BNE  = (D.f[dirTSW ])[ktsw ];
		//f_BNW  = (D.f[dirTSE ])[ktse ];
		//f_BSE  = (D.f[dirTNW ])[ktnw ];
		//f_TSW  = (D.f[dirBNE ])[kbne ];
		//f_TNE  = (D.f[dirBSW ])[kbsw ];
		//f_TNW  = (D.f[dirBSE ])[kbse ];
		//f_TSE  = (D.f[dirBNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////
		doubflo /*vx1, vx2,*/ vx3, drho, feq, q, cu_sq;
		//drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
		//	f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
		//	f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

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
		if (evenOrOdd==false)
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
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Test
		//(D.f[dirZERO])[k]=c1o10;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//ToDo anders Klammern

		//q = q_dirE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        )-cu_sq); 
		//	(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
		//	//(D.f[dirW])[kw]=zero;
		//}

		//q = q_dirW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
		//	(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
		//	//(D.f[dirE])[ke]=zero;
		//}

		//q = q_dirN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
		//	(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
		//	//(D.f[dirS])[ks]=zero;
		//}

		//q = q_dirS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
		//	(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
		//	//(D.f[dirN])[kn]=zero;
		//}

		//q = q_dirT[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3)-cu_sq); 
		//	(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
		//	//(D.f[dirB])[kb]=one;
		//}

		//q = q_dirB[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
		//	(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
		//	//(D.f[dirT])[kt]=zero;
		//}

		//q = q_dirNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
		//	(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
		//	//(D.f[dirSW])[ksw]=zero;
		//}

		//q = q_dirSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
		//	(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
		//	//(D.f[dirNE])[kne]=zero;
		//}

		//q = q_dirSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
		//	(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
		//	//(D.f[dirNW])[knw]=zero;
		//}

		//q = q_dirNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
		//	(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
		//	//(D.f[dirSE])[kse]=zero;
		//}

		//q = q_dirTE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
		//	(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
		//	//(D.f[dirBW])[kbw]=zero;
		//}

		//q = q_dirBW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
		//	(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
		//	//(D.f[dirTE])[kte]=zero;
		//}

		//q = q_dirBE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
		//	(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
		//	//(D.f[dirTW])[ktw]=zero;
		//}

		//q = q_dirTW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
		//	(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
		//	//(D.f[dirBE])[kbe]=zero;
		//}

		//q = q_dirTN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
		//	(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBS])[kbs]=zero;
		//}

		//q = q_dirBS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
		//	(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTN])[ktn]=zero;
		//}

		//q = q_dirBN[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
		//	(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTS])[kts]=zero;
		//}

		//q = q_dirTS[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
		//	(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBN])[kbn]=zero;
		//}

		//q = q_dirTNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
		//	(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBSW])[kbsw]=zero;
		//}

		//q = q_dirBSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
		//	(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTNE])[ktne]=zero;
		//}

		//q = q_dirBNE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
		//	(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTSW])[ktsw]=zero;
		//}

		//q = q_dirTSW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
		//	(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBNE])[kbne]=zero;
		//}

		//q = q_dirTSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
		//	(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBNW])[kbnw]=zero;
		//}

		//q = q_dirBNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
		//	(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTSE])[ktse]=zero;
		//}

		//q = q_dirBSE[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
		//	(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
		//	//(D.f[dirTNW])[ktnw]=zero;
		//}

		//q = q_dirTNW[k];
		//if (q>=zero && q<=one)
		//{
		//	feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
		//	(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
		//	//(D.f[dirBSE])[kbse]=zero;
		//}

		///////// equilibrium BC
		cu_sq=c3o2*(VeloX*VeloX +VeloY*VeloY);
		VeloXpur*=-one;
		VeloYpur*=-one;
		vx3=zero;
		drho=zero;
		q = q_dirE[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*( VeloXpur        )+c9over2*( VeloX        )*( VeloX        )-cu_sq); 
			(D.f[dirW])[kw]=feq;
			//(D.f[dirW])[kw]=zero;
		}

		q = q_dirW[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*(-VeloXpur        )+c9over2*(-VeloX        )*(-VeloX        )-cu_sq); 
			(D.f[dirE])[ke]=feq;
			//(D.f[dirE])[ke]=zero;
		}

		q = q_dirN[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*(    VeloYpur     )+c9over2*(     VeloY    )*(     VeloY    )-cu_sq); 
			(D.f[dirS])[ks]=feq;
			//(D.f[dirS])[ks]=zero;
		}

		q = q_dirS[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*(   -VeloYpur     )+c9over2*(    -VeloY    )*(    -VeloY    )-cu_sq); 
			(D.f[dirN])[kn]=feq;
			//(D.f[dirN])[kn]=zero;
		}

		q = q_dirT[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
			(D.f[dirB])[kb]=feq;
			//(D.f[dirB])[kb]=one;
		}

		q = q_dirB[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
			(D.f[dirT])[kt]=feq;
			//(D.f[dirT])[kt]=zero;
		}

		q = q_dirNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*( VeloXpur+VeloYpur    )+c9over2*( VeloX+VeloY    )*( VeloX+VeloY    )-cu_sq); 
			(D.f[dirSW])[ksw]=feq;
			//(D.f[dirSW])[ksw]=zero;
		}

		q = q_dirSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(-VeloXpur-VeloYpur    )+c9over2*(-VeloX-VeloY    )*(-VeloX-VeloY    )-cu_sq); 
			(D.f[dirNE])[kne]=feq;
			//(D.f[dirNE])[kne]=zero;
		}

		q = q_dirSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*( VeloXpur-VeloYpur    )+c9over2*( VeloX-VeloY    )*( VeloX-VeloY    )-cu_sq); 
			(D.f[dirNW])[knw]=feq;
			//(D.f[dirNW])[knw]=zero;
		}

		q = q_dirNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(-VeloXpur+VeloYpur    )+c9over2*(-VeloX+VeloY    )*(-VeloX+VeloY    )-cu_sq); 
			(D.f[dirSE])[kse]=feq;
			//(D.f[dirSE])[kse]=zero;
		}

		q = q_dirTE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*( VeloXpur    +vx3)+c9over2*( VeloX    +vx3)*( VeloX    +vx3)-cu_sq); 
			(D.f[dirBW])[kbw]=feq;
			//(D.f[dirBW])[kbw]=zero;
		}

		q = q_dirBW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(-VeloXpur    -vx3)+c9over2*(-VeloX    -vx3)*(-VeloX    -vx3)-cu_sq); 
			(D.f[dirTE])[kte]=feq;
			//(D.f[dirTE])[kte]=zero;
		}

		q = q_dirBE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*( VeloXpur    -vx3)+c9over2*( VeloX    -vx3)*( VeloX    -vx3)-cu_sq); 
			(D.f[dirTW])[ktw]=feq;
			//(D.f[dirTW])[ktw]=zero;
		}

		q = q_dirTW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(-VeloXpur    +vx3)+c9over2*(-VeloX    +vx3)*(-VeloX    +vx3)-cu_sq); 
			(D.f[dirBE])[kbe]=feq;
			//(D.f[dirBE])[kbe]=zero;
		}

		q = q_dirTN[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(     VeloYpur+vx3)+c9over2*(     VeloY+vx3)*(     VeloY+vx3)-cu_sq); 
			(D.f[dirBS])[kbs]=feq;
			//(D.f[dirBS])[kbs]=zero;
		}

		q = q_dirBS[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(    -VeloYpur-vx3)+c9over2*(    -VeloY-vx3)*(    -VeloY-vx3)-cu_sq); 
			(D.f[dirTN])[ktn]=feq;
			//(D.f[dirTN])[ktn]=zero;
		}

		q = q_dirBN[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(     VeloYpur-vx3)+c9over2*(     VeloY-vx3)*(     VeloY-vx3)-cu_sq); 
			(D.f[dirTS])[kts]=feq;
			//(D.f[dirTS])[kts]=zero;
		}

		q = q_dirTS[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+three*(    -VeloYpur+vx3)+c9over2*(    -VeloY+vx3)*(    -VeloY+vx3)-cu_sq); 
			(D.f[dirBN])[kbn]=feq;
			//(D.f[dirBN])[kbn]=zero;
		}

		q = q_dirTNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*( VeloXpur+VeloYpur+vx3)+c9over2*( VeloX+VeloY+vx3)*( VeloX+VeloY+vx3)-cu_sq); 
			(D.f[dirBSW])[kbsw]=feq;
			//(D.f[dirBSW])[kbsw]=zero;
		}

		q = q_dirBSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*(-VeloXpur-VeloYpur-vx3)+c9over2*(-VeloX-VeloY-vx3)*(-VeloX-VeloY-vx3)-cu_sq); 
			(D.f[dirTNE])[ktne]=feq;
			//(D.f[dirTNE])[ktne]=zero;
		}

		q = q_dirBNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*( VeloXpur+VeloYpur-vx3)+c9over2*( VeloX+VeloY-vx3)*( VeloX+VeloY-vx3)-cu_sq); 
			(D.f[dirTSW])[ktsw]=feq;
			//(D.f[dirTSW])[ktsw]=zero;
		}

		q = q_dirTSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*(-VeloXpur-VeloYpur+vx3)+c9over2*(-VeloX-VeloY+vx3)*(-VeloX-VeloY+vx3)-cu_sq); 
			(D.f[dirBNE])[kbne]=feq;
			//(D.f[dirBNE])[kbne]=zero;
		}

		q = q_dirTSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*( VeloXpur-VeloYpur+vx3)+c9over2*( VeloX-VeloY+vx3)*( VeloX-VeloY+vx3)-cu_sq); 
			(D.f[dirBNW])[kbnw]=feq;
			//(D.f[dirBNW])[kbnw]=zero;
		}

		q = q_dirBNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*(-VeloXpur+VeloYpur-vx3)+c9over2*(-VeloX+VeloY-vx3)*(-VeloX+VeloY-vx3)-cu_sq); 
			(D.f[dirTSE])[ktse]=feq;
			//(D.f[dirTSE])[ktse]=zero;
		}

		q = q_dirBSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*( VeloXpur-VeloYpur-vx3)+c9over2*( VeloX-VeloY-vx3)*( VeloX-VeloY-vx3)-cu_sq); 
			(D.f[dirTNW])[ktnw]=feq;
			//(D.f[dirTNW])[ktnw]=zero;
		}

		q = q_dirTNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+three*(-VeloXpur+VeloYpur+vx3)+c9over2*(-VeloX+VeloY+vx3)*(-VeloX+VeloY+vx3)-cu_sq); 
			(D.f[dirBSE])[kbse]=feq;
			//(D.f[dirBSE])[kbse]=zero;
		}
	
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDeviceComp27(int inx,
											int iny,
											doubflo* vx,
											doubflo* vy,
											doubflo* vz,
											doubflo* DD, 
											int* k_Q, 
											doubflo* QQ,
											unsigned int sizeQ,
											int kQ, 
											doubflo om1, 
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat, 
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

   if(k<kQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      doubflo VeloX = vx[k];
      doubflo VeloY = vy[k];
      doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (one + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (one + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (one + drho); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (one + drho);

      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ) /** (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        ) * (one + drho)-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ) /** (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ) /** (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ) /** (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3) * (one + drho)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     )/* * (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3) * (one + drho)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ) /** (one + drho)*/)/(one+q);// - c2over27 * drho;
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    ) * (one + drho)-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3) * (one + drho)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3) * (one + drho)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over54 * drho;
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (one + drho)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (one + drho)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ) /** (one + drho)*/)/(one+q);// - c1over216 * drho;
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevice27(int inx,
                                        int iny,
                                        doubflo* vx,
                                        doubflo* vy,
                                        doubflo* vz,
                                        doubflo* DD, 
                                        int* k_Q, 
                                        doubflo* QQ,
                                        unsigned int sizeQ,
                                        int kQ, 
                                        doubflo om1, 
                                        unsigned int* neighborX,
                                        unsigned int* neighborY,
                                        unsigned int* neighborZ,
                                        unsigned int size_Mat, 
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

   if(k<kQ)
   {
      ////////////////////////////////////////////////////////////////////////////////
      doubflo VeloX = vx[k];
      doubflo VeloY = vy[k];
      doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
      ////////////////////////////////////////////////////////////////////////////////
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W); 
         

      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[dirZERO])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //ToDo anders Klammern

      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*( vx1        )*/+c9over2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
         //(D.f[dirW])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(-vx1        )*/+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
         //(D.f[dirE])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(    vx2     )*/+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
         //(D.f[dirS])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(   -vx2     )*/+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
         //(D.f[dirN])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(         vx3)*/+c9over2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
         //(D.f[dirB])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         feq=c2over27* (drho/*+three*(        -vx3)*/+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
         //(D.f[dirT])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1+vx2    )*/+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
         //(D.f[dirSW])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1-vx2    )*/+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
         //(D.f[dirNE])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1-vx2    )*/+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
         //(D.f[dirNW])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1+vx2    )*/+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
         //(D.f[dirSE])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    +vx3)*/+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
         //(D.f[dirBW])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    -vx3)*/+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
         //(D.f[dirTE])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*( vx1    -vx3)*/+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
         //(D.f[dirTW])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(-vx1    +vx3)*/+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
         //(D.f[dirBE])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2+vx3)*/+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
         //(D.f[dirBS])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2-vx3)*/+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
         //(D.f[dirTN])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(     vx2-vx3)*/+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
         //(D.f[dirTS])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         feq=c1over54* (drho/*+three*(    -vx2+vx3)*/+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
         //(D.f[dirBN])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2+vx3)*/+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSW])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2-vx3)*/+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNE])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1+vx2-vx3)*/+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSW])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1-vx2+vx3)*/+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNE])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2+vx3)*/+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
         //(D.f[dirBNW])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2-vx3)*/+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
         //(D.f[dirTSE])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*( vx1-vx2-vx3)*/+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
         //(D.f[dirTNW])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         feq=c1over216*(drho/*+three*(-vx1+vx2+vx3)*/+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
         //(D.f[dirBSE])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void PropellerBC(unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       doubflo* rho,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       int* k_Q, 
									   unsigned int size_Prop,
                                       unsigned int size_Mat,
                                       unsigned int* bcMatD,
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

   if(k<size_Prop)
   {
    ////////////////////////////////////////////////////////////////////////////////
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
		unsigned int KQK = k_Q[k];
		unsigned int BC  = bcMatD[KQK];
		if( (BC != GEO_SOLID) && (BC != GEO_VOID))
		{		
		//////////////////////////////////////////////////////////////////////////
        doubflo  vx1 = ux[k];
        doubflo  vx2 = uy[k];
        doubflo  vx3 = uz[k];
        //doubflo  vx1 = -c1o100;
        //doubflo  vx2 = zero;
        //doubflo  vx3 = zero;
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
		doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
		f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

		f_ZERO= (D.f[dirZERO])[kzero];
		f_E   = (D.f[dirE   ])[ke   ];
		f_W   = (D.f[dirW   ])[kw   ];
		f_N   = (D.f[dirN   ])[kn   ];
		f_S   = (D.f[dirS   ])[ks   ];
		f_T   = (D.f[dirT   ])[kt   ];
		f_B   = (D.f[dirB   ])[kb   ];
		f_NE  = (D.f[dirNE  ])[kne  ];
		f_SW  = (D.f[dirSW  ])[ksw  ];
		f_SE  = (D.f[dirSE  ])[kse  ];
		f_NW  = (D.f[dirNW  ])[knw  ];
		f_TE  = (D.f[dirTE  ])[kte  ];
		f_BW  = (D.f[dirBW  ])[kbw  ];
		f_BE  = (D.f[dirBE  ])[kbe  ];
		f_TW  = (D.f[dirTW  ])[ktw  ];
		f_TN  = (D.f[dirTN  ])[ktn  ];
		f_BS  = (D.f[dirBS  ])[kbs  ];
		f_BN  = (D.f[dirBN  ])[kbn  ];
		f_TS  = (D.f[dirTS  ])[kts  ];
		f_TNE = (D.f[dirTNE ])[ktne ];
		f_BSW = (D.f[dirBSW ])[kbsw ];
		f_BNE = (D.f[dirBNE ])[kbne ];
		f_TSW = (D.f[dirTSW ])[ktsw ];
		f_TSE = (D.f[dirTSE ])[ktse ];
		f_BNW = (D.f[dirBNW ])[kbnw ];
		f_BSE = (D.f[dirBSE ])[kbse ];
		f_TNW = (D.f[dirTNW ])[ktnw ];
		//f_W    = (D.f[dirE   ])[ke   ];
		//f_E    = (D.f[dirW   ])[kw   ];
		//f_S    = (D.f[dirN   ])[kn   ];
		//f_N    = (D.f[dirS   ])[ks   ];
		//f_B    = (D.f[dirT   ])[kt   ];
		//f_T    = (D.f[dirB   ])[kb   ];
		//f_SW   = (D.f[dirNE  ])[kne  ];
		//f_NE   = (D.f[dirSW  ])[ksw  ];
		//f_NW   = (D.f[dirSE  ])[kse  ];
		//f_SE   = (D.f[dirNW  ])[knw  ];
		//f_BW   = (D.f[dirTE  ])[kte  ];
		//f_TE   = (D.f[dirBW  ])[kbw  ];
		//f_TW   = (D.f[dirBE  ])[kbe  ];
		//f_BE   = (D.f[dirTW  ])[ktw  ];
		//f_BS   = (D.f[dirTN  ])[ktn  ];
		//f_TN   = (D.f[dirBS  ])[kbs  ];
		//f_TS   = (D.f[dirBN  ])[kbn  ];
		//f_BN   = (D.f[dirTS  ])[kts  ];
		//f_BSW  = (D.f[dirTNE ])[ktne ];
		//f_TNE  = (D.f[dirBSW ])[kbsw ];
		//f_TSW  = (D.f[dirBNE ])[kbne ];
		//f_BNE  = (D.f[dirTSW ])[ktsw ];
		//f_BNW  = (D.f[dirTSE ])[ktse ];
		//f_TSE  = (D.f[dirBNW ])[kbnw ];
		//f_TNW  = (D.f[dirBSE ])[kbse ];
		//f_BSE  = (D.f[dirTNW ])[ktnw ];
		//////////////////////////////////////////////////////////////////////////////////
		doubflo vxo1, vxo2, vxo3, drho;
		drho   =  /*zero;*/f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
				  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
				  f_T + f_B + f_N + f_S + f_E + f_W + f_ZERO; 

		vxo1   =   (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
					(f_E - f_W) )/ (one + drho); 
        

		vxo2   =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
					((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
					(f_N - f_S) )/ (one + drho); 

		vxo3   =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
		 			(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
					(f_T - f_B) )/ (one + drho); 

		doubflo cusq=c3o2*(vxo1*vxo1+vxo2*vxo2+vxo3*vxo3);
		//vx1 = vx1 * two - vxo1;
		//vx2 = vx2 * two - vxo2;
		//vx3 = vx3 * two - vxo3;
		doubflo cusq2=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

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
         f_ZERO = f_ZERO + ((one+drho) * (-  c8over27* (-cusq)																   +   c8over27* (-cusq2)));
         f_E    = f_E    + ((one+drho) * (-  c2over27* (three*( vxo1          )+c9over2*( vxo1          )*( vxo1          )-cusq) +   c2over27* (three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq2)));
         f_W    = f_W    + ((one+drho) * (-  c2over27* (three*(-vxo1          )+c9over2*(-vxo1          )*(-vxo1          )-cusq) +   c2over27* (three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq2)));
         f_N    = f_N    + ((one+drho) * (-  c2over27* (three*(      vxo2     )+c9over2*(      vxo2     )*(      vxo2     )-cusq) +   c2over27* (three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq2)));
         f_S    = f_S    + ((one+drho) * (-  c2over27* (three*(     -vxo2     )+c9over2*(     -vxo2     )*(     -vxo2     )-cusq) +   c2over27* (three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq2)));
         f_T    = f_T    + ((one+drho) * (-  c2over27* (three*(           vxo3)+c9over2*(           vxo3)*(           vxo3)-cusq) +   c2over27* (three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq2)));
         f_B    = f_B    + ((one+drho) * (-  c2over27* (three*(          -vxo3)+c9over2*(          -vxo3)*(          -vxo3)-cusq) +   c2over27* (three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq2)));
         f_NE   = f_NE   + ((one+drho) * (-  c1over54* (three*( vxo1+vxo2     )+c9over2*( vxo1+vxo2     )*( vxo1+vxo2     )-cusq) +   c1over54* (three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq2)));
         f_SW   = f_SW   + ((one+drho) * (-  c1over54* (three*(-vxo1-vxo2     )+c9over2*(-vxo1-vxo2     )*(-vxo1-vxo2     )-cusq) +   c1over54* (three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq2)));
         f_SE   = f_SE   + ((one+drho) * (-  c1over54* (three*( vxo1-vxo2     )+c9over2*( vxo1-vxo2     )*( vxo1-vxo2     )-cusq) +   c1over54* (three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq2)));
         f_NW   = f_NW   + ((one+drho) * (-  c1over54* (three*(-vxo1+vxo2     )+c9over2*(-vxo1+vxo2     )*(-vxo1+vxo2     )-cusq) +   c1over54* (three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq2)));
         f_TE   = f_TE   + ((one+drho) * (-  c1over54* (three*( vxo1     +vxo3)+c9over2*( vxo1     +vxo3)*( vxo1     +vxo3)-cusq) +   c1over54* (three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq2)));
         f_BW   = f_BW   + ((one+drho) * (-  c1over54* (three*(-vxo1     -vxo3)+c9over2*(-vxo1     -vxo3)*(-vxo1     -vxo3)-cusq) +   c1over54* (three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq2)));
         f_BE   = f_BE   + ((one+drho) * (-  c1over54* (three*( vxo1     -vxo3)+c9over2*( vxo1     -vxo3)*( vxo1     -vxo3)-cusq) +   c1over54* (three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq2)));
         f_TW   = f_TW   + ((one+drho) * (-  c1over54* (three*(-vxo1     +vxo3)+c9over2*(-vxo1     +vxo3)*(-vxo1     +vxo3)-cusq) +   c1over54* (three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq2)));
         f_TN   = f_TN   + ((one+drho) * (-  c1over54* (three*(      vxo2+vxo3)+c9over2*(      vxo2+vxo3)*(      vxo2+vxo3)-cusq) +   c1over54* (three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq2)));
         f_BS   = f_BS   + ((one+drho) * (-  c1over54* (three*(     -vxo2-vxo3)+c9over2*(     -vxo2-vxo3)*(     -vxo2-vxo3)-cusq) +   c1over54* (three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq2)));
         f_BN   = f_BN   + ((one+drho) * (-  c1over54* (three*(      vxo2-vxo3)+c9over2*(      vxo2-vxo3)*(      vxo2-vxo3)-cusq) +   c1over54* (three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq2)));
         f_TS   = f_TS   + ((one+drho) * (-  c1over54* (three*(     -vxo2+vxo3)+c9over2*(     -vxo2+vxo3)*(     -vxo2+vxo3)-cusq) +   c1over54* (three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq2)));
         f_TNE  = f_TNE  + ((one+drho) * (-  c1over216*(three*( vxo1+vxo2+vxo3)+c9over2*( vxo1+vxo2+vxo3)*( vxo1+vxo2+vxo3)-cusq) +   c1over216*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq2)));
         f_BSW  = f_BSW  + ((one+drho) * (-  c1over216*(three*(-vxo1-vxo2-vxo3)+c9over2*(-vxo1-vxo2-vxo3)*(-vxo1-vxo2-vxo3)-cusq) +   c1over216*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq2)));
         f_BNE  = f_BNE  + ((one+drho) * (-  c1over216*(three*( vxo1+vxo2-vxo3)+c9over2*( vxo1+vxo2-vxo3)*( vxo1+vxo2-vxo3)-cusq) +   c1over216*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq2)));
         f_TSW  = f_TSW  + ((one+drho) * (-  c1over216*(three*(-vxo1-vxo2+vxo3)+c9over2*(-vxo1-vxo2+vxo3)*(-vxo1-vxo2+vxo3)-cusq) +   c1over216*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq2)));
         f_TSE  = f_TSE  + ((one+drho) * (-  c1over216*(three*( vxo1-vxo2+vxo3)+c9over2*( vxo1-vxo2+vxo3)*( vxo1-vxo2+vxo3)-cusq) +   c1over216*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq2)));
         f_BNW  = f_BNW  + ((one+drho) * (-  c1over216*(three*(-vxo1+vxo2-vxo3)+c9over2*(-vxo1+vxo2-vxo3)*(-vxo1+vxo2-vxo3)-cusq) +   c1over216*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq2)));
         f_BSE  = f_BSE  + ((one+drho) * (-  c1over216*(three*( vxo1-vxo2-vxo3)+c9over2*( vxo1-vxo2-vxo3)*( vxo1-vxo2-vxo3)-cusq) +   c1over216*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq2)));
         f_TNW  = f_TNW  + ((one+drho) * (-  c1over216*(three*(-vxo1+vxo2+vxo3)+c9over2*(-vxo1+vxo2+vxo3)*(-vxo1+vxo2+vxo3)-cusq) +   c1over216*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq2)));

		(D.f[dirZERO])[kzero] =  f_ZERO;
        (D.f[dirE   ])[ke   ] =  f_E   ;	// f_W   ;//    	
        (D.f[dirW   ])[kw   ] =  f_W   ;	// f_E   ;//    	
        (D.f[dirN   ])[kn   ] =  f_N   ;	// f_S   ;//    	
        (D.f[dirS   ])[ks   ] =  f_S   ;	// f_N   ;//    	
        (D.f[dirT   ])[kt   ] =  f_T   ;	// f_B   ;//    	
        (D.f[dirB   ])[kb   ] =  f_B   ;	// f_T   ;//    	
        (D.f[dirNE  ])[kne  ] =  f_NE  ;	// f_SW  ;//    	
        (D.f[dirSW  ])[ksw  ] =  f_SW  ;	// f_NE  ;//    	
        (D.f[dirSE  ])[kse  ] =  f_SE  ;	// f_NW  ;//    	
        (D.f[dirNW  ])[knw  ] =  f_NW  ;	// f_SE  ;//    	
        (D.f[dirTE  ])[kte  ] =  f_TE  ;	// f_BW  ;//    	
        (D.f[dirBW  ])[kbw  ] =  f_BW  ;	// f_TE  ;//    	
        (D.f[dirBE  ])[kbe  ] =  f_BE  ;	// f_TW  ;//    	
        (D.f[dirTW  ])[ktw  ] =  f_TW  ;	// f_BE  ;//    	
        (D.f[dirTN  ])[ktn  ] =  f_TN  ;	// f_BS  ;//    	
        (D.f[dirBS  ])[kbs  ] =  f_BS  ;	// f_TN  ;//    	
        (D.f[dirBN  ])[kbn  ] =  f_BN  ;	// f_TS  ;//    	
        (D.f[dirTS  ])[kts  ] =  f_TS  ;	// f_BN  ;//    	
        (D.f[dirTNE ])[ktne ] =  f_TNE ;	// f_BSW ;//    	
        (D.f[dirBSW ])[kbsw ] =  f_BSW ;	// f_BNE ;//    	
        (D.f[dirBNE ])[kbne ] =  f_BNE ;	// f_BNW ;//    	
        (D.f[dirTSW ])[ktsw ] =  f_TSW ;	// f_BSE ;//    	
        (D.f[dirTSE ])[ktse ] =  f_TSE ;	// f_TSW ;//    	
        (D.f[dirBNW ])[kbnw ] =  f_BNW ;	// f_TNE ;//    	
        (D.f[dirBSE ])[kbse ] =  f_BSE ;	// f_TNW ;//    	
        (D.f[dirTNW ])[ktnw ] =  f_TNW ;	// f_TSE ;//    	

		//////////////////////////////////////////////////////////////////////////
        ////(D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[dirE   ])[ke   ] =   three*c2over27* ( vx1        );		//six
        //(D.f[dirW   ])[kw   ] =   three*c2over27* (-vx1        );		//six
        //(D.f[dirN   ])[kn   ] =   three*c2over27* (     vx2    );		//six
        //(D.f[dirS   ])[ks   ] =   three*c2over27* (    -vx2    );		//six
        //(D.f[dirT   ])[kt   ] =   three*c2over27* (         vx3);		//six
        //(D.f[dirB   ])[kb   ] =   three*c2over27* (        -vx3);		//six
        //(D.f[dirNE  ])[kne  ] =   three*c1over54* ( vx1+vx2    );		//six
        //(D.f[dirSW  ])[ksw  ] =   three*c1over54* (-vx1-vx2    );		//six
        //(D.f[dirSE  ])[kse  ] =   three*c1over54* ( vx1-vx2    );		//six
        //(D.f[dirNW  ])[knw  ] =   three*c1over54* (-vx1+vx2    );		//six
        //(D.f[dirTE  ])[kte  ] =   three*c1over54* ( vx1    +vx3);		//six
        //(D.f[dirBW  ])[kbw  ] =   three*c1over54* (-vx1    -vx3);		//six
        //(D.f[dirBE  ])[kbe  ] =   three*c1over54* ( vx1    -vx3);		//six
        //(D.f[dirTW  ])[ktw  ] =   three*c1over54* (-vx1    +vx3);		//six
        //(D.f[dirTN  ])[ktn  ] =   three*c1over54* (     vx2+vx3);		//six
        //(D.f[dirBS  ])[kbs  ] =   three*c1over54* (    -vx2-vx3);		//six
        //(D.f[dirBN  ])[kbn  ] =   three*c1over54* (     vx2-vx3);		//six
        //(D.f[dirTS  ])[kts  ] =   three*c1over54* (    -vx2+vx3);		//six
        //(D.f[dirTNE ])[ktne ] =   three*c1over216*( vx1+vx2+vx3);		//six
        //(D.f[dirBSW ])[kbsw ] =   three*c1over216*(-vx1-vx2-vx3);		//six
        //(D.f[dirBNE ])[kbne ] =   three*c1over216*( vx1+vx2-vx3);		//six
        //(D.f[dirTSW ])[ktsw ] =   three*c1over216*(-vx1-vx2+vx3);		//six
        //(D.f[dirTSE ])[ktse ] =   three*c1over216*( vx1-vx2+vx3);		//six
        //(D.f[dirBNW ])[kbnw ] =   three*c1over216*(-vx1+vx2-vx3);		//six
        //(D.f[dirBSE ])[kbse ] =   three*c1over216*( vx1-vx2-vx3);		//six
        //(D.f[dirTNW ])[ktnw ] =   three*c1over216*(-vx1+vx2+vx3);		//six
        //(D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq);
        //(D.f[dirE   ])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
        //(D.f[dirW   ])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
        //(D.f[dirN   ])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
        //(D.f[dirS   ])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
        //(D.f[dirT   ])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
        //(D.f[dirB   ])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
        //(D.f[dirNE  ])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
        //(D.f[dirSW  ])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
        //(D.f[dirSE  ])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
        //(D.f[dirNW  ])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
        //(D.f[dirTE  ])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
        //(D.f[dirBW  ])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
        //(D.f[dirBE  ])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
        //(D.f[dirTW  ])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
        //(D.f[dirTN  ])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
        //(D.f[dirBS  ])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
        //(D.f[dirBN  ])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
        //(D.f[dirTS  ])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
        //(D.f[dirTNE ])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
        //(D.f[dirBSW ])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
        //(D.f[dirBNE ])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
        //(D.f[dirTSW ])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
        //(D.f[dirTSE ])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
        //(D.f[dirBNW ])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
        //(D.f[dirBSE ])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
        //(D.f[dirTNW ])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
		}
    }
}
//////////////////////////////////////////////////////////////////////////


