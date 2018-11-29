//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

//////////////////////////////////////////////////////////////////////////
extern "C" __global__ void scaleFC_comp_D3Q27F3_2018(real* DC,
													 real* DF,
													 real* G6,
													 unsigned int* neighborCX,
													 unsigned int* neighborCY,
													 unsigned int* neighborCZ,
													 unsigned int* neighborFX,
													 unsigned int* neighborFY,
													 unsigned int* neighborFZ,
													 unsigned int size_MatC, 
													 unsigned int size_MatF, 
													 bool evenOrOdd,
													 unsigned int* posC, 
													 unsigned int* posFSWB, 
													 unsigned int kFC, 
													 real omCoarse, 
													 real omFine, 
													 real nu, 
													 unsigned int nxC, 
													 unsigned int nyC, 
													 unsigned int nxF, 
													 unsigned int nyF,
													 OffFC offFC)
{
   real 
	   *fP00source, *fM00source, *f0P0source, *f0M0source, *f00Psource, *f00Msource, *fPP0source, *fMM0source, *fPM0source,
	   *fMP0source, *fP0Psource, *fM0Msource, *fP0Msource, *fM0Psource, *f0PPsource, *f0MMsource, *f0PMsource, *f0MPsource,
	   *f000source, *fMMMsource, *fMMPsource, *fMPPsource, *fMPMsource, *fPPMsource, *fPPPsource, *fPMPsource, *fPMMsource;


   fP00source = &DF[dirE   *size_MatF];
   fM00source = &DF[dirW   *size_MatF];
   f0P0source = &DF[dirN   *size_MatF];
   f0M0source = &DF[dirS   *size_MatF];
   f00Psource = &DF[dirT   *size_MatF];
   f00Msource = &DF[dirB   *size_MatF];
   fPP0source = &DF[dirNE  *size_MatF];
   fMM0source = &DF[dirSW  *size_MatF];
   fPM0source = &DF[dirSE  *size_MatF];
   fMP0source = &DF[dirNW  *size_MatF];
   fP0Psource = &DF[dirTE  *size_MatF];
   fM0Msource = &DF[dirBW  *size_MatF];
   fP0Msource = &DF[dirBE  *size_MatF];
   fM0Psource = &DF[dirTW  *size_MatF];
   f0PPsource = &DF[dirTN  *size_MatF];
   f0MMsource = &DF[dirBS  *size_MatF];
   f0PMsource = &DF[dirBN  *size_MatF];
   f0MPsource = &DF[dirTS  *size_MatF];
   f000source = &DF[dirZERO*size_MatF];
   fMMMsource = &DF[dirBSW *size_MatF];
   fMMPsource = &DF[dirTSW *size_MatF];
   fMPPsource = &DF[dirTNW *size_MatF];
   fMPMsource = &DF[dirBNW *size_MatF];
   fPPMsource = &DF[dirBNE *size_MatF];
   fPPPsource = &DF[dirTNE *size_MatF];
   fPMPsource = &DF[dirTSE *size_MatF];
   fPMMsource = &DF[dirBSE *size_MatF];

   real
	   *fP00dest, *fM00dest, *f0P0dest, *f0M0dest, *f00Pdest, *f00Mdest, *fPP0dest, *fMM0dest, *fPM0dest,
	   *fMP0dest, *fP0Pdest, *fM0Mdest, *fP0Mdest, *fM0Pdest, *f0PPdest, *f0MMdest, *f0PMdest, *f0MPdest,
	   *f000dest, *fMMMdest, *fMMPdest, *fMPPdest, *fMPMdest, *fPPMdest, *fPPPdest, *fPMPdest, *fPMMdest;

   if (evenOrOdd==true)
   {
	   fP00dest = &DC[dirE   *size_MatC];
	   fM00dest = &DC[dirW   *size_MatC];
	   f0P0dest = &DC[dirN   *size_MatC];
	   f0M0dest = &DC[dirS   *size_MatC];
	   f00Pdest = &DC[dirT   *size_MatC];
	   f00Mdest = &DC[dirB   *size_MatC];
	   fPP0dest = &DC[dirNE  *size_MatC];
	   fMM0dest = &DC[dirSW  *size_MatC];
	   fPM0dest = &DC[dirSE  *size_MatC];
	   fMP0dest = &DC[dirNW  *size_MatC];
	   fP0Pdest = &DC[dirTE  *size_MatC];
	   fM0Mdest = &DC[dirBW  *size_MatC];
	   fP0Mdest = &DC[dirBE  *size_MatC];
	   fM0Pdest = &DC[dirTW  *size_MatC];
	   f0PPdest = &DC[dirTN  *size_MatC];
	   f0MMdest = &DC[dirBS  *size_MatC];
	   f0PMdest = &DC[dirBN  *size_MatC];
	   f0MPdest = &DC[dirTS  *size_MatC];
	   f000dest = &DC[dirZERO*size_MatC];
	   fMMMdest = &DC[dirBSW *size_MatC];
	   fMMPdest = &DC[dirTSW *size_MatC];
	   fMPPdest = &DC[dirTNW *size_MatC];
	   fMPMdest = &DC[dirBNW *size_MatC];
	   fPPMdest = &DC[dirBNE *size_MatC];
	   fPPPdest = &DC[dirTNE *size_MatC];
	   fPMPdest = &DC[dirTSE *size_MatC];
	   fPMMdest = &DC[dirBSE *size_MatC];
   } 
   else
   {
	   fP00dest = &DC[dirW   *size_MatC];
	   fM00dest = &DC[dirE   *size_MatC];
	   f0P0dest = &DC[dirS   *size_MatC];
	   f0M0dest = &DC[dirN   *size_MatC];
	   f00Pdest = &DC[dirB   *size_MatC];
	   f00Mdest = &DC[dirT   *size_MatC];
	   fPP0dest = &DC[dirSW  *size_MatC];
	   fMM0dest = &DC[dirNE  *size_MatC];
	   fPM0dest = &DC[dirNW  *size_MatC];
	   fMP0dest = &DC[dirSE  *size_MatC];
	   fP0Pdest = &DC[dirBW  *size_MatC];
	   fM0Mdest = &DC[dirTE  *size_MatC];
	   fP0Mdest = &DC[dirTW  *size_MatC];
	   fM0Pdest = &DC[dirBE  *size_MatC];
	   f0PPdest = &DC[dirBS  *size_MatC];
	   f0MMdest = &DC[dirTN  *size_MatC];
	   f0PMdest = &DC[dirTS  *size_MatC];
	   f0MPdest = &DC[dirBN  *size_MatC];
	   f000dest = &DC[dirZERO*size_MatC];
	   fMMMdest = &DC[dirTNE *size_MatC];
	   fMMPdest = &DC[dirBNE *size_MatC];
	   fMPPdest = &DC[dirBSE *size_MatC];
	   fMPMdest = &DC[dirTSE *size_MatC];
	   fPPMdest = &DC[dirTSW *size_MatC];
	   fPPPdest = &DC[dirBSW *size_MatC];
	   fPMPdest = &DC[dirBNW *size_MatC];
	   fPMMdest = &DC[dirTNW *size_MatC];
   }

   Distributions6 G;
   if (evenOrOdd == true)
   {
	   G.g[dirE] = &G6[dirE   *size_MatC];
	   G.g[dirW] = &G6[dirW   *size_MatC];
	   G.g[dirN] = &G6[dirN   *size_MatC];
	   G.g[dirS] = &G6[dirS   *size_MatC];
	   G.g[dirT] = &G6[dirT   *size_MatC];
	   G.g[dirB] = &G6[dirB   *size_MatC];
   }
   else
   {
	   G.g[dirW] = &G6[dirE   *size_MatC];
	   G.g[dirE] = &G6[dirW   *size_MatC];
	   G.g[dirS] = &G6[dirN   *size_MatC];
	   G.g[dirN] = &G6[dirS   *size_MatC];
	   G.g[dirB] = &G6[dirT   *size_MatC];
	   G.g[dirT] = &G6[dirB   *size_MatC];
   }

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  ix = threadIdx.x;  // Globaler x-Index 
   const unsigned  iy = blockIdx.x;   // Globaler y-Index 
   const unsigned  iz = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*iz + iy) + ix;
   //////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   real eps_new = two;
   real omegaS = omFine;//-omFine;
   real o  = omCoarse;//-omCoarse;
   //real op = one;
   //real cu_sq;

   real xoff,    yoff,    zoff;
   real xoff_sq, yoff_sq, zoff_sq;

   real        vvx, vvy, vvz, vx2, vy2, vz2, drho;
   real        press;//,drho,vx1,vx2,vx3;
   real        /*pressMMM,*/drhoMMM,vx1MMM,vx2MMM,vx3MMM;
   real        /*pressMMP,*/drhoMMP,vx1MMP,vx2MMP,vx3MMP;
   real        /*pressMPP,*/drhoMPP,vx1MPP,vx2MPP,vx3MPP;
   real        /*pressMPM,*/drhoMPM,vx1MPM,vx2MPM,vx3MPM;
   real        /*pressPPM,*/drhoPPM,vx1PPM,vx2PPM,vx3PPM;
   real        /*pressPPP,*/drhoPPP,vx1PPP,vx2PPP,vx3PPP;
   real        /*pressPMP,*/drhoPMP,vx1PMP,vx2PMP,vx3PMP;
   real        /*pressPMM,*/drhoPMM,vx1PMM,vx2PMM,vx3PMM;
   real        fP00,fM00,f0P0,f0M0,f00P,f00M,fPP0,fMM0,fPM0,fMP0,fP0P,fM0M,fP0M,fM0P,f0PP,f0MM,f0PM,f0MP,f000,fPPP, fMMP, fPMP, fMPP, fPPM, fMMM, fPMM, fMPM;
   real        kxyFromfcNEQMMM, kyzFromfcNEQMMM, kxzFromfcNEQMMM, kxxMyyFromfcNEQMMM, kxxMzzFromfcNEQMMM, kyyMzzFromfcNEQMMM;
   real        kxyFromfcNEQMMP, kyzFromfcNEQMMP, kxzFromfcNEQMMP, kxxMyyFromfcNEQMMP, kxxMzzFromfcNEQMMP, kyyMzzFromfcNEQMMP;
   real        kxyFromfcNEQMPP, kyzFromfcNEQMPP, kxzFromfcNEQMPP, kxxMyyFromfcNEQMPP, kxxMzzFromfcNEQMPP, kyyMzzFromfcNEQMPP;
   real        kxyFromfcNEQMPM, kyzFromfcNEQMPM, kxzFromfcNEQMPM, kxxMyyFromfcNEQMPM, kxxMzzFromfcNEQMPM, kyyMzzFromfcNEQMPM;
   real        kxyFromfcNEQPPM, kyzFromfcNEQPPM, kxzFromfcNEQPPM, kxxMyyFromfcNEQPPM, kxxMzzFromfcNEQPPM, kyyMzzFromfcNEQPPM;
   real        kxyFromfcNEQPPP, kyzFromfcNEQPPP, kxzFromfcNEQPPP, kxxMyyFromfcNEQPPP, kxxMzzFromfcNEQPPP, kyyMzzFromfcNEQPPP;
   real        kxyFromfcNEQPMP, kyzFromfcNEQPMP, kxzFromfcNEQPMP, kxxMyyFromfcNEQPMP, kxxMzzFromfcNEQPMP, kyyMzzFromfcNEQPMP;
   real        kxyFromfcNEQPMM, kyzFromfcNEQPMM, kxzFromfcNEQPMM, kxxMyyFromfcNEQPMM, kxxMzzFromfcNEQPMM, kyyMzzFromfcNEQPMM;
   real        a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz, b0, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz, c0, cx, cy, cz, cxx, cyy, czz, cxy, cxz, cyz/*, axyz, bxyz, cxyz*/;
   real        d0, dx, dy, dz, dxy, dxz, dyz/*, dxyz*/;

   if(k<kFC)
   {
      //////////////////////////////////////////////////////////////////////////
      xoff = offFC.xOffFC[k];
      yoff = offFC.yOffFC[k];
      zoff = offFC.zOffFC[k];      
      xoff_sq = xoff * xoff;
      yoff_sq = yoff * yoff;
      zoff_sq = zoff * zoff;
      //////////////////////////////////////////////////////////////////////////
      //SWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      unsigned int k000base = posFSWB[k];
      unsigned int kM00base = neighborFX[k000base];
      unsigned int k0M0base = neighborFY[k000base];
      unsigned int k00Mbase = neighborFZ[k000base];
      unsigned int kMM0base = neighborFY[kM00base];
      unsigned int kM0Mbase = neighborFZ[kM00base];
      unsigned int k0MMbase = neighborFZ[k0M0base];
      unsigned int kMMMbase = neighborFZ[kMM0base];
      //////////////////////////////////////////////////////////////////////////
      //index 
      unsigned int k000 = k000base;
      unsigned int kM00 = kM00base;   
      unsigned int k0M0 = k0M0base;   
      unsigned int k00M = k00Mbase;   
      unsigned int kMM0 = kMM0base;  
      unsigned int kM0M = kM0Mbase;  
      unsigned int k0MM = k0MMbase;  
      unsigned int kMMM = kMMMbase; 
      ////////////////////////////////////////////////////////////////////////////////
      fP00 = fP00source[k000];
      fM00 = fM00source[kM00];
      f0P0 = f0P0source[k000];
      f0M0 = f0M0source[k0M0];
      f00P = f00Psource[k000];
      f00M = f00Msource[k00M];
      fPP0 = fPP0source[k000];
      fMM0 = fMM0source[kMM0];
      fPM0 = fPM0source[k0M0];
      fMP0 = fMP0source[kM00];
      fP0P = fP0Psource[k000];
      fM0M = fM0Msource[kM0M];
      fP0M = fP0Msource[k00M];
      fM0P = fM0Psource[kM00];
      f0PP = f0PPsource[k000];
      f0MM = f0MMsource[k0MM];
      f0PM = f0PMsource[k00M];
      f0MP = f0MPsource[k0M0];
      f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
      fPMP = fPMPsource[k0M0];
      fPMM = fPMMsource[k0MM];

      drhoMMM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MMM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMMM);
	  vx2MMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMMM);
	  vx3MMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMMM);

	  kxyFromfcNEQMMM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMMM) - ((vx1MMM*vx2MMM)));
	  kyzFromfcNEQMMM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMMM) - ((vx2MMM*vx3MMM)));
	  kxzFromfcNEQMMM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMMM) - ((vx1MMM*vx3MMM)));
	  kxxMyyFromfcNEQMMM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMMM) - ((vx1MMM*vx1MMM - vx2MMM*vx2MMM)));
	  kxxMzzFromfcNEQMMM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMMM) - ((vx1MMM*vx1MMM - vx3MMM*vx3MMM)));
	  kyyMzzFromfcNEQMMM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMMM) - ((vx2MMM*vx2MMM - vx3MMM*vx3MMM)));

      //////////////////////////////////////////////////////////////////////////
      //SWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k00M;
      kM00 = kM0M;   
      k0M0 = k0MM;   
      k00M = neighborFZ[k00M];   
      kMM0 = kMMM;  
      kM0M = neighborFZ[kM0M];  
      k0MM = neighborFZ[k0MM];  
      kMMM = neighborFZ[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMMP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MMP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMMP);
	  vx2MMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMMP);
	  vx3MMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMMP);

	  kxyFromfcNEQMMP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMMP) - ((vx1MMP*vx2MMP)));
	  kyzFromfcNEQMMP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMMP) - ((vx2MMP*vx3MMP)));
	  kxzFromfcNEQMMP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMMP) - ((vx1MMP*vx3MMP)));
	  kxxMyyFromfcNEQMMP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMMP) - ((vx1MMP*vx1MMP - vx2MMP*vx2MMP)));
	  kxxMzzFromfcNEQMMP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMMP) - ((vx1MMP*vx1MMP - vx3MMP*vx3MMP)));
	  kyyMzzFromfcNEQMMP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMMP) - ((vx2MMP*vx2MMP - vx3MMP*vx3MMP)));

      //////////////////////////////////////////////////////////////////////////
      //SET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = kM00;
      kM00 = neighborFX[kM00];   
      k0M0 = kMM0;   
      k00M = kM0M;   
      kMM0 = neighborFX[kMM0];  
      kM0M = neighborFX[kM0M];  
      k0MM = kMMM;  
      kMMM = neighborFX[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPMP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PMP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPMP);
	  vx2PMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPMP);
	  vx3PMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPMP);

	  kxyFromfcNEQPMP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPMP) - ((vx1PMP*vx2PMP)));
	  kyzFromfcNEQPMP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPMP) - ((vx2PMP*vx3PMP)));
	  kxzFromfcNEQPMP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPMP) - ((vx1PMP*vx3PMP)));
	  kxxMyyFromfcNEQPMP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPMP) - ((vx1PMP*vx1PMP - vx2PMP*vx2PMP)));
	  kxxMzzFromfcNEQPMP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPMP) - ((vx1PMP*vx1PMP - vx3PMP*vx3PMP)));
	  kyyMzzFromfcNEQPMP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPMP) - ((vx2PMP*vx2PMP - vx3PMP*vx3PMP)));

      //////////////////////////////////////////////////////////////////////////
      //SEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k00M = k000;   
      kM0M = kM00;  
      k0MM = k0M0;  
      kMMM = kMM0; 
      k000 = kM00base;
      kM00 = neighborFX[kM00base];   
      k0M0 = kMM0base;   
      kMM0 = neighborFX[kMM0base];  
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPMM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PMM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPMM);
	  vx2PMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPMM);
	  vx3PMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPMM);

	  kxyFromfcNEQPMM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPMM) - ((vx1PMM*vx2PMM)));
	  kyzFromfcNEQPMM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPMM) - ((vx2PMM*vx3PMM)));
	  kxzFromfcNEQPMM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPMM) - ((vx1PMM*vx3PMM)));
	  kxxMyyFromfcNEQPMM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPMM) - ((vx1PMM*vx1PMM - vx2PMM*vx2PMM)));
	  kxxMzzFromfcNEQPMM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPMM) - ((vx1PMM*vx1PMM - vx3PMM*vx3PMM)));
	  kyyMzzFromfcNEQPMM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPMM) - ((vx2PMM*vx2PMM - vx3PMM*vx3PMM)));

      //////////////////////////////////////////////////////////////////////////
      //NWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      k000base = k0M0base;
      kM00base = kMM0base;
      k0M0base = neighborFY[k0M0base];
      k00Mbase = k0MMbase;
      kMM0base = neighborFY[kMM0base];
      kM0Mbase = kMMMbase;
      k0MMbase = neighborFY[k0MMbase];
      kMMMbase = neighborFY[kMMMbase];
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k000base;
      kM00 = kM00base;   
      k0M0 = k0M0base;   
      k00M = k00Mbase;   
      kMM0 = kMM0base;  
      kM0M = kM0Mbase;  
      k0MM = k0MMbase;  
      kMMM = kMMMbase; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMPM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MPM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMPM);
	  vx2MPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMPM);
	  vx3MPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMPM);

	  kxyFromfcNEQMPM    = -three*omegaS*   ((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMPM) - ((vx1MPM*vx2MPM)));
	  kyzFromfcNEQMPM    = -three*omegaS*   ((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMPM) - ((vx2MPM*vx3MPM)));
	  kxzFromfcNEQMPM    = -three*omegaS*   ((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMPM) - ((vx1MPM*vx3MPM)));
	  kxxMyyFromfcNEQMPM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMPM) - ((vx1MPM*vx1MPM - vx2MPM*vx2MPM)));
	  kxxMzzFromfcNEQMPM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMPM) - ((vx1MPM*vx1MPM - vx3MPM*vx3MPM)));
	  kyyMzzFromfcNEQMPM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMPM) - ((vx2MPM*vx2MPM - vx3MPM*vx3MPM)));

	  //////////////////////////////////////////////////////////////////////////
      //NWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k00M;
      kM00 = kM0M;   
      k0M0 = k0MM;   
      k00M = neighborFZ[k00M];   
      kMM0 = kMMM;  
      kM0M = neighborFZ[kM0M];  
      k0MM = neighborFZ[k0MM];  
      kMMM = neighborFZ[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMPP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MPP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMPP);
	  vx2MPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMPP);
	  vx3MPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMPP);

	  kxyFromfcNEQMPP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMPP) - ((vx1MPP*vx2MPP)));
	  kyzFromfcNEQMPP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMPP) - ((vx2MPP*vx3MPP)));
	  kxzFromfcNEQMPP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMPP) - ((vx1MPP*vx3MPP)));
	  kxxMyyFromfcNEQMPP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMPP) - ((vx1MPP*vx1MPP - vx2MPP*vx2MPP)));
	  kxxMzzFromfcNEQMPP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMPP) - ((vx1MPP*vx1MPP - vx3MPP*vx3MPP)));
	  kyyMzzFromfcNEQMPP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMPP) - ((vx2MPP*vx2MPP - vx3MPP*vx3MPP)));

      //////////////////////////////////////////////////////////////////////////
      //NET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = kM00;
      kM00 = neighborFX[kM00];   
      k0M0 = kMM0;   
      k00M = kM0M;   
      kMM0 = neighborFX[kMM0];  
      kM0M = neighborFX[kM0M];  
      k0MM = kMMM;  
      kMMM = neighborFX[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPPP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PPP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPPP);
	  vx2PPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPPP);
	  vx3PPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPPP);

	  kxyFromfcNEQPPP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPPP) - ((vx1PPP*vx2PPP)));
	  kyzFromfcNEQPPP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPPP) - ((vx2PPP*vx3PPP)));
	  kxzFromfcNEQPPP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPPP) - ((vx1PPP*vx3PPP)));
	  kxxMyyFromfcNEQPPP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPPP) - ((vx1PPP*vx1PPP - vx2PPP*vx2PPP)));
	  kxxMzzFromfcNEQPPP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPPP) - ((vx1PPP*vx1PPP - vx3PPP*vx3PPP)));
	  kyyMzzFromfcNEQPPP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPPP) - ((vx2PPP*vx2PPP - vx3PPP*vx3PPP)));

      //////////////////////////////////////////////////////////////////////////
      //NEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k00M = k000;   
      kM0M = kM00;  
      k0MM = k0M0;  
      kMMM = kMM0; 
      k000 = kM00base;
      kM00 = neighborFX[kM00base];   
      k0M0 = kMM0base;   
      kMM0 = neighborFX[kMM0base];  
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPPM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PPM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPPM);
	  vx2PPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPPM);
	  vx3PPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPPM);

	  kxyFromfcNEQPPM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPPM) - ((vx1PPM*vx2PPM)));
	  kyzFromfcNEQPPM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPPM) - ((vx2PPM*vx3PPM)));
	  kxzFromfcNEQPPM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPPM) - ((vx1PPM*vx3PPM)));
	  kxxMyyFromfcNEQPPM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPPM) - ((vx1PPM*vx1PPM - vx2PPM*vx2PPM)));
	  kxxMzzFromfcNEQPPM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPPM) - ((vx1PPM*vx1PPM - vx3PPM*vx3PPM)));
	  kyyMzzFromfcNEQPPM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPPM) - ((vx2PPM*vx2PPM - vx3PPM*vx3PPM)));

      //////////////////////////////////////////////////////////////////////////
      //3
      //////////////////////////////////////////////////////////////////////////
	  a0  = c1o8*(((vx1PPM + vx1MMP) + (vx1MPM + vx1PMP)) + ((vx1PMM + vx1MPP) + (vx1MMM + vx1PPP)));
	  ax  = c1o4*(((vx1PPM - vx1MMP) + (vx1PPP - vx1MMM)) + ((vx1PMP - vx1MPM) + (vx1PMM - vx1MPP)));
	  ay  = c1o4*(((vx1PPM - vx1MMP) + (vx1PPP - vx1MMM)) + ((vx1MPM - vx1PMP) + (vx1MPP - vx1PMM)));
	  az  = c1o4*(((vx1MMP - vx1PPM) + (vx1PPP - vx1MMM)) + ((vx1PMP - vx1MPM) + (vx1MPP - vx1PMM)));
	  axy = c1o2*(((vx1PPM - vx1PMP) + (vx1MMM - vx1MPP)) + ((vx1MMP - vx1MPM) + (vx1PPP - vx1PMM)));
	  axz = c1o2*(((vx1PMP - vx1PPM) + (vx1MMM - vx1MPP)) + ((vx1MPM - vx1MMP) + (vx1PPP - vx1PMM)));
	  ayz = c1o2*(((vx1PPP - vx1MPM) + (vx1PMM - vx1MMP)) + ((vx1MPP - vx1PPM) + (vx1MMM - vx1PMP)));
	  //axyz=		  ((vx1MMP - vx1PPM) + (vx1PPP - vx1MMM)) + ((vx1MPM - vx1PMP) + (vx1PMM - vx1MPP)) ;
	  b0  = c1o8*(((vx2PPM + vx2MMP) + (vx2MPM + vx2PMP)) + ((vx2PMM + vx2MPP) + (vx2MMM + vx2PPP)));
	  bx  = c1o4*(((vx2PPM - vx2MMP) + (vx2PPP - vx2MMM)) + ((vx2PMP - vx2MPM) + (vx2PMM - vx2MPP)));
	  by  = c1o4*(((vx2PPM - vx2MMP) + (vx2PPP - vx2MMM)) + ((vx2MPM - vx2PMP) + (vx2MPP - vx2PMM)));
	  bz  = c1o4*(((vx2MMP - vx2PPM) + (vx2PPP - vx2MMM)) + ((vx2PMP - vx2MPM) + (vx2MPP - vx2PMM)));
	  bxy = c1o2*(((vx2PPM - vx2MPP) + (vx2MMM - vx2PMP)) + ((vx2MMP - vx2PMM) + (vx2PPP - vx2MPM)));
	  bxz = c1o2*(((vx2MMM - vx2PPM) + (vx2PMP - vx2MPP)) + ((vx2MPM - vx2PMM) + (vx2PPP - vx2MMP)));
	  byz = c1o2*(((vx2MPP - vx2PPM) + (vx2MMM - vx2PMP)) + ((vx2PMM - vx2MMP) + (vx2PPP - vx2MPM)));
	  //bxyz=		  ((vx2MMP - vx2PPM) + (vx2PPP - vx2MMM)) + ((vx2MPM - vx2PMP) + (vx2PMM - vx2MPP)) ;
	  c0  = c1o8*(((vx3PPM + vx3MMP) + (vx3MPM + vx3PMP)) + ((vx3PMM + vx3MPP) + (vx3MMM + vx3PPP)));
	  cx  = c1o4*(((vx3PPM - vx3MMP) + (vx3PPP - vx3MMM)) + ((vx3PMP - vx3MPM) + (vx3PMM - vx3MPP)));
	  cy  = c1o4*(((vx3PPM - vx3MMP) + (vx3PPP - vx3MMM)) + ((vx3MPM - vx3PMP) + (vx3MPP - vx3PMM)));
	  cz  = c1o4*(((vx3MMP - vx3PPM) + (vx3PPP - vx3MMM)) + ((vx3PMP - vx3MPM) + (vx3MPP - vx3PMM)));
	  cxy = c1o2*(((vx3PPM - vx3PMP) + (vx3MMM - vx3MPP)) + ((vx3MMP - vx3MPM) + (vx3PPP - vx3PMM)));
	  cxz = c1o2*(((vx3MMM - vx3PPM) + (vx3PMP - vx3MPP)) + ((vx3MPM - vx3PMM) + (vx3PPP - vx3MMP)));
	  cyz = c1o2*(((vx3MMM - vx3PPM) + (vx3MPP - vx3PMP)) + ((vx3PMM - vx3MPM) + (vx3PPP - vx3MMP)));
	  //cxyz=		  ((vx3MMP - vx3PPM) + (vx3PPP - vx3MMM)) + ((vx3MPM - vx3PMP) + (vx3PMM - vx3MPP)) ;

	  //bxx = zero;
	  //cxx = zero;
	  //ayy = zero;
	  //cyy = zero;
	  //azz = zero;
	  //bzz = zero;
	  //axx = zero;
	  //byy = zero;
	  //czz = zero;

	  bxx = c1o8*(((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPMM - kxyFromfcNEQMPP)) + ((kxyFromfcNEQPMP - kxyFromfcNEQMPM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP))) - c1o2*axy;
	  cxx = c1o8*(((kxzFromfcNEQPPP - kxzFromfcNEQMMM) + (kxzFromfcNEQPMM - kxzFromfcNEQMPP)) + ((kxzFromfcNEQPMP - kxzFromfcNEQMPM) + (kxzFromfcNEQPPM - kxzFromfcNEQMMP))) - c1o2*axz;

	  ayy = c1o8*(((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP)) - ((kxyFromfcNEQPMM - kxyFromfcNEQMPP) + (kxyFromfcNEQPMP - kxyFromfcNEQMPM))) - c1o2*bxy;
	  cyy = c1o8*(((kyzFromfcNEQPPP - kyzFromfcNEQMMM) + (kyzFromfcNEQPPM - kyzFromfcNEQMMP)) - ((kyzFromfcNEQPMM - kyzFromfcNEQMPP) + (kyzFromfcNEQPMP - kyzFromfcNEQMPM))) - c1o2*byz;

	  azz = c1o8*(((kxzFromfcNEQPPP - kxzFromfcNEQMMM) - (kxzFromfcNEQPMM - kxzFromfcNEQMPP)) + ((kxzFromfcNEQPMP - kxzFromfcNEQMPM) - (kxzFromfcNEQPPM - kxzFromfcNEQMMP))) - c1o2*cxz;
	  bzz = c1o8*(((kyzFromfcNEQPPP - kyzFromfcNEQMMM) - (kyzFromfcNEQPMM - kyzFromfcNEQMPP)) + ((kyzFromfcNEQPMP - kyzFromfcNEQMPM) - (kyzFromfcNEQPPM - kyzFromfcNEQMMP))) - c1o2*cyz;

	  axx = ( c1o16*(((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) + (kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP)) + ((kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM) + (kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP))) - c1o4*bxy)
		  + ( c1o16*(((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) + (kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP)) + ((kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM) + (kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP))) - c1o4*cxz);

	  byy = (-c1o16*(((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) - (kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP)) + ((kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP) - (kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM))) - c1o4*axy)
		  + ( c1o16*(((kyyMzzFromfcNEQPPP - kyyMzzFromfcNEQMMM) - (kyyMzzFromfcNEQPMM - kyyMzzFromfcNEQMPP)) + ((kyyMzzFromfcNEQPPM - kyyMzzFromfcNEQMMP) - (kyyMzzFromfcNEQPMP - kyyMzzFromfcNEQMPM))) - c1o4*cyz);

	  czz = (-c1o16*(((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) - (kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP)) - ((kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP) - (kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM))) - c1o4*axz)
		  + ( c1o16*(((kyyMzzFromfcNEQPPP - kyyMzzFromfcNEQMMM) - (kyyMzzFromfcNEQPMM - kyyMzzFromfcNEQMPP)) - ((kyyMzzFromfcNEQPPM - kyyMzzFromfcNEQMMP) - (kyyMzzFromfcNEQPMP - kyyMzzFromfcNEQMPM))) - c1o4*byz);

	  a0 -= c1o4*(axx + ayy + azz);
	  b0 -= c1o4*(bxx + byy + bzz);
	  c0 -= c1o4*(cxx + cyy + czz);

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real kxyAverage = zero;
	  real kyzAverage = zero;
	  real kxzAverage = zero;
	  real kxxMyyAverage = zero;
	  real kxxMzzAverage = zero;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////Press
	  //d0   = ( pressPPM + pressPPP + pressMPM + pressMPP + pressPMM + pressPMP + pressMMM + pressMMP) * c1o8;
	  //dx   = ( pressPPM + pressPPP - pressMPM - pressMPP + pressPMM + pressPMP - pressMMM - pressMMP) * c1o4;
	  //dy   = ( pressPPM + pressPPP + pressMPM + pressMPP - pressPMM - pressPMP - pressMMM - pressMMP) * c1o4;
	  //dz   = (-pressPPM + pressPPP - pressMPM + pressMPP - pressPMM + pressPMP - pressMMM + pressMMP) * c1o4;
	  //dxy  = ( pressPPM + pressPPP - pressMPM - pressMPP - pressPMM - pressPMP + pressMMM + pressMMP) * c1o2;
	  //dxz  = (-pressPPM + pressPPP + pressMPM - pressMPP - pressPMM + pressPMP + pressMMM - pressMMP) * c1o2;
	  //dyz  = (-pressPPM + pressPPP - pressMPM + pressMPP + pressPMM - pressPMP + pressMMM - pressMMP) * c1o2;
	  //dxyz =  -pressPPM + pressPPP + pressMPM - pressMPP + pressPMM - pressPMP - pressMMM + pressMMP;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //drho
	  d0   = ( ((drhoPPM + drhoMMP) + (drhoMPM + drhoPMP)) + ((drhoPMM + drhoMPP) + (drhoMMM + drhoPPP))) * c1o8;
	  dx   = ( ((drhoPPM - drhoMMP) + (drhoPMM - drhoMPP)) + ((drhoPMP - drhoMPM) + (drhoPPP - drhoMMM))) * c1o4;
	  dy   = ( ((drhoPPM - drhoMMP) + (drhoMPP - drhoPMM)) + ((drhoMPM - drhoPMP) + (drhoPPP - drhoMMM))) * c1o4;
	  dz   = ( ((drhoMMP - drhoPPM) + (drhoPPP - drhoMMM)) + ((drhoPMP - drhoMPM) + (drhoMPP - drhoPMM))) * c1o4;
	  dxy  = ( ((drhoPPM - drhoPMP) + (drhoPPP - drhoPMM)) + ((drhoMMP - drhoMPM) + (drhoMMM - drhoMPP))) * c1o2;
	  dxz  = ( ((drhoMMM - drhoPPM) + (drhoPPP - drhoMMP)) + ((drhoMPM - drhoPMM) + (drhoPMP - drhoMPP))) * c1o2;
	  dyz  = ( ((drhoMPP - drhoPPM) + (drhoPPP - drhoMPM)) + ((drhoPMM - drhoMMP) + (drhoMMM - drhoPMP))) * c1o2;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // Bernd das Brot 
	  //
      //
	  // x------x
	  // |      |
	  // |	 ---+--->X
	  // |		|  \
	  // x------x   \
	  //			off-vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      a0 = a0 + xoff * ax + yoff * ay + zoff * az + xoff_sq * axx + yoff_sq * ayy + zoff_sq * azz + xoff*yoff*axy + xoff*zoff*axz + yoff*zoff*ayz;
      ax = ax + two * xoff * axx + yoff * axy + zoff * axz;
      ay = ay + two * yoff * ayy + xoff * axy + zoff * ayz;
      az = az + two * zoff * azz + xoff * axz + yoff * ayz;
      b0 = b0 + xoff * bx + yoff * by + zoff * bz + xoff_sq * bxx + yoff_sq * byy + zoff_sq * bzz + xoff*yoff*bxy + xoff*zoff*bxz + yoff*zoff*byz;
      bx = bx + two * xoff * bxx + yoff * bxy + zoff * bxz;
      by = by + two * yoff * byy + xoff * bxy + zoff * byz;
      bz = bz + two * zoff * bzz + xoff * bxz + yoff * byz;
      c0 = c0 + xoff * cx + yoff * cy + zoff * cz + xoff_sq * cxx + yoff_sq * cyy + zoff_sq * czz + xoff*yoff*cxy + xoff*zoff*cxz + yoff*zoff*cyz;
      cx = cx + two * xoff * cxx + yoff * cxy + zoff * cxz;
      cy = cy + two * yoff * cyy + xoff * cxy + zoff * cyz;
      cz = cz + two * zoff * czz + xoff * cxz + yoff * cyz;
	  d0 = d0 + xoff * dx + yoff * dy + zoff * dz + xoff*yoff*dxy + xoff*zoff*dxz + yoff*zoff*dyz;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real mfcbb = zero;
	  real mfabb = zero;
	  real mfbcb = zero;
	  real mfbab = zero;
	  real mfbbc = zero;
	  real mfbba = zero;
	  real mfccb = zero;
	  real mfaab = zero;
	  real mfcab = zero;
	  real mfacb = zero;
	  real mfcbc = zero;
	  real mfaba = zero;
	  real mfcba = zero;
	  real mfabc = zero;
	  real mfbcc = zero;
	  real mfbaa = zero;
	  real mfbca = zero;
	  real mfbac = zero;
	  real mfbbb = zero;
	  real mfccc = zero;
	  real mfaac = zero;
	  real mfcac = zero;
	  real mfacc = zero;
	  real mfcca = zero;
	  real mfaaa = zero;
	  real mfcaa = zero;
	  real mfaca = zero;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real mgcbb = zero;
	  real mgabb = zero;
	  real mgbcb = zero;
	  real mgbab = zero;
	  real mgbbc = zero;
	  real mgbba = zero;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real m0, m1, m2, oMdrho;
	  real mxxPyyPzz, mxxMyy, mxxMzz, mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
	  real qudricLimit = c1o100;//ganz schlechte Idee -> muss global sein
	  real O3 = two - o;
	  real residu, residutmp;
	  residutmp = zero;///*-*/ c2o9 * (1./o - c1o2) * eps_new * eps_new;
	  real NeqOn = one;//zero;//one;   //.... one = on ..... zero = off 
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //
	  //Position C 0., 0., 0.
	  //
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //x = 0.;
	  //y = 0.;
	  //z = 0.;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real mxoff = -xoff;
	  //real myoff = -yoff;
	  //real mzoff = -zoff;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //press = pressPPP * (c1o8 - c1o4 * mxoff - c1o4 * myoff - c1o4 * mzoff) + 
			//  pressMPP * (c1o8 + c1o4 * mxoff - c1o4 * myoff - c1o4 * mzoff) + 
			//  pressPMP * (c1o8 - c1o4 * mxoff + c1o4 * myoff - c1o4 * mzoff) + 
			//  pressMMP * (c1o8 + c1o4 * mxoff + c1o4 * myoff - c1o4 * mzoff) + 
			//  pressPPM * (c1o8 - c1o4 * mxoff - c1o4 * myoff + c1o4 * mzoff) + 
			//  pressMPM * (c1o8 + c1o4 * mxoff - c1o4 * myoff + c1o4 * mzoff) + 
			//  pressPMM * (c1o8 - c1o4 * mxoff + c1o4 * myoff + c1o4 * mzoff) + 
			//  pressMMM * (c1o8 + c1o4 * mxoff + c1o4 * myoff + c1o4 * mzoff);
	  //drho  = drhoPPP * (c1o8 - c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) + 
			//  drhoMPP * (c1o8 + c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) + 
			//  drhoPMP * (c1o8 - c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) + 
			//  drhoMMP * (c1o8 + c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) + 
			//  drhoPPM * (c1o8 - c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) + 
			//  drhoMPM * (c1o8 + c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) + 
			//  drhoPMM * (c1o8 - c1o4 * xoff + c1o4 * yoff + c1o4 * zoff) + 
			//  drhoMMM * (c1o8 + c1o4 * xoff + c1o4 * yoff + c1o4 * zoff);
	  press = d0;
	  vvx   = a0;
	  vvy   = b0;
	  vvz   = c0;

	  //mfaaa = drho;
	  //mfaaa = press + (ax+by+cz)/three;  //  1/3 = 2/3*(1/op-1/2)
	  mfaaa = press; // if drho is interpolated directly

	  vx2 = vvx*vvx;
	  vy2 = vvy*vvy;
	  vz2 = vvz*vvz;
	  oMdrho = one;
	  //oMdrho = one - mfaaa;

	  //two
	  // linear combinations
	  real oP = o;// :(
	  mxxPyyPzz = mfaaa    -c2o3*(ax+by+cz)*eps_new/oP*(one+press); 
	  //mxxMyy    = -c2o3*(ax - by)*eps_new/o;
	  //mxxMzz    = -c2o3*(ax - cz)*eps_new/o;

	  //mfabb     = -c1o3 * (bz + cy)*eps_new/o;
	  //mfbab     = -c1o3 * (az + cx)*eps_new/o;
	  //mfbba     = -c1o3 * (ay + bx)*eps_new/o;
	  mxxMyy    = -c2o3*((ax - by)+kxxMyyAverage)*eps_new/o * (one + press);
	  mxxMzz    = -c2o3*((ax - cz)+kxxMzzAverage)*eps_new/o * (one + press);

	  mfabb     = -c1o3 * ((bz + cy)+kyzAverage)*eps_new/o * (one + press);
	  mfbab     = -c1o3 * ((az + cx)+kxzAverage)*eps_new/o * (one + press);
	  mfbba     = -c1o3 * ((ay + bx)+kxyAverage)*eps_new/o * (one + press);

	  
	  // linear combinations back
	  mfcaa = c1o3 * (       mxxMyy +       mxxMzz + mxxPyyPzz) * NeqOn;
	  mfaca = c1o3 * (-two * mxxMyy +       mxxMzz + mxxPyyPzz) * NeqOn;
	  mfaac = c1o3 * (       mxxMyy - two * mxxMzz + mxxPyyPzz) * NeqOn;

	  //3.
	  // linear combinations
	  //residu = residutmp * (ayz + bxz + cxy );
	  //mfbbb = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mfbbb = zero;

	  //residu = residutmp * (axy + two*bxx + two*bzz + cyz );
	  //residu = -(c1o9*(axy - 2*bxx - 2*bzz + cyz ));
	  //mxxyPyzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxyPyzz = zero;

	  //residu = residutmp * (axy + two*bxx - two*bzz - cyz );
	  //residu = c1o9*(axy - 2*bxx + 2*bzz - cyz );
	  //mxxyMyzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxyMyzz = zero;

	  //residu = residutmp * (axz + byz + two*cxx + two*cyy );
	  //residu = -(c1o9*(axz + byz - 2*cxx - 2*cyy ));
	  //mxxzPyyz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxzPyyz = zero;

	  //residu = residutmp * (axz - byz + two*cxx - two*cyy );
	  //residu = c1o9*(axz - byz - 2*cxx + 2*cyy );
	  //mxxzMyyz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxzMyyz = zero;

	  //residu = residutmp * (two*ayy + two*azz + bxy + cxz );
	  //residu = c1o9*(2*ayy + 2*azz - bxy - cxz );
	  //mxyyPxzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxyyPxzz = zero;

	  //residu = residutmp * (two*ayy - two*azz + bxy - cxz );
	  //residu = c1o9*(-2*ayy + 2*azz + bxy - cxz );
	  //mxyyMxzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxyyMxzz = zero;

	  ////////////////////////////////////////////////////////////////////////////////////
	  // D3Q27F 
	  mgcbb =  (vvx * axx + ax * ax) * (eps_new * eps_new) * (one + press);
	  mgabb = -(vvx * axx + ax * ax) * (eps_new * eps_new) * (one + press);
	  mgbcb =  (vvy * byy + by * by) * (eps_new * eps_new) * (one + press);
	  mgbab = -(vvy * byy + by * by) * (eps_new * eps_new) * (one + press);
	  mgbbc =  (vvz * czz + cz * cz) * (eps_new * eps_new) * (one + press);
	  mgbba = -(vvz * czz + cz * cz) * (eps_new * eps_new) * (one + press);
	  //mgcbb = zero;
	  //mgabb = zero;
	  //mgbcb = zero;
	  //mgbab = zero;
	  //mgbbc = zero;
	  //mgbba = zero;
	  ////////////////////////////////////////////////////////////////////////////////////

	  // linear combinations back
	  mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
	  mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
	  mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
	  mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
	  mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
	  mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

	  //4.
	  mfacc = mfaaa*c1o9; 
	  mfcac = mfacc; 
	  mfcca = mfacc; 
	  //5.

	  //6.
	  mfccc = mfaaa*c1o27;
	  ////////////////////////////////////////////////////////////////////////////////////
	  //back
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // Z - Dir
	  m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfaac        - two * mfaab *  vvz         +  mfaaa                * (one - vz2)              - one * oMdrho * vz2; 
	  m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfaaa = m0;
	  mfaab = m1;
	  mfaac = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfabc        - two * mfabb *  vvz         + mfaba * (one - vz2); 
	  m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
	  mfaba = m0;
	  mfabb = m1;
	  mfabc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfacc        - two * mfacb *  vvz         +  mfaca                  * (one - vz2)              - c1o3 * oMdrho * vz2; 
	  m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfaca = m0;
	  mfacb = m1;
	  mfacc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbac        - two * mfbab *  vvz         + mfbaa * (one - vz2); 
	  m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
	  mfbaa = m0;
	  mfbab = m1;
	  mfbac = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbbc        - two * mfbbb *  vvz         + mfbba * (one - vz2); 
	  m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
	  mfbba = m0;
	  mfbbb = m1;
	  mfbbc = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbcc        - two * mfbcb *  vvz         + mfbca * (one - vz2); 
	  m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
	  mfbca = m0;
	  mfbcb = m1;
	  mfbcc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfcac        - two * mfcab *  vvz         +  mfcaa                  * (one - vz2)              - c1o3 * oMdrho * vz2; 
	  m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfcaa = m0;
	  mfcab = m1;
	  mfcac = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfcbc        - two * mfcbb *  vvz         + mfcba * (one - vz2); 
	  m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
	  mfcba = m0;
	  mfcbb = m1;
	  mfcbc = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfccc        - two * mfccb *  vvz         +  mfcca                  * (one - vz2)              - c1o9 * oMdrho * vz2; 
	  m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfcca = m0;
	  mfccb = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // Y - Dir
	  m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfaca        - two * mfaba *  vvy         +  mfaaa                  * (one - vy2)              - c1o6 * oMdrho * vy2; 
	  m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaaa = m0;
	  mfaba = m1;
	  mfaca = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfacb        - two * mfabb *  vvy         +  mfaab                  * (one - vy2)              - c2o3 * oMdrho * vy2; 
	  m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaab = m0;
	  mfabb = m1;
	  mfacb = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfacc        - two * mfabc *  vvy         +  mfaac                  * (one - vy2)              - c1o6 * oMdrho * vy2; 
	  m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaac = m0;
	  mfabc = m1;
	  mfacc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbca        - two * mfbba *  vvy         + mfbaa * (one - vy2); 
	  m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
	  mfbaa = m0;
	  mfbba = m1;
	  mfbca = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbcb        - two * mfbbb *  vvy         + mfbab * (one - vy2); 
	  m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
	  mfbab = m0;
	  mfbbb = m1;
	  mfbcb = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbcc        - two * mfbbc *  vvy         + mfbac * (one - vy2); 
	  m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
	  mfbac = m0;
	  mfbbc = m1;
	  mfbcc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfcca        - two * mfcba *  vvy         +  mfcaa                   * (one - vy2)              - c1o18 * oMdrho * vy2; 
	  m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcaa = m0;
	  mfcba = m1;
	  mfcca = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfccb        - two * mfcbb *  vvy         +  mfcab                  * (one - vy2)              - c2o9 * oMdrho * vy2; 
	  m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcab = m0;
	  mfcbb = m1;
	  mfccb = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfccc        - two * mfcbc *  vvy         +  mfcac                   * (one - vy2)              - c1o18 * oMdrho * vy2; 
	  m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcac = m0;
	  mfcbc = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // X - Dir
	  m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcaa        - two * mfbaa *  vvx         +  mfaaa                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaaa = m0;
	  mfbaa = m1;
	  mfcaa = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcba        - two * mfbba *  vvx         +  mfaba                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaba = m0;
	  mfbba = m1;
	  mfcba = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcca        - two * mfbca *  vvx         +  mfaca                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaca = m0;
	  mfbca = m1;
	  mfcca = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcab        - two * mfbab *  vvx         +  mfaab                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaab = m0;
	  mfbab = m1;
	  mfcab = m2;
	  ///////////b////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcbb        - two * mfbbb *  vvx         +  mfabb                  * (one - vx2)              - c4o9 * oMdrho * vx2; 
	  m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfabb = m0;
	  mfbbb = m1;
	  mfcbb = m2;
	  ///////////b////////////////////////////////////////////////////////////////////////
	  m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfccb        - two * mfbcb *  vvx         +  mfacb                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfacb = m0;
	  mfbcb = m1;
	  mfccb = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcac        - two * mfbac *  vvx         +  mfaac                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaac = m0;
	  mfbac = m1;
	  mfcac = m2;
	  ///////////c////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcbc        - two * mfbbc *  vvx         +  mfabc                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfabc = m0;
	  mfbbc = m1;
	  mfcbc = m2;
	  ///////////c////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfccc        - two * mfbcc *  vvx         +  mfacc                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfacc = m0;
	  mfbcc = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////
	  //index 0
	  k000 = posC[k];
	  kM00 = neighborCX[k000];
	  k0M0 = neighborCY[k000];
	  k00M = neighborCZ[k000];
	  kMM0 = neighborCY[kM00];
	  kM0M = neighborCZ[kM00];
	  k0MM = neighborCZ[k0M0];
	  kMMM = neighborCZ[kMM0];
	  ////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////
	  (G.g[dirE])[k000] = mgcbb;
	  (G.g[dirW])[kM00] = mgabb;
	  (G.g[dirN])[k000] = mgbcb;
	  (G.g[dirS])[k0M0] = mgbab;
	  (G.g[dirT])[k000] = mgbbc;
	  (G.g[dirB])[k00M] = mgbba;
	  ////////////////////////////////////////////////////////////////////////////////////
	  fP00dest[k000] = mfcbb;                                                                 
	  fM00dest[kM00] = mfabb;                                                               
	  f0P0dest[k000] = mfbcb;
	  f0M0dest[k0M0] = mfbab;
	  f00Pdest[k000] = mfbbc;
	  f00Mdest[k00M] = mfbba;
	  fPP0dest[k000] = mfccb;
	  fMM0dest[kMM0] = mfaab;
	  fPM0dest[k0M0] = mfcab;
	  fMP0dest[kM00] = mfacb;
	  fP0Pdest[k000] = mfcbc;
	  fM0Mdest[kM0M] = mfaba;
	  fP0Mdest[k00M] = mfcba;
	  fM0Pdest[kM00] = mfabc;
	  f0PPdest[k000] = mfbcc;
	  f0MMdest[k0MM] = mfbaa;
	  f0PMdest[k00M] = mfbca;
	  f0MPdest[k0M0] = mfbac;
	  f000dest[k000] = mfbbb;
	  fMMMdest[kMMM] = mfaaa;
	  fMMPdest[kMM0] = mfaac;
	  fMPPdest[kM00] = mfacc;
	  fMPMdest[kM0M] = mfaca;
	  fPPMdest[k00M] = mfcca;
	  fPPPdest[k000] = mfccc;
	  fPMPdest[k0M0] = mfcac;
	  fPMMdest[k0MM] = mfcaa;
	  ////////////////////////////////////////////////////////////////////////////////////
   }
}
//////////////////////////////////////////////////////////////////////////






















































//////////////////////////////////////////////////////////////////////////
extern "C" __global__ void scaleFC_comp_D3Q27F3( real* DC,
												 real* DF,
												 real* G6,
												 unsigned int* neighborCX,
												 unsigned int* neighborCY,
												 unsigned int* neighborCZ,
												 unsigned int* neighborFX,
												 unsigned int* neighborFY,
												 unsigned int* neighborFZ,
												 unsigned int size_MatC, 
												 unsigned int size_MatF, 
												 bool evenOrOdd,
												 unsigned int* posC, 
												 unsigned int* posFSWB, 
												 unsigned int kFC, 
												 real omCoarse, 
												 real omFine, 
												 real nu, 
												 unsigned int nxC, 
												 unsigned int nyC, 
												 unsigned int nxF, 
												 unsigned int nyF,
												 OffFC offFC)
{
   real 
	   *fP00source, *fM00source, *f0P0source, *f0M0source, *f00Psource, *f00Msource, *fPP0source, *fMM0source, *fPM0source,
	   *fMP0source, *fP0Psource, *fM0Msource, *fP0Msource, *fM0Psource, *f0PPsource, *f0MMsource, *f0PMsource, *f0MPsource,
	   *f000source, *fMMMsource, *fMMPsource, *fMPPsource, *fMPMsource, *fPPMsource, *fPPPsource, *fPMPsource, *fPMMsource;


   fP00source = &DF[dirE   *size_MatF];
   fM00source = &DF[dirW   *size_MatF];
   f0P0source = &DF[dirN   *size_MatF];
   f0M0source = &DF[dirS   *size_MatF];
   f00Psource = &DF[dirT   *size_MatF];
   f00Msource = &DF[dirB   *size_MatF];
   fPP0source = &DF[dirNE  *size_MatF];
   fMM0source = &DF[dirSW  *size_MatF];
   fPM0source = &DF[dirSE  *size_MatF];
   fMP0source = &DF[dirNW  *size_MatF];
   fP0Psource = &DF[dirTE  *size_MatF];
   fM0Msource = &DF[dirBW  *size_MatF];
   fP0Msource = &DF[dirBE  *size_MatF];
   fM0Psource = &DF[dirTW  *size_MatF];
   f0PPsource = &DF[dirTN  *size_MatF];
   f0MMsource = &DF[dirBS  *size_MatF];
   f0PMsource = &DF[dirBN  *size_MatF];
   f0MPsource = &DF[dirTS  *size_MatF];
   f000source = &DF[dirZERO*size_MatF];
   fMMMsource = &DF[dirBSW *size_MatF];
   fMMPsource = &DF[dirTSW *size_MatF];
   fMPPsource = &DF[dirTNW *size_MatF];
   fMPMsource = &DF[dirBNW *size_MatF];
   fPPMsource = &DF[dirBNE *size_MatF];
   fPPPsource = &DF[dirTNE *size_MatF];
   fPMPsource = &DF[dirTSE *size_MatF];
   fPMMsource = &DF[dirBSE *size_MatF];

   real
	   *fP00dest, *fM00dest, *f0P0dest, *f0M0dest, *f00Pdest, *f00Mdest, *fPP0dest, *fMM0dest, *fPM0dest,
	   *fMP0dest, *fP0Pdest, *fM0Mdest, *fP0Mdest, *fM0Pdest, *f0PPdest, *f0MMdest, *f0PMdest, *f0MPdest,
	   *f000dest, *fMMMdest, *fMMPdest, *fMPPdest, *fMPMdest, *fPPMdest, *fPPPdest, *fPMPdest, *fPMMdest;

   if (evenOrOdd==true)
   {
	   fP00dest = &DC[dirE   *size_MatC];
	   fM00dest = &DC[dirW   *size_MatC];
	   f0P0dest = &DC[dirN   *size_MatC];
	   f0M0dest = &DC[dirS   *size_MatC];
	   f00Pdest = &DC[dirT   *size_MatC];
	   f00Mdest = &DC[dirB   *size_MatC];
	   fPP0dest = &DC[dirNE  *size_MatC];
	   fMM0dest = &DC[dirSW  *size_MatC];
	   fPM0dest = &DC[dirSE  *size_MatC];
	   fMP0dest = &DC[dirNW  *size_MatC];
	   fP0Pdest = &DC[dirTE  *size_MatC];
	   fM0Mdest = &DC[dirBW  *size_MatC];
	   fP0Mdest = &DC[dirBE  *size_MatC];
	   fM0Pdest = &DC[dirTW  *size_MatC];
	   f0PPdest = &DC[dirTN  *size_MatC];
	   f0MMdest = &DC[dirBS  *size_MatC];
	   f0PMdest = &DC[dirBN  *size_MatC];
	   f0MPdest = &DC[dirTS  *size_MatC];
	   f000dest = &DC[dirZERO*size_MatC];
	   fMMMdest = &DC[dirBSW *size_MatC];
	   fMMPdest = &DC[dirTSW *size_MatC];
	   fMPPdest = &DC[dirTNW *size_MatC];
	   fMPMdest = &DC[dirBNW *size_MatC];
	   fPPMdest = &DC[dirBNE *size_MatC];
	   fPPPdest = &DC[dirTNE *size_MatC];
	   fPMPdest = &DC[dirTSE *size_MatC];
	   fPMMdest = &DC[dirBSE *size_MatC];
   } 
   else
   {
	   fP00dest = &DC[dirW   *size_MatC];
	   fM00dest = &DC[dirE   *size_MatC];
	   f0P0dest = &DC[dirS   *size_MatC];
	   f0M0dest = &DC[dirN   *size_MatC];
	   f00Pdest = &DC[dirB   *size_MatC];
	   f00Mdest = &DC[dirT   *size_MatC];
	   fPP0dest = &DC[dirSW  *size_MatC];
	   fMM0dest = &DC[dirNE  *size_MatC];
	   fPM0dest = &DC[dirNW  *size_MatC];
	   fMP0dest = &DC[dirSE  *size_MatC];
	   fP0Pdest = &DC[dirBW  *size_MatC];
	   fM0Mdest = &DC[dirTE  *size_MatC];
	   fP0Mdest = &DC[dirTW  *size_MatC];
	   fM0Pdest = &DC[dirBE  *size_MatC];
	   f0PPdest = &DC[dirBS  *size_MatC];
	   f0MMdest = &DC[dirTN  *size_MatC];
	   f0PMdest = &DC[dirTS  *size_MatC];
	   f0MPdest = &DC[dirBN  *size_MatC];
	   f000dest = &DC[dirZERO*size_MatC];
	   fMMMdest = &DC[dirTNE *size_MatC];
	   fMMPdest = &DC[dirBNE *size_MatC];
	   fMPPdest = &DC[dirBSE *size_MatC];
	   fMPMdest = &DC[dirTSE *size_MatC];
	   fPPMdest = &DC[dirTSW *size_MatC];
	   fPPPdest = &DC[dirBSW *size_MatC];
	   fPMPdest = &DC[dirBNW *size_MatC];
	   fPMMdest = &DC[dirTNW *size_MatC];
   }

   Distributions6 G;
   if (evenOrOdd == true)
   {
	   G.g[dirE] = &G6[dirE   *size_MatC];
	   G.g[dirW] = &G6[dirW   *size_MatC];
	   G.g[dirN] = &G6[dirN   *size_MatC];
	   G.g[dirS] = &G6[dirS   *size_MatC];
	   G.g[dirT] = &G6[dirT   *size_MatC];
	   G.g[dirB] = &G6[dirB   *size_MatC];
   }
   else
   {
	   G.g[dirW] = &G6[dirE   *size_MatC];
	   G.g[dirE] = &G6[dirW   *size_MatC];
	   G.g[dirS] = &G6[dirN   *size_MatC];
	   G.g[dirN] = &G6[dirS   *size_MatC];
	   G.g[dirB] = &G6[dirT   *size_MatC];
	   G.g[dirT] = &G6[dirB   *size_MatC];
   }

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  ix = threadIdx.x;  // Globaler x-Index 
   const unsigned  iy = blockIdx.x;   // Globaler y-Index 
   const unsigned  iz = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*iz + iy) + ix;
   //////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   real eps_new = two;
   real omegaS = omFine;//-omFine;
   real o  = omCoarse;//-omCoarse;
   //real op = one;
   //real cu_sq;

   real xoff,    yoff,    zoff;
   real xoff_sq, yoff_sq, zoff_sq;

   real        vvx, vvy, vvz, vx2, vy2, vz2, drho;
   real        press;//,drho,vx1,vx2,vx3;
   real        /*pressMMM,*/drhoMMM,vx1MMM,vx2MMM,vx3MMM;
   real        /*pressMMP,*/drhoMMP,vx1MMP,vx2MMP,vx3MMP;
   real        /*pressMPP,*/drhoMPP,vx1MPP,vx2MPP,vx3MPP;
   real        /*pressMPM,*/drhoMPM,vx1MPM,vx2MPM,vx3MPM;
   real        /*pressPPM,*/drhoPPM,vx1PPM,vx2PPM,vx3PPM;
   real        /*pressPPP,*/drhoPPP,vx1PPP,vx2PPP,vx3PPP;
   real        /*pressPMP,*/drhoPMP,vx1PMP,vx2PMP,vx3PMP;
   real        /*pressPMM,*/drhoPMM,vx1PMM,vx2PMM,vx3PMM;
   real        fP00,fM00,f0P0,f0M0,f00P,f00M,fPP0,fMM0,fPM0,fMP0,fP0P,fM0M,fP0M,fM0P,f0PP,f0MM,f0PM,f0MP,f000,fPPP, fMMP, fPMP, fMPP, fPPM, fMMM, fPMM, fMPM;
   real        kxyFromfcNEQMMM, kyzFromfcNEQMMM, kxzFromfcNEQMMM, kxxMyyFromfcNEQMMM, kxxMzzFromfcNEQMMM, kyyMzzFromfcNEQMMM;
   real        kxyFromfcNEQMMP, kyzFromfcNEQMMP, kxzFromfcNEQMMP, kxxMyyFromfcNEQMMP, kxxMzzFromfcNEQMMP, kyyMzzFromfcNEQMMP;
   real        kxyFromfcNEQMPP, kyzFromfcNEQMPP, kxzFromfcNEQMPP, kxxMyyFromfcNEQMPP, kxxMzzFromfcNEQMPP, kyyMzzFromfcNEQMPP;
   real        kxyFromfcNEQMPM, kyzFromfcNEQMPM, kxzFromfcNEQMPM, kxxMyyFromfcNEQMPM, kxxMzzFromfcNEQMPM, kyyMzzFromfcNEQMPM;
   real        kxyFromfcNEQPPM, kyzFromfcNEQPPM, kxzFromfcNEQPPM, kxxMyyFromfcNEQPPM, kxxMzzFromfcNEQPPM, kyyMzzFromfcNEQPPM;
   real        kxyFromfcNEQPPP, kyzFromfcNEQPPP, kxzFromfcNEQPPP, kxxMyyFromfcNEQPPP, kxxMzzFromfcNEQPPP, kyyMzzFromfcNEQPPP;
   real        kxyFromfcNEQPMP, kyzFromfcNEQPMP, kxzFromfcNEQPMP, kxxMyyFromfcNEQPMP, kxxMzzFromfcNEQPMP, kyyMzzFromfcNEQPMP;
   real        kxyFromfcNEQPMM, kyzFromfcNEQPMM, kxzFromfcNEQPMM, kxxMyyFromfcNEQPMM, kxxMzzFromfcNEQPMM, kyyMzzFromfcNEQPMM;
   real        a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz, b0, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz, c0, cx, cy, cz, cxx, cyy, czz, cxy, cxz, cyz/*, axyz, bxyz, cxyz*/;
   real        d0, dx, dy, dz, dxy, dxz, dyz/*, dxyz*/;

   if(k<kFC)
   {
      //////////////////////////////////////////////////////////////////////////
      xoff = offFC.xOffFC[k];
      yoff = offFC.yOffFC[k];
      zoff = offFC.zOffFC[k];      
      xoff_sq = xoff * xoff;
      yoff_sq = yoff * yoff;
      zoff_sq = zoff * zoff;
      //////////////////////////////////////////////////////////////////////////
      //SWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      unsigned int k000base = posFSWB[k];
      unsigned int kM00base = neighborFX[k000base];
      unsigned int k0M0base = neighborFY[k000base];
      unsigned int k00Mbase = neighborFZ[k000base];
      unsigned int kMM0base = neighborFY[kM00base];
      unsigned int kM0Mbase = neighborFZ[kM00base];
      unsigned int k0MMbase = neighborFZ[k0M0base];
      unsigned int kMMMbase = neighborFZ[kMM0base];
      //////////////////////////////////////////////////////////////////////////
      //index 
      unsigned int k000 = k000base;
      unsigned int kM00 = kM00base;   
      unsigned int k0M0 = k0M0base;   
      unsigned int k00M = k00Mbase;   
      unsigned int kMM0 = kMM0base;  
      unsigned int kM0M = kM0Mbase;  
      unsigned int k0MM = k0MMbase;  
      unsigned int kMMM = kMMMbase; 
      ////////////////////////////////////////////////////////////////////////////////
      fP00 = fP00source[k000];
      fM00 = fM00source[kM00];
      f0P0 = f0P0source[k000];
      f0M0 = f0M0source[k0M0];
      f00P = f00Psource[k000];
      f00M = f00Msource[k00M];
      fPP0 = fPP0source[k000];
      fMM0 = fMM0source[kMM0];
      fPM0 = fPM0source[k0M0];
      fMP0 = fMP0source[kM00];
      fP0P = fP0Psource[k000];
      fM0M = fM0Msource[kM0M];
      fP0M = fP0Msource[k00M];
      fM0P = fM0Psource[kM00];
      f0PP = f0PPsource[k000];
      f0MM = f0MMsource[k0MM];
      f0PM = f0PMsource[k00M];
      f0MP = f0MPsource[k0M0];
      f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
      fPMP = fPMPsource[k0M0];
      fPMM = fPMMsource[k0MM];

      drhoMMM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MMM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMMM);
	  vx2MMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMMM);
	  vx3MMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMMM);

	  kxyFromfcNEQMMM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMMM) - ((vx1MMM*vx2MMM)));
	  kyzFromfcNEQMMM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMMM) - ((vx2MMM*vx3MMM)));
	  kxzFromfcNEQMMM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMMM) - ((vx1MMM*vx3MMM)));
	  kxxMyyFromfcNEQMMM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMMM) - ((vx1MMM*vx1MMM - vx2MMM*vx2MMM)));
	  kxxMzzFromfcNEQMMM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMMM) - ((vx1MMM*vx1MMM - vx3MMM*vx3MMM)));
	  kyyMzzFromfcNEQMMM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMMM) - ((vx2MMM*vx2MMM - vx3MMM*vx3MMM)));

      //////////////////////////////////////////////////////////////////////////
      //SWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k00M;
      kM00 = kM0M;   
      k0M0 = k0MM;   
      k00M = neighborFZ[k00M];   
      kMM0 = kMMM;  
      kM0M = neighborFZ[kM0M];  
      k0MM = neighborFZ[k0MM];  
      kMMM = neighborFZ[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMMP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MMP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMMP);
	  vx2MMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMMP);
	  vx3MMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMMP);

	  kxyFromfcNEQMMP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMMP) - ((vx1MMP*vx2MMP)));
	  kyzFromfcNEQMMP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMMP) - ((vx2MMP*vx3MMP)));
	  kxzFromfcNEQMMP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMMP) - ((vx1MMP*vx3MMP)));
	  kxxMyyFromfcNEQMMP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMMP) - ((vx1MMP*vx1MMP - vx2MMP*vx2MMP)));
	  kxxMzzFromfcNEQMMP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMMP) - ((vx1MMP*vx1MMP - vx3MMP*vx3MMP)));
	  kyyMzzFromfcNEQMMP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMMP) - ((vx2MMP*vx2MMP - vx3MMP*vx3MMP)));

      //////////////////////////////////////////////////////////////////////////
      //SET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = kM00;
      kM00 = neighborFX[kM00];   
      k0M0 = kMM0;   
      k00M = kM0M;   
      kMM0 = neighborFX[kMM0];  
      kM0M = neighborFX[kM0M];  
      k0MM = kMMM;  
      kMMM = neighborFX[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPMP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PMP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPMP);
	  vx2PMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPMP);
	  vx3PMP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPMP);

	  kxyFromfcNEQPMP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPMP) - ((vx1PMP*vx2PMP)));
	  kyzFromfcNEQPMP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPMP) - ((vx2PMP*vx3PMP)));
	  kxzFromfcNEQPMP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPMP) - ((vx1PMP*vx3PMP)));
	  kxxMyyFromfcNEQPMP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPMP) - ((vx1PMP*vx1PMP - vx2PMP*vx2PMP)));
	  kxxMzzFromfcNEQPMP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPMP) - ((vx1PMP*vx1PMP - vx3PMP*vx3PMP)));
	  kyyMzzFromfcNEQPMP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPMP) - ((vx2PMP*vx2PMP - vx3PMP*vx3PMP)));

      //////////////////////////////////////////////////////////////////////////
      //SEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k00M = k000;   
      kM0M = kM00;  
      k0MM = k0M0;  
      kMMM = kMM0; 
      k000 = kM00base;
      kM00 = neighborFX[kM00base];   
      k0M0 = kMM0base;   
      kMM0 = neighborFX[kMM0base];  
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPMM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PMM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPMM);
	  vx2PMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPMM);
	  vx3PMM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPMM);

	  kxyFromfcNEQPMM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPMM) - ((vx1PMM*vx2PMM)));
	  kyzFromfcNEQPMM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPMM) - ((vx2PMM*vx3PMM)));
	  kxzFromfcNEQPMM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPMM) - ((vx1PMM*vx3PMM)));
	  kxxMyyFromfcNEQPMM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPMM) - ((vx1PMM*vx1PMM - vx2PMM*vx2PMM)));
	  kxxMzzFromfcNEQPMM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPMM) - ((vx1PMM*vx1PMM - vx3PMM*vx3PMM)));
	  kyyMzzFromfcNEQPMM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPMM) - ((vx2PMM*vx2PMM - vx3PMM*vx3PMM)));

      //////////////////////////////////////////////////////////////////////////
      //NWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      k000base = k0M0base;
      kM00base = kMM0base;
      k0M0base = neighborFY[k0M0base];
      k00Mbase = k0MMbase;
      kMM0base = neighborFY[kMM0base];
      kM0Mbase = kMMMbase;
      k0MMbase = neighborFY[k0MMbase];
      kMMMbase = neighborFY[kMMMbase];
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k000base;
      kM00 = kM00base;   
      k0M0 = k0M0base;   
      k00M = k00Mbase;   
      kMM0 = kMM0base;  
      kM0M = kM0Mbase;  
      k0MM = k0MMbase;  
      kMMM = kMMMbase; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMPM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MPM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMPM);
	  vx2MPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMPM);
	  vx3MPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMPM);

	  kxyFromfcNEQMPM    = -three*omegaS*   ((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMPM) - ((vx1MPM*vx2MPM)));
	  kyzFromfcNEQMPM    = -three*omegaS*   ((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMPM) - ((vx2MPM*vx3MPM)));
	  kxzFromfcNEQMPM    = -three*omegaS*   ((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMPM) - ((vx1MPM*vx3MPM)));
	  kxxMyyFromfcNEQMPM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMPM) - ((vx1MPM*vx1MPM - vx2MPM*vx2MPM)));
	  kxxMzzFromfcNEQMPM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMPM) - ((vx1MPM*vx1MPM - vx3MPM*vx3MPM)));
	  kyyMzzFromfcNEQMPM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMPM) - ((vx2MPM*vx2MPM - vx3MPM*vx3MPM)));

	  //////////////////////////////////////////////////////////////////////////
      //NWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = k00M;
      kM00 = kM0M;   
      k0M0 = k0MM;   
      k00M = neighborFZ[k00M];   
      kMM0 = kMMM;  
      kM0M = neighborFZ[kM0M];  
      k0MM = neighborFZ[k0MM];  
      kMMM = neighborFZ[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoMPP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1MPP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoMPP);
	  vx2MPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoMPP);
	  vx3MPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoMPP);

	  kxyFromfcNEQMPP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoMPP) - ((vx1MPP*vx2MPP)));
	  kyzFromfcNEQMPP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoMPP) - ((vx2MPP*vx3MPP)));
	  kxzFromfcNEQMPP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoMPP) - ((vx1MPP*vx3MPP)));
	  kxxMyyFromfcNEQMPP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoMPP) - ((vx1MPP*vx1MPP - vx2MPP*vx2MPP)));
	  kxxMzzFromfcNEQMPP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoMPP) - ((vx1MPP*vx1MPP - vx3MPP*vx3MPP)));
	  kyyMzzFromfcNEQMPP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoMPP) - ((vx2MPP*vx2MPP - vx3MPP*vx3MPP)));

      //////////////////////////////////////////////////////////////////////////
      //NET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k000 = kM00;
      kM00 = neighborFX[kM00];   
      k0M0 = kMM0;   
      k00M = kM0M;   
      kMM0 = neighborFX[kMM0];  
      kM0M = neighborFX[kM0M];  
      k0MM = kMMM;  
      kMMM = neighborFX[kMMM]; 
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPPP = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PPP  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPPP);
	  vx2PPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPPP);
	  vx3PPP  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPPP);

	  kxyFromfcNEQPPP    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPPP) - ((vx1PPP*vx2PPP)));
	  kyzFromfcNEQPPP    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPPP) - ((vx2PPP*vx3PPP)));
	  kxzFromfcNEQPPP    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPPP) - ((vx1PPP*vx3PPP)));
	  kxxMyyFromfcNEQPPP = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPPP) - ((vx1PPP*vx1PPP - vx2PPP*vx2PPP)));
	  kxxMzzFromfcNEQPPP = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPPP) - ((vx1PPP*vx1PPP - vx3PPP*vx3PPP)));
	  kyyMzzFromfcNEQPPP = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPPP) - ((vx2PPP*vx2PPP - vx3PPP*vx3PPP)));

      //////////////////////////////////////////////////////////////////////////
      //NEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      k00M = k000;   
      kM0M = kM00;  
      k0MM = k0M0;  
      kMMM = kMM0; 
      k000 = kM00base;
      kM00 = neighborFX[kM00base];   
      k0M0 = kMM0base;   
      kMM0 = neighborFX[kMM0base];  
      ////////////////////////////////////////////////////////////////////////////////
	  fP00 = fP00source[k000];
	  fM00 = fM00source[kM00];
	  f0P0 = f0P0source[k000];
	  f0M0 = f0M0source[k0M0];
	  f00P = f00Psource[k000];
	  f00M = f00Msource[k00M];
	  fPP0 = fPP0source[k000];
	  fMM0 = fMM0source[kMM0];
	  fPM0 = fPM0source[k0M0];
	  fMP0 = fMP0source[kM00];
	  fP0P = fP0Psource[k000];
	  fM0M = fM0Msource[kM0M];
	  fP0M = fP0Msource[k00M];
	  fM0P = fM0Psource[kM00];
	  f0PP = f0PPsource[k000];
	  f0MM = f0MMsource[k0MM];
	  f0PM = f0PMsource[k00M];
	  f0MP = f0MPsource[k0M0];
	  f000 = f000source[k000];
	  fMMM = fMMMsource[kMMM];
	  fMMP = fMMPsource[kMM0];
	  fMPP = fMPPsource[kM00];
	  fMPM = fMPMsource[kM0M];
	  fPPM = fPPMsource[k00M];
	  fPPP = fPPPsource[k000];
	  fPMP = fPMPsource[k0M0];
	  fPMM = fPMMsource[k0MM];

      drhoPPM = fP00+fM00+f0P0+f0M0+f00P+f00M+fPP0+fMM0+fPM0+fMP0+fP0P+fM0M+fP0M+fM0P+f0PP+f0MM+f0PM+f0MP+f000+fPPP+fMMP+fPMP+fMPP+fPPM+fMMM+fPMM+fMPM;
      vx1PPM  = (((fPPP-fMMM)+(fPMP-fMPM)+(fPPM-fMMP)+(fPMM-fMPP)) + (((fPP0-fMM0)+(fP0P-fM0M))+((fPM0-fMP0)+(fP0M-fM0P))) + (fP00-fM00))/(one + drhoPPM);
	  vx2PPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPPM-fMMP)+(fMPM-fPMP)) + (((fPP0-fMM0)+(f0PP-f0MM))+((f0PM-f0MP)+(fMP0-fPM0))) + (f0P0-f0M0))/(one + drhoPPM);
	  vx3PPM  = (((fPPP-fMMM)+(fMPP-fPMM)+(fPMP-fMPM)+(fMMP-fPPM)) + (((fP0P-fM0M)+(f0PP-f0MM))+((fM0P-fP0M)+(f0MP-f0PM))) + (f00P-f00M))/(one + drhoPPM);

	  kxyFromfcNEQPPM    = -three*omegaS*((((fMM0 - fPM0) + (fPP0 - fMP0)) + (((fMMM - fPMM) + (fPPM - fMPM)) + ((fMMP - fPMP) + (fPPP - fMPP)))) / (one + drhoPPM) - ((vx1PPM*vx2PPM)));
	  kyzFromfcNEQPPM    = -three*omegaS*((((f0MM - f0MP) + (f0PP - f0PM)) + (((fPMM - fPMP) + (fMMM - fMPM)) + ((fPPP - fPPM) + (fMPP - fMMP)))) / (one + drhoPPM) - ((vx2PPM*vx3PPM)));
	  kxzFromfcNEQPPM    = -three*omegaS*((((fM0M - fP0M) + (fP0P - fM0P)) + (((fMMM - fPMM) + (fMPM - fPPM)) + ((fPMP - fMMP) + (fPPP - fMPP)))) / (one + drhoPPM) - ((vx1PPM*vx3PPM)));
	  kxxMyyFromfcNEQPPM = -c3o2*omegaS *(((((fM0M - f0MM) + (fM0P - f0MP)) + ((fP0M - f0PM) + (fP0P - f0PP))) + ((fM00 - f0M0) + (fP00 - f0P0))) / (one + drhoPPM) - ((vx1PPM*vx1PPM - vx2PPM*vx2PPM)));
	  kxxMzzFromfcNEQPPM = -c3o2*omegaS *(((((fMM0 - f0MM) + (fMP0 - f0PM)) + ((fPM0 - f0MP) + (fPP0 - f0PP))) + ((fM00 - f00M) + (fP00 - f00P))) / (one + drhoPPM) - ((vx1PPM*vx1PPM - vx3PPM*vx3PPM)));
	  kyyMzzFromfcNEQPPM = -c3o2*omegaS *(((((fPM0 - fP0M) + (fMM0 - fM0M)) + ((fPP0 - fP0P) + (fMP0 - fM0P))) + ((f0M0 - f00M) + (f0P0 - f00P))) / (one + drhoPPM) - ((vx2PPM*vx2PPM - vx3PPM*vx3PPM)));

      //////////////////////////////////////////////////////////////////////////
      //3
      //////////////////////////////////////////////////////////////////////////
	  a0  = c1o8*(((vx1PPM + vx1MMP) + (vx1MPM + vx1PMP)) + ((vx1PMM + vx1MPP) + (vx1MMM + vx1PPP)));
	  ax  = c1o4*(((vx1PPM - vx1MMP) + (vx1PPP - vx1MMM)) + ((vx1PMP - vx1MPM) + (vx1PMM - vx1MPP)));
	  ay  = c1o4*(((vx1PPM - vx1MMP) + (vx1PPP - vx1MMM)) + ((vx1MPM - vx1PMP) + (vx1MPP - vx1PMM)));
	  az  = c1o4*(((vx1MMP - vx1PPM) + (vx1PPP - vx1MMM)) + ((vx1PMP - vx1MPM) + (vx1MPP - vx1PMM)));
	  axy = c1o2*(((vx1PPM - vx1PMP) + (vx1MMM - vx1MPP)) + ((vx1MMP - vx1MPM) + (vx1PPP - vx1PMM)));
	  axz = c1o2*(((vx1PMP - vx1PPM) + (vx1MMM - vx1MPP)) + ((vx1MPM - vx1MMP) + (vx1PPP - vx1PMM)));
	  ayz = c1o2*(((vx1PPP - vx1MPM) + (vx1PMM - vx1MMP)) + ((vx1MPP - vx1PPM) + (vx1MMM - vx1PMP)));
	  //axyz=		  ((vx1MMP - vx1PPM) + (vx1PPP - vx1MMM)) + ((vx1MPM - vx1PMP) + (vx1PMM - vx1MPP)) ;
	  b0  = c1o8*(((vx2PPM + vx2MMP) + (vx2MPM + vx2PMP)) + ((vx2PMM + vx2MPP) + (vx2MMM + vx2PPP)));
	  bx  = c1o4*(((vx2PPM - vx2MMP) + (vx2PPP - vx2MMM)) + ((vx2PMP - vx2MPM) + (vx2PMM - vx2MPP)));
	  by  = c1o4*(((vx2PPM - vx2MMP) + (vx2PPP - vx2MMM)) + ((vx2MPM - vx2PMP) + (vx2MPP - vx2PMM)));
	  bz  = c1o4*(((vx2MMP - vx2PPM) + (vx2PPP - vx2MMM)) + ((vx2PMP - vx2MPM) + (vx2MPP - vx2PMM)));
	  bxy = c1o2*(((vx2PPM - vx2MPP) + (vx2MMM - vx2PMP)) + ((vx2MMP - vx2PMM) + (vx2PPP - vx2MPM)));
	  bxz = c1o2*(((vx2MMM - vx2PPM) + (vx2PMP - vx2MPP)) + ((vx2MPM - vx2PMM) + (vx2PPP - vx2MMP)));
	  byz = c1o2*(((vx2MPP - vx2PPM) + (vx2MMM - vx2PMP)) + ((vx2PMM - vx2MMP) + (vx2PPP - vx2MPM)));
	  //bxyz=		  ((vx2MMP - vx2PPM) + (vx2PPP - vx2MMM)) + ((vx2MPM - vx2PMP) + (vx2PMM - vx2MPP)) ;
	  c0  = c1o8*(((vx3PPM + vx3MMP) + (vx3MPM + vx3PMP)) + ((vx3PMM + vx3MPP) + (vx3MMM + vx3PPP)));
	  cx  = c1o4*(((vx3PPM - vx3MMP) + (vx3PPP - vx3MMM)) + ((vx3PMP - vx3MPM) + (vx3PMM - vx3MPP)));
	  cy  = c1o4*(((vx3PPM - vx3MMP) + (vx3PPP - vx3MMM)) + ((vx3MPM - vx3PMP) + (vx3MPP - vx3PMM)));
	  cz  = c1o4*(((vx3MMP - vx3PPM) + (vx3PPP - vx3MMM)) + ((vx3PMP - vx3MPM) + (vx3MPP - vx3PMM)));
	  cxy = c1o2*(((vx3PPM - vx3PMP) + (vx3MMM - vx3MPP)) + ((vx3MMP - vx3MPM) + (vx3PPP - vx3PMM)));
	  cxz = c1o2*(((vx3MMM - vx3PPM) + (vx3PMP - vx3MPP)) + ((vx3MPM - vx3PMM) + (vx3PPP - vx3MMP)));
	  cyz = c1o2*(((vx3MMM - vx3PPM) + (vx3MPP - vx3PMP)) + ((vx3PMM - vx3MPM) + (vx3PPP - vx3MMP)));
	  //cxyz=		  ((vx3MMP - vx3PPM) + (vx3PPP - vx3MMM)) + ((vx3MPM - vx3PMP) + (vx3PMM - vx3MPP)) ;

	  //bxx = zero;
	  //cxx = zero;
	  //ayy = zero;
	  //cyy = zero;
	  //azz = zero;
	  //bzz = zero;
	  //axx = zero;
	  //byy = zero;
	  //czz = zero;

	  bxx = c1o8*(((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPMM - kxyFromfcNEQMPP)) + ((kxyFromfcNEQPMP - kxyFromfcNEQMPM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP))) - c1o2*axy;
	  cxx = c1o8*(((kxzFromfcNEQPPP - kxzFromfcNEQMMM) + (kxzFromfcNEQPMM - kxzFromfcNEQMPP)) + ((kxzFromfcNEQPMP - kxzFromfcNEQMPM) + (kxzFromfcNEQPPM - kxzFromfcNEQMMP))) - c1o2*axz;

	  ayy = c1o8*(((kxyFromfcNEQPPP - kxyFromfcNEQMMM) + (kxyFromfcNEQPPM - kxyFromfcNEQMMP)) - ((kxyFromfcNEQPMM - kxyFromfcNEQMPP) + (kxyFromfcNEQPMP - kxyFromfcNEQMPM))) - c1o2*bxy;
	  cyy = c1o8*(((kyzFromfcNEQPPP - kyzFromfcNEQMMM) + (kyzFromfcNEQPPM - kyzFromfcNEQMMP)) - ((kyzFromfcNEQPMM - kyzFromfcNEQMPP) + (kyzFromfcNEQPMP - kyzFromfcNEQMPM))) - c1o2*byz;

	  azz = c1o8*(((kxzFromfcNEQPPP - kxzFromfcNEQMMM) - (kxzFromfcNEQPMM - kxzFromfcNEQMPP)) + ((kxzFromfcNEQPMP - kxzFromfcNEQMPM) - (kxzFromfcNEQPPM - kxzFromfcNEQMMP))) - c1o2*cxz;
	  bzz = c1o8*(((kyzFromfcNEQPPP - kyzFromfcNEQMMM) - (kyzFromfcNEQPMM - kyzFromfcNEQMPP)) + ((kyzFromfcNEQPMP - kyzFromfcNEQMPM) - (kyzFromfcNEQPPM - kyzFromfcNEQMMP))) - c1o2*cyz;

	  axx = ( c1o16*(((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) + (kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP)) + ((kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM) + (kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP))) - c1o4*bxy)
		  + ( c1o16*(((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) + (kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP)) + ((kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM) + (kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP))) - c1o4*cxz);

	  byy = (-c1o16*(((kxxMyyFromfcNEQPPP - kxxMyyFromfcNEQMMM) - (kxxMyyFromfcNEQPMM - kxxMyyFromfcNEQMPP)) + ((kxxMyyFromfcNEQPPM - kxxMyyFromfcNEQMMP) - (kxxMyyFromfcNEQPMP - kxxMyyFromfcNEQMPM))) - c1o4*axy)
		  + ( c1o16*(((kyyMzzFromfcNEQPPP - kyyMzzFromfcNEQMMM) - (kyyMzzFromfcNEQPMM - kyyMzzFromfcNEQMPP)) + ((kyyMzzFromfcNEQPPM - kyyMzzFromfcNEQMMP) - (kyyMzzFromfcNEQPMP - kyyMzzFromfcNEQMPM))) - c1o4*cyz);

	  czz = (-c1o16*(((kxxMzzFromfcNEQPPP - kxxMzzFromfcNEQMMM) - (kxxMzzFromfcNEQPMM - kxxMzzFromfcNEQMPP)) - ((kxxMzzFromfcNEQPPM - kxxMzzFromfcNEQMMP) - (kxxMzzFromfcNEQPMP - kxxMzzFromfcNEQMPM))) - c1o4*axz)
		  + ( c1o16*(((kyyMzzFromfcNEQPPP - kyyMzzFromfcNEQMMM) - (kyyMzzFromfcNEQPMM - kyyMzzFromfcNEQMPP)) - ((kyyMzzFromfcNEQPPM - kyyMzzFromfcNEQMMP) - (kyyMzzFromfcNEQPMP - kyyMzzFromfcNEQMPM))) - c1o4*byz);

	  a0 -= c1o4*(axx + ayy + azz);
	  b0 -= c1o4*(bxx + byy + bzz);
	  c0 -= c1o4*(cxx + cyy + czz);

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real kxyAverage = zero;
	  real kyzAverage = zero;
	  real kxzAverage = zero;
	  real kxxMyyAverage = zero;
	  real kxxMzzAverage = zero;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////Press
	  //d0   = ( pressPPM + pressPPP + pressMPM + pressMPP + pressPMM + pressPMP + pressMMM + pressMMP) * c1o8;
	  //dx   = ( pressPPM + pressPPP - pressMPM - pressMPP + pressPMM + pressPMP - pressMMM - pressMMP) * c1o4;
	  //dy   = ( pressPPM + pressPPP + pressMPM + pressMPP - pressPMM - pressPMP - pressMMM - pressMMP) * c1o4;
	  //dz   = (-pressPPM + pressPPP - pressMPM + pressMPP - pressPMM + pressPMP - pressMMM + pressMMP) * c1o4;
	  //dxy  = ( pressPPM + pressPPP - pressMPM - pressMPP - pressPMM - pressPMP + pressMMM + pressMMP) * c1o2;
	  //dxz  = (-pressPPM + pressPPP + pressMPM - pressMPP - pressPMM + pressPMP + pressMMM - pressMMP) * c1o2;
	  //dyz  = (-pressPPM + pressPPP - pressMPM + pressMPP + pressPMM - pressPMP + pressMMM - pressMMP) * c1o2;
	  //dxyz =  -pressPPM + pressPPP + pressMPM - pressMPP + pressPMM - pressPMP - pressMMM + pressMMP;
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //drho
	  d0   = ( ((drhoPPM + drhoMMP) + (drhoMPM + drhoPMP)) + ((drhoPMM + drhoMPP) + (drhoMMM + drhoPPP))) * c1o8;
	  dx   = ( ((drhoPPM - drhoMMP) + (drhoPMM - drhoMPP)) + ((drhoPMP - drhoMPM) + (drhoPPP - drhoMMM))) * c1o4;
	  dy   = ( ((drhoPPM - drhoMMP) + (drhoMPP - drhoPMM)) + ((drhoMPM - drhoPMP) + (drhoPPP - drhoMMM))) * c1o4;
	  dz   = ( ((drhoMMP - drhoPPM) + (drhoPPP - drhoMMM)) + ((drhoPMP - drhoMPM) + (drhoMPP - drhoPMM))) * c1o4;
	  dxy  = ( ((drhoPPM - drhoPMP) + (drhoPPP - drhoPMM)) + ((drhoMMP - drhoMPM) + (drhoMMM - drhoMPP))) * c1o2;
	  dxz  = ( ((drhoMMM - drhoPPM) + (drhoPPP - drhoMMP)) + ((drhoMPM - drhoPMM) + (drhoPMP - drhoMPP))) * c1o2;
	  dyz  = ( ((drhoMPP - drhoPPM) + (drhoPPP - drhoMPM)) + ((drhoPMM - drhoMMP) + (drhoMMM - drhoPMP))) * c1o2;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // Bernd das Brot 
	  //
      //
	  // x------x
	  // |      |
	  // |	 ---+--->X
	  // |		|  \
	  // x------x   \
	  //			off-vector
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      a0 = a0 + xoff * ax + yoff * ay + zoff * az + xoff_sq * axx + yoff_sq * ayy + zoff_sq * azz + xoff*yoff*axy + xoff*zoff*axz + yoff*zoff*ayz;
      ax = ax + two * xoff * axx + yoff * axy + zoff * axz;
      ay = ay + two * yoff * ayy + xoff * axy + zoff * ayz;
      az = az + two * zoff * azz + xoff * axz + yoff * ayz;
      b0 = b0 + xoff * bx + yoff * by + zoff * bz + xoff_sq * bxx + yoff_sq * byy + zoff_sq * bzz + xoff*yoff*bxy + xoff*zoff*bxz + yoff*zoff*byz;
      bx = bx + two * xoff * bxx + yoff * bxy + zoff * bxz;
      by = by + two * yoff * byy + xoff * bxy + zoff * byz;
      bz = bz + two * zoff * bzz + xoff * bxz + yoff * byz;
      c0 = c0 + xoff * cx + yoff * cy + zoff * cz + xoff_sq * cxx + yoff_sq * cyy + zoff_sq * czz + xoff*yoff*cxy + xoff*zoff*cxz + yoff*zoff*cyz;
      cx = cx + two * xoff * cxx + yoff * cxy + zoff * cxz;
      cy = cy + two * yoff * cyy + xoff * cxy + zoff * cyz;
      cz = cz + two * zoff * czz + xoff * cxz + yoff * cyz;
	  d0 = d0 + xoff * dx + yoff * dy + zoff * dz + xoff*yoff*dxy + xoff*zoff*dxz + yoff*zoff*dyz;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real mfcbb = zero;
	  real mfabb = zero;
	  real mfbcb = zero;
	  real mfbab = zero;
	  real mfbbc = zero;
	  real mfbba = zero;
	  real mfccb = zero;
	  real mfaab = zero;
	  real mfcab = zero;
	  real mfacb = zero;
	  real mfcbc = zero;
	  real mfaba = zero;
	  real mfcba = zero;
	  real mfabc = zero;
	  real mfbcc = zero;
	  real mfbaa = zero;
	  real mfbca = zero;
	  real mfbac = zero;
	  real mfbbb = zero;
	  real mfccc = zero;
	  real mfaac = zero;
	  real mfcac = zero;
	  real mfacc = zero;
	  real mfcca = zero;
	  real mfaaa = zero;
	  real mfcaa = zero;
	  real mfaca = zero;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real mgcbb = zero;
	  real mgabb = zero;
	  real mgbcb = zero;
	  real mgbab = zero;
	  real mgbbc = zero;
	  real mgbba = zero;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  real m0, m1, m2, oMdrho;
	  real mxxPyyPzz, mxxMyy, mxxMzz, mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
	  real qudricLimit = c1o100;//ganz schlechte Idee -> muss global sein
	  real O3 = two - o;
	  real residu, residutmp;
	  residutmp = zero;///*-*/ c2o9 * (1./o - c1o2) * eps_new * eps_new;
	  real NeqOn = one;//zero;//one;   //.... one = on ..... zero = off 
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //
	  //Position C 0., 0., 0.
	  //
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //x = 0.;
	  //y = 0.;
	  //z = 0.;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //real mxoff = -xoff;
	  //real myoff = -yoff;
	  //real mzoff = -zoff;
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //press = pressPPP * (c1o8 - c1o4 * mxoff - c1o4 * myoff - c1o4 * mzoff) + 
			//  pressMPP * (c1o8 + c1o4 * mxoff - c1o4 * myoff - c1o4 * mzoff) + 
			//  pressPMP * (c1o8 - c1o4 * mxoff + c1o4 * myoff - c1o4 * mzoff) + 
			//  pressMMP * (c1o8 + c1o4 * mxoff + c1o4 * myoff - c1o4 * mzoff) + 
			//  pressPPM * (c1o8 - c1o4 * mxoff - c1o4 * myoff + c1o4 * mzoff) + 
			//  pressMPM * (c1o8 + c1o4 * mxoff - c1o4 * myoff + c1o4 * mzoff) + 
			//  pressPMM * (c1o8 - c1o4 * mxoff + c1o4 * myoff + c1o4 * mzoff) + 
			//  pressMMM * (c1o8 + c1o4 * mxoff + c1o4 * myoff + c1o4 * mzoff);
	  //drho  = drhoPPP * (c1o8 - c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) + 
			//  drhoMPP * (c1o8 + c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) + 
			//  drhoPMP * (c1o8 - c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) + 
			//  drhoMMP * (c1o8 + c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) + 
			//  drhoPPM * (c1o8 - c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) + 
			//  drhoMPM * (c1o8 + c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) + 
			//  drhoPMM * (c1o8 - c1o4 * xoff + c1o4 * yoff + c1o4 * zoff) + 
			//  drhoMMM * (c1o8 + c1o4 * xoff + c1o4 * yoff + c1o4 * zoff);
	  press = d0;
	  vvx   = a0;
	  vvy   = b0;
	  vvz   = c0;

	  //mfaaa = drho;
	  //mfaaa = press + (ax+by+cz)/three;  //  1/3 = 2/3*(1/op-1/2)
	  mfaaa = press; // if drho is interpolated directly

	  vx2 = vvx*vvx;
	  vy2 = vvy*vvy;
	  vz2 = vvz*vvz;
	  oMdrho = one;
	  //oMdrho = one - mfaaa;

	  //two
	  // linear combinations
	  real oP = o;// :(
	  mxxPyyPzz = mfaaa    -c2o3*(ax+by+cz)*eps_new/oP*(one+press); 
	  //mxxMyy    = -c2o3*(ax - by)*eps_new/o;
	  //mxxMzz    = -c2o3*(ax - cz)*eps_new/o;

	  //mfabb     = -c1o3 * (bz + cy)*eps_new/o;
	  //mfbab     = -c1o3 * (az + cx)*eps_new/o;
	  //mfbba     = -c1o3 * (ay + bx)*eps_new/o;
	  mxxMyy    = -c2o3*((ax - by)+kxxMyyAverage)*eps_new/o * (one + press);
	  mxxMzz    = -c2o3*((ax - cz)+kxxMzzAverage)*eps_new/o * (one + press);

	  mfabb     = -c1o3 * ((bz + cy)+kyzAverage)*eps_new/o * (one + press);
	  mfbab     = -c1o3 * ((az + cx)+kxzAverage)*eps_new/o * (one + press);
	  mfbba     = -c1o3 * ((ay + bx)+kxyAverage)*eps_new/o * (one + press);

	  
	  // linear combinations back
	  mfcaa = c1o3 * (       mxxMyy +       mxxMzz + mxxPyyPzz) * NeqOn;
	  mfaca = c1o3 * (-two * mxxMyy +       mxxMzz + mxxPyyPzz) * NeqOn;
	  mfaac = c1o3 * (       mxxMyy - two * mxxMzz + mxxPyyPzz) * NeqOn;

	  //3.
	  // linear combinations
	  //residu = residutmp * (ayz + bxz + cxy );
	  //mfbbb = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mfbbb = zero;

	  //residu = residutmp * (axy + two*bxx + two*bzz + cyz );
	  //residu = -(c1o9*(axy - 2*bxx - 2*bzz + cyz ));
	  //mxxyPyzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxyPyzz = zero;

	  //residu = residutmp * (axy + two*bxx - two*bzz - cyz );
	  //residu = c1o9*(axy - 2*bxx + 2*bzz - cyz );
	  //mxxyMyzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxyMyzz = zero;

	  //residu = residutmp * (axz + byz + two*cxx + two*cyy );
	  //residu = -(c1o9*(axz + byz - 2*cxx - 2*cyy ));
	  //mxxzPyyz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxzPyyz = zero;

	  //residu = residutmp * (axz - byz + two*cxx - two*cyy );
	  //residu = c1o9*(axz - byz - 2*cxx + 2*cyy );
	  //mxxzMyyz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxxzMyyz = zero;

	  //residu = residutmp * (two*ayy + two*azz + bxy + cxz );
	  //residu = c1o9*(2*ayy + 2*azz - bxy - cxz );
	  //mxyyPxzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxyyPxzz = zero;

	  //residu = residutmp * (two*ayy - two*azz + bxy - cxz );
	  //residu = c1o9*(-2*ayy + 2*azz + bxy - cxz );
	  //mxyyMxzz = (abs(residu)+qudricLimit) * residu / (qudricLimit * O3 + abs(residu));
	  mxyyMxzz = zero;

	  ////////////////////////////////////////////////////////////////////////////////////
	  // D3Q27F 
	  mgcbb = (ax - four * axx) * eps_new;
	  mgabb = (ax + four * axx) * eps_new;
	  mgbcb = (by - four * byy) * eps_new;
	  mgbab = (by + four * byy) * eps_new;
	  mgbbc = (cz - four * czz) * eps_new;
	  mgbba = (cz + four * czz) * eps_new;
	  ////////////////////////////////////////////////////////////////////////////////////

	  // linear combinations back
	  mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
	  mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
	  mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
	  mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
	  mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
	  mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

	  //4.
	  mfacc = mfaaa*c1o9; 
	  mfcac = mfacc; 
	  mfcca = mfacc; 
	  //5.

	  //6.
	  mfccc = mfaaa*c1o27;
	  ////////////////////////////////////////////////////////////////////////////////////
	  //back
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // Z - Dir
	  m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfaac        - two * mfaab *  vvz         +  mfaaa                * (one - vz2)              - one * oMdrho * vz2; 
	  m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfaaa = m0;
	  mfaab = m1;
	  mfaac = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfabc        - two * mfabb *  vvz         + mfaba * (one - vz2); 
	  m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
	  mfaba = m0;
	  mfabb = m1;
	  mfabc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfacc        - two * mfacb *  vvz         +  mfaca                  * (one - vz2)              - c1o3 * oMdrho * vz2; 
	  m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfaca = m0;
	  mfacb = m1;
	  mfacc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbac        - two * mfbab *  vvz         + mfbaa * (one - vz2); 
	  m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
	  mfbaa = m0;
	  mfbab = m1;
	  mfbac = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbbc        - two * mfbbb *  vvz         + mfbba * (one - vz2); 
	  m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
	  mfbba = m0;
	  mfbbb = m1;
	  mfbbc = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
	  m1 = -mfbcc        - two * mfbcb *  vvz         + mfbca * (one - vz2); 
	  m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
	  mfbca = m0;
	  mfbcb = m1;
	  mfbcc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfcac        - two * mfcab *  vvz         +  mfcaa                  * (one - vz2)              - c1o3 * oMdrho * vz2; 
	  m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfcaa = m0;
	  mfcab = m1;
	  mfcac = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
	  m1 = -mfcbc        - two * mfcbb *  vvz         + mfcba * (one - vz2); 
	  m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
	  mfcba = m0;
	  mfcbb = m1;
	  mfcbc = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
	  m1 = -mfccc        - two * mfccb *  vvz         +  mfcca                  * (one - vz2)              - c1o9 * oMdrho * vz2; 
	  m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
	  mfcca = m0;
	  mfccb = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // Y - Dir
	  m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfaca        - two * mfaba *  vvy         +  mfaaa                  * (one - vy2)              - c1o6 * oMdrho * vy2; 
	  m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaaa = m0;
	  mfaba = m1;
	  mfaca = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfacb        - two * mfabb *  vvy         +  mfaab                  * (one - vy2)              - c2o3 * oMdrho * vy2; 
	  m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaab = m0;
	  mfabb = m1;
	  mfacb = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfacc        - two * mfabc *  vvy         +  mfaac                  * (one - vy2)              - c1o6 * oMdrho * vy2; 
	  m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfaac = m0;
	  mfabc = m1;
	  mfacc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbca        - two * mfbba *  vvy         + mfbaa * (one - vy2); 
	  m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
	  mfbaa = m0;
	  mfbba = m1;
	  mfbca = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbcb        - two * mfbbb *  vvy         + mfbab * (one - vy2); 
	  m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
	  mfbab = m0;
	  mfbbb = m1;
	  mfbcb = m2;
	  /////////b//////////////////////////////////////////////////////////////////////////
	  m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
	  m1 = -mfbcc        - two * mfbbc *  vvy         + mfbac * (one - vy2); 
	  m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
	  mfbac = m0;
	  mfbbc = m1;
	  mfbcc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfcca        - two * mfcba *  vvy         +  mfcaa                   * (one - vy2)              - c1o18 * oMdrho * vy2; 
	  m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcaa = m0;
	  mfcba = m1;
	  mfcca = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfccb        - two * mfcbb *  vvy         +  mfcab                  * (one - vy2)              - c2o9 * oMdrho * vy2; 
	  m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcab = m0;
	  mfcbb = m1;
	  mfccb = m2;
	  /////////c//////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
	  m1 = -mfccc        - two * mfcbc *  vvy         +  mfcac                   * (one - vy2)              - c1o18 * oMdrho * vy2; 
	  m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
	  mfcac = m0;
	  mfcbc = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
	  ////////////////////////////////////////////////////////////////////////////////////
	  // X - Dir
	  m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcaa        - two * mfbaa *  vvx         +  mfaaa                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaaa = m0;
	  mfbaa = m1;
	  mfcaa = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcba        - two * mfbba *  vvx         +  mfaba                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaba = m0;
	  mfbba = m1;
	  mfcba = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcca        - two * mfbca *  vvx         +  mfaca                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaca = m0;
	  mfbca = m1;
	  mfcca = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcab        - two * mfbab *  vvx         +  mfaab                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaab = m0;
	  mfbab = m1;
	  mfcab = m2;
	  ///////////b////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcbb        - two * mfbbb *  vvx         +  mfabb                  * (one - vx2)              - c4o9 * oMdrho * vx2; 
	  m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfabb = m0;
	  mfbbb = m1;
	  mfcbb = m2;
	  ///////////b////////////////////////////////////////////////////////////////////////
	  m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfccb        - two * mfbcb *  vvx         +  mfacb                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfacb = m0;
	  mfbcb = m1;
	  mfccb = m2;
	  ////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	  m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcac        - two * mfbac *  vvx         +  mfaac                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfaac = m0;
	  mfbac = m1;
	  mfcac = m2;
	  ///////////c////////////////////////////////////////////////////////////////////////
	  m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfcbc        - two * mfbbc *  vvx         +  mfabc                  * (one - vx2)              - c1o9 * oMdrho * vx2; 
	  m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfabc = m0;
	  mfbbc = m1;
	  mfcbc = m2;
	  ///////////c////////////////////////////////////////////////////////////////////////
	  m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
	  m1 = -mfccc        - two * mfbcc *  vvx         +  mfacc                   * (one - vx2)              - c1o36 * oMdrho * vx2; 
	  m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
	  mfacc = m0;
	  mfbcc = m1;
	  mfccc = m2;
	  ////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////
	  //index 0
	  k000 = posC[k];
	  kM00 = neighborCX[k000];
	  k0M0 = neighborCY[k000];
	  k00M = neighborCZ[k000];
	  kMM0 = neighborCY[kM00];
	  kM0M = neighborCZ[kM00];
	  k0MM = neighborCZ[k0M0];
	  kMMM = neighborCZ[kMM0];
	  ////////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////////
	  (G.g[dirE])[k000] = mgcbb;
	  (G.g[dirW])[kM00] = mgabb;
	  (G.g[dirN])[k000] = mgbcb;
	  (G.g[dirS])[k0M0] = mgbab;
	  (G.g[dirT])[k000] = mgbbc;
	  (G.g[dirB])[k00M] = mgbba;
	  ////////////////////////////////////////////////////////////////////////////////////
	  fP00dest[k000] = mfcbb;                                                                 
	  fM00dest[kM00] = mfabb;                                                               
	  f0P0dest[k000] = mfbcb;
	  f0M0dest[k0M0] = mfbab;
	  f00Pdest[k000] = mfbbc;
	  f00Mdest[k00M] = mfbba;
	  fPP0dest[k000] = mfccb;
	  fMM0dest[kMM0] = mfaab;
	  fPM0dest[k0M0] = mfcab;
	  fMP0dest[kM00] = mfacb;
	  fP0Pdest[k000] = mfcbc;
	  fM0Mdest[kM0M] = mfaba;
	  fP0Mdest[k00M] = mfcba;
	  fM0Pdest[kM00] = mfabc;
	  f0PPdest[k000] = mfbcc;
	  f0MMdest[k0MM] = mfbaa;
	  f0PMdest[k00M] = mfbca;
	  f0MPdest[k0M0] = mfbac;
	  f000dest[k000] = mfbbb;
	  fMMMdest[kMMM] = mfaaa;
	  fMMPdest[kMM0] = mfaac;
	  fMPPdest[kM00] = mfacc;
	  fMPMdest[kM0M] = mfaca;
	  fPPMdest[k00M] = mfcca;
	  fPPPdest[k000] = mfccc;
	  fPMPdest[k0M0] = mfcac;
	  fPMMdest[k0MM] = mfcaa;
	  ////////////////////////////////////////////////////////////////////////////////////
   }
}
//////////////////////////////////////////////////////////////////////////






















































