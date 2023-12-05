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
//! \author Martin Schoenherr
//=======================================================================================
#include "Calculation/Calculation.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__ void scaleCoarseToFineAdvectionDiffusion_Device(
    real* DC, 
    real* DF,
    real* DD27C, 
    real* DD27F, 
    uint* neighborCX,
    uint* neighborCY,
    uint* neighborCZ,
    uint* neighborFX,
    uint* neighborFY,
    uint* neighborFZ,
    unsigned long long numberOfLBnodesCoarse, 
    unsigned long long numberOfLBnodesFine, 
    bool isEvenTimestep,
    uint* posCSWB, 
    uint* posFSWB, 
    uint kCF, 
    real nu,
    real diffusivity_fine,
    ICellNeigh offCF)
{
   real *feF, *fwF, *fnF, *fsF, *ftF, *fbF, *fneF, *fswF, *fseF, *fnwF, *fteF, *fbwF, *fbeF, *ftwF, *ftnF, *fbsF, *fbnF, *ftsF, /**fzeroF,*/ *ftneF, *ftswF, *ftseF, *ftnwF, *fbneF, *fbswF, *fbseF, *fbnwF;

   feF    = &DF[dP00 * numberOfLBnodesFine];
   fwF    = &DF[dM00 * numberOfLBnodesFine];
   fnF    = &DF[d0P0 * numberOfLBnodesFine];
   fsF    = &DF[d0M0 * numberOfLBnodesFine];
   ftF    = &DF[d00P * numberOfLBnodesFine];
   fbF    = &DF[d00M * numberOfLBnodesFine];
   fneF   = &DF[dPP0 * numberOfLBnodesFine];
   fswF   = &DF[dMM0 * numberOfLBnodesFine];
   fseF   = &DF[dPM0 * numberOfLBnodesFine];
   fnwF   = &DF[dMP0 * numberOfLBnodesFine];
   fteF   = &DF[dP0P * numberOfLBnodesFine];
   fbwF   = &DF[dM0M * numberOfLBnodesFine];
   fbeF   = &DF[dP0M * numberOfLBnodesFine];
   ftwF   = &DF[dM0P * numberOfLBnodesFine];
   ftnF   = &DF[d0PP * numberOfLBnodesFine];
   fbsF   = &DF[d0MM * numberOfLBnodesFine];
   fbnF   = &DF[d0PM * numberOfLBnodesFine];
   ftsF   = &DF[d0MP * numberOfLBnodesFine];
   ftneF  = &DF[dPPP * numberOfLBnodesFine];
   ftswF  = &DF[dMMP * numberOfLBnodesFine];
   ftseF  = &DF[dPMP * numberOfLBnodesFine];
   ftnwF  = &DF[dMPP * numberOfLBnodesFine];
   fbneF  = &DF[dPPM * numberOfLBnodesFine];
   fbswF  = &DF[dMMM * numberOfLBnodesFine];
   fbseF  = &DF[dPMM * numberOfLBnodesFine];
   fbnwF  = &DF[dMPM * numberOfLBnodesFine];

   real *feC, *fwC, *fnC, *fsC, *ftC, *fbC, *fneC, *fswC, *fseC, *fnwC, *fteC, *fbwC, *fbeC, *ftwC, *ftnC, *fbsC, *fbnC, *ftsC, //*fzeroC,
      *ftneC, *ftswC, *ftseC, *ftnwC, *fbneC, *fbswC, *fbseC, *fbnwC;

   if (isEvenTimestep==true)
   {
      feC    = &DC[dP00 * numberOfLBnodesCoarse];
      fwC    = &DC[dM00 * numberOfLBnodesCoarse];
      fnC    = &DC[d0P0 * numberOfLBnodesCoarse];
      fsC    = &DC[d0M0 * numberOfLBnodesCoarse];
      ftC    = &DC[d00P * numberOfLBnodesCoarse];
      fbC    = &DC[d00M * numberOfLBnodesCoarse];
      fneC   = &DC[dPP0 * numberOfLBnodesCoarse];
      fswC   = &DC[dMM0 * numberOfLBnodesCoarse];
      fseC   = &DC[dPM0 * numberOfLBnodesCoarse];
      fnwC   = &DC[dMP0 * numberOfLBnodesCoarse];
      fteC   = &DC[dP0P * numberOfLBnodesCoarse];
      fbwC   = &DC[dM0M * numberOfLBnodesCoarse];
      fbeC   = &DC[dP0M * numberOfLBnodesCoarse];
      ftwC   = &DC[dM0P * numberOfLBnodesCoarse];
      ftnC   = &DC[d0PP * numberOfLBnodesCoarse];
      fbsC   = &DC[d0MM * numberOfLBnodesCoarse];
      fbnC   = &DC[d0PM * numberOfLBnodesCoarse];
      ftsC   = &DC[d0MP * numberOfLBnodesCoarse];
      ftneC  = &DC[dPPP * numberOfLBnodesCoarse];
      ftswC  = &DC[dMMP * numberOfLBnodesCoarse];
      ftseC  = &DC[dPMP * numberOfLBnodesCoarse];
      ftnwC  = &DC[dMPP * numberOfLBnodesCoarse];
      fbneC  = &DC[dPPM * numberOfLBnodesCoarse];
      fbswC  = &DC[dMMM * numberOfLBnodesCoarse];
      fbseC  = &DC[dPMM * numberOfLBnodesCoarse];
      fbnwC  = &DC[dMPM * numberOfLBnodesCoarse];
   } 
   else
   {
      fwC    = &DC[dP00 * numberOfLBnodesCoarse];
      feC    = &DC[dM00 * numberOfLBnodesCoarse];
      fsC    = &DC[d0P0 * numberOfLBnodesCoarse];
      fnC    = &DC[d0M0 * numberOfLBnodesCoarse];
      fbC    = &DC[d00P * numberOfLBnodesCoarse];
      ftC    = &DC[d00M * numberOfLBnodesCoarse];
      fswC   = &DC[dPP0 * numberOfLBnodesCoarse];
      fneC   = &DC[dMM0 * numberOfLBnodesCoarse];
      fnwC   = &DC[dPM0 * numberOfLBnodesCoarse];
      fseC   = &DC[dMP0 * numberOfLBnodesCoarse];
      fbwC   = &DC[dP0P * numberOfLBnodesCoarse];
      fteC   = &DC[dM0M * numberOfLBnodesCoarse];
      ftwC   = &DC[dP0M * numberOfLBnodesCoarse];
      fbeC   = &DC[dM0P * numberOfLBnodesCoarse];
      fbsC   = &DC[d0PP * numberOfLBnodesCoarse];
      ftnC   = &DC[d0MM * numberOfLBnodesCoarse];
      ftsC   = &DC[d0PM * numberOfLBnodesCoarse];
      fbnC   = &DC[d0MP * numberOfLBnodesCoarse];
      fbswC  = &DC[dPPP * numberOfLBnodesCoarse];
      fbneC  = &DC[dMMP * numberOfLBnodesCoarse];
      fbnwC  = &DC[dPMP * numberOfLBnodesCoarse];
      fbseC  = &DC[dMPP * numberOfLBnodesCoarse];
      ftswC  = &DC[dPPM * numberOfLBnodesCoarse];
      ftneC  = &DC[dMMM * numberOfLBnodesCoarse];
      ftnwC  = &DC[dPMM * numberOfLBnodesCoarse];
      ftseC  = &DC[dMPM * numberOfLBnodesCoarse];
   }


   Distributions27 D27F;
   D27F.f[dP00] = &DD27F[dP00 * numberOfLBnodesFine];
   D27F.f[dM00] = &DD27F[dM00 * numberOfLBnodesFine];
   D27F.f[d0P0] = &DD27F[d0P0 * numberOfLBnodesFine];
   D27F.f[d0M0] = &DD27F[d0M0 * numberOfLBnodesFine];
   D27F.f[d00P] = &DD27F[d00P * numberOfLBnodesFine];
   D27F.f[d00M] = &DD27F[d00M * numberOfLBnodesFine];
   D27F.f[dPP0] = &DD27F[dPP0 * numberOfLBnodesFine];
   D27F.f[dMM0] = &DD27F[dMM0 * numberOfLBnodesFine];
   D27F.f[dPM0] = &DD27F[dPM0 * numberOfLBnodesFine];
   D27F.f[dMP0] = &DD27F[dMP0 * numberOfLBnodesFine];
   D27F.f[dP0P] = &DD27F[dP0P * numberOfLBnodesFine];
   D27F.f[dM0M] = &DD27F[dM0M * numberOfLBnodesFine];
   D27F.f[dP0M] = &DD27F[dP0M * numberOfLBnodesFine];
   D27F.f[dM0P] = &DD27F[dM0P * numberOfLBnodesFine];
   D27F.f[d0PP] = &DD27F[d0PP * numberOfLBnodesFine];
   D27F.f[d0MM] = &DD27F[d0MM * numberOfLBnodesFine];
   D27F.f[d0PM] = &DD27F[d0PM * numberOfLBnodesFine];
   D27F.f[d0MP] = &DD27F[d0MP * numberOfLBnodesFine];
   D27F.f[d000] = &DD27F[d000 * numberOfLBnodesFine];
   D27F.f[dPPP] = &DD27F[dPPP * numberOfLBnodesFine];
   D27F.f[dMMP] = &DD27F[dMMP * numberOfLBnodesFine];
   D27F.f[dPMP] = &DD27F[dPMP * numberOfLBnodesFine];
   D27F.f[dMPP] = &DD27F[dMPP * numberOfLBnodesFine];
   D27F.f[dPPM] = &DD27F[dPPM * numberOfLBnodesFine];
   D27F.f[dMMM] = &DD27F[dMMM * numberOfLBnodesFine];
   D27F.f[dPMM] = &DD27F[dPMM * numberOfLBnodesFine];
   D27F.f[dMPM] = &DD27F[dMPM * numberOfLBnodesFine];

   Distributions27 D27C;
   if (isEvenTimestep==true)
   {
      D27C.f[dP00] = &DD27C[dP00 * numberOfLBnodesCoarse];
      D27C.f[dM00] = &DD27C[dM00 * numberOfLBnodesCoarse];
      D27C.f[d0P0] = &DD27C[d0P0 * numberOfLBnodesCoarse];
      D27C.f[d0M0] = &DD27C[d0M0 * numberOfLBnodesCoarse];
      D27C.f[d00P] = &DD27C[d00P * numberOfLBnodesCoarse];
      D27C.f[d00M] = &DD27C[d00M * numberOfLBnodesCoarse];
      D27C.f[dPP0] = &DD27C[dPP0 * numberOfLBnodesCoarse];
      D27C.f[dMM0] = &DD27C[dMM0 * numberOfLBnodesCoarse];
      D27C.f[dPM0] = &DD27C[dPM0 * numberOfLBnodesCoarse];
      D27C.f[dMP0] = &DD27C[dMP0 * numberOfLBnodesCoarse];
      D27C.f[dP0P] = &DD27C[dP0P * numberOfLBnodesCoarse];
      D27C.f[dM0M] = &DD27C[dM0M * numberOfLBnodesCoarse];
      D27C.f[dP0M] = &DD27C[dP0M * numberOfLBnodesCoarse];
      D27C.f[dM0P] = &DD27C[dM0P * numberOfLBnodesCoarse];
      D27C.f[d0PP] = &DD27C[d0PP * numberOfLBnodesCoarse];
      D27C.f[d0MM] = &DD27C[d0MM * numberOfLBnodesCoarse];
      D27C.f[d0PM] = &DD27C[d0PM * numberOfLBnodesCoarse];
      D27C.f[d0MP] = &DD27C[d0MP * numberOfLBnodesCoarse];
      D27C.f[d000] = &DD27C[d000 * numberOfLBnodesCoarse];
      D27C.f[dPPP] = &DD27C[dPPP * numberOfLBnodesCoarse];
      D27C.f[dMMP] = &DD27C[dMMP * numberOfLBnodesCoarse];
      D27C.f[dPMP] = &DD27C[dPMP * numberOfLBnodesCoarse];
      D27C.f[dMPP] = &DD27C[dMPP * numberOfLBnodesCoarse];
      D27C.f[dPPM] = &DD27C[dPPM * numberOfLBnodesCoarse];
      D27C.f[dMMM] = &DD27C[dMMM * numberOfLBnodesCoarse];
      D27C.f[dPMM] = &DD27C[dPMM * numberOfLBnodesCoarse];
      D27C.f[dMPM] = &DD27C[dMPM * numberOfLBnodesCoarse];
   }
   else
   {
      D27C.f[dM00] = &DD27C[dP00 * numberOfLBnodesCoarse];
      D27C.f[dP00] = &DD27C[dM00 * numberOfLBnodesCoarse];
      D27C.f[d0M0] = &DD27C[d0P0 * numberOfLBnodesCoarse];
      D27C.f[d0P0] = &DD27C[d0M0 * numberOfLBnodesCoarse];
      D27C.f[d00M] = &DD27C[d00P * numberOfLBnodesCoarse];
      D27C.f[d00P] = &DD27C[d00M * numberOfLBnodesCoarse];
      D27C.f[dMM0] = &DD27C[dPP0 * numberOfLBnodesCoarse];
      D27C.f[dPP0] = &DD27C[dMM0 * numberOfLBnodesCoarse];
      D27C.f[dMP0] = &DD27C[dPM0 * numberOfLBnodesCoarse];
      D27C.f[dPM0] = &DD27C[dMP0 * numberOfLBnodesCoarse];
      D27C.f[dM0M] = &DD27C[dP0P * numberOfLBnodesCoarse];
      D27C.f[dP0P] = &DD27C[dM0M * numberOfLBnodesCoarse];
      D27C.f[dM0P] = &DD27C[dP0M * numberOfLBnodesCoarse];
      D27C.f[dP0M] = &DD27C[dM0P * numberOfLBnodesCoarse];
      D27C.f[d0MM] = &DD27C[d0PP * numberOfLBnodesCoarse];
      D27C.f[d0PP] = &DD27C[d0MM * numberOfLBnodesCoarse];
      D27C.f[d0MP] = &DD27C[d0PM * numberOfLBnodesCoarse];
      D27C.f[d0PM] = &DD27C[d0MP * numberOfLBnodesCoarse];
      D27C.f[d000] = &DD27C[d000 * numberOfLBnodesCoarse];
      D27C.f[dMMM] = &DD27C[dPPP * numberOfLBnodesCoarse];
      D27C.f[dPPM] = &DD27C[dMMP * numberOfLBnodesCoarse];
      D27C.f[dMPM] = &DD27C[dPMP * numberOfLBnodesCoarse];
      D27C.f[dPMM] = &DD27C[dMPP * numberOfLBnodesCoarse];
      D27C.f[dMMP] = &DD27C[dPPM * numberOfLBnodesCoarse];
      D27C.f[dPPP] = &DD27C[dMMM * numberOfLBnodesCoarse];
      D27C.f[dMPP] = &DD27C[dPMM * numberOfLBnodesCoarse];
      D27C.f[dPMP] = &DD27C[dMPM * numberOfLBnodesCoarse];
   }

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  ix = threadIdx.x;
   const unsigned  iy = blockIdx.x; 
   const unsigned  iz = blockIdx.y; 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*iz + iy) + ix;
   //////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   real vx1,vx2,vx3, cu_sq;
   real f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_TNE,f_TSW,f_TSE,f_TNW,f_BNE,f_BSW,f_BSE,f_BNW;

   real f27E,f27W,f27N,f27S,f27T,f27B,f27NE,f27SW,f27SE,f27NW,f27TE,f27BW,f27BE,f27TW,f27TN,f27BS,f27BN,f27TS,f27ZERO,f27TNE,f27TSW,f27TSE,f27TNW,f27BNE,f27BSW,f27BSE,f27BNW;
   real Mx,My,Mz; 
   real Conc_C_SWB, Conc_C_SWT, Conc_C_SET, Conc_C_SEB, Conc_C_NWB, Conc_C_NWT, Conc_C_NET, Conc_C_NEB;

   real omegaD_C     = c2o1 / (c6o1 * diffusivity_fine/c2o1 + c1o1);
   real omegaD_F     = c2o1 / (c6o1 * diffusivity_fine + c1o1);

   real x,       y,       z;
   real xoff,    yoff,    zoff;
   real xoff_sq, yoff_sq, zoff_sq;

   if(k<kCF)
   {
      //////////////////////////////////////////////////////////////////////////
      xoff    = offCF.x[k];
      yoff    = offCF.y[k];
      zoff    = offCF.z[k];
      xoff_sq = xoff * xoff;
      yoff_sq = yoff * yoff;
      zoff_sq = zoff * zoff;
      //////////////////////////////////////////////////////////////////////////
      //SWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      unsigned int k0zero= posCSWB[k];
      unsigned int k0w   = neighborCX[k0zero];
      unsigned int k0s   = neighborCY[k0zero];
      unsigned int k0b   = neighborCZ[k0zero];
      unsigned int k0sw  = neighborCY[k0w];
      unsigned int k0bw  = neighborCZ[k0w];
      unsigned int k0bs  = neighborCZ[k0s];
      unsigned int k0bsw = neighborCZ[k0sw];
      //////////////////////////////////////////////////////////////////////////
      //index 
      unsigned int kzero= k0zero;
      unsigned int kw   = k0w;   
      unsigned int ks   = k0s;   
      unsigned int kb   = k0b;   
      unsigned int ksw  = k0sw;  
      unsigned int kbw  = k0bw;  
      unsigned int kbs  = k0bs;  
      unsigned int kbsw = k0bsw; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_SWB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CSWB = (Conc_C_SWB * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CSWB = (Conc_C_SWB * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CSWB = (Conc_C_SWB * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //SWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborCZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborCZ[kbw];  
      kbs  = neighborCZ[kbs];  
      kbsw = neighborCZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_SWT = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CSWT = (Conc_C_SWT * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CSWT = (Conc_C_SWT * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CSWT = (Conc_C_SWT * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //SET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborCX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborCX[ksw];  
      kbw  = neighborCX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborCX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_SET = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CSET = (Conc_C_SET * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CSET = (Conc_C_SET * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CSET = (Conc_C_SET * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //SEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborCX[k0w];   
      ks   = k0sw;   
      ksw  = neighborCX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_SEB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CSEB = (Conc_C_SEB * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CSEB = (Conc_C_SEB * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CSEB = (Conc_C_SEB * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //NWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      k0zero= k0s;
      k0w   = k0sw;
      k0s   = neighborCY[k0s];
      k0b   = k0bs;
      k0sw  = neighborCY[k0sw];
      k0bw  = k0bsw;
      k0bs  = neighborCY[k0bs];
      k0bsw = neighborCY[k0bsw];
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= k0zero;
      kw   = k0w;   
      ks   = k0s;   
      kb   = k0b;   
      ksw  = k0sw;  
      kbw  = k0bw;  
      kbs  = k0bs;  
      kbsw = k0bsw; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_NWB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CNWB = (Conc_C_NWB * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CNWB = (Conc_C_NWB * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CNWB = (Conc_C_NWB * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //NWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborCZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborCZ[kbw];  
      kbs  = neighborCZ[kbs];  
      kbsw = neighborCZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_NWT = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CNWT = (Conc_C_NWT * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CNWT = (Conc_C_NWT * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CNWT = (Conc_C_NWT * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //NET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborCX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborCX[ksw];  
      kbw  = neighborCX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborCX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_NET = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CNET = (Conc_C_NET * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CNET = (Conc_C_NET * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CNET = (Conc_C_NET * vx3 - Mz) * (c3o1*omegaD_C);




      //////////////////////////////////////////////////////////////////////////
      //NEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborCX[k0w];   
      ks   = k0sw;   
      ksw  = neighborCX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feC[kzero];
      f_W    = fwC[kw];
      f_N    = fnC[kzero];
      f_S    = fsC[ks];
      f_T    = ftC[kzero];
      f_B    = fbC[kb];
      f_NE   = fneC[kzero];
      f_SW   = fswC[ksw];
      f_SE   = fseC[ks];
      f_NW   = fnwC[kw];
      f_TE   = fteC[kzero];
      f_BW   = fbwC[kbw];
      f_BE   = fbeC[kb];
      f_TW   = ftwC[kw];
      f_TN   = ftnC[kzero];
      f_BS   = fbsC[kbs];
      f_BN   = fbnC[kb];
      f_TS   = ftsC[ks];
      f_TNE  = ftneC[kzero];
      f_TSW  = ftswC[ksw];
      f_TSE  = ftseC[ks];
      f_TNW  = ftnwC[kw];
      f_BNE  = fbneC[kb];
      f_BSW  = fbswC[kbsw];
      f_BSE  = fbseC[kbs];
      f_BNW  = fbnwC[kbw];
      ////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27C.f[dP00])[kzero];//ke
      f27W    =  (D27C.f[dM00])[kw   ];
      f27N    =  (D27C.f[d0P0])[kzero];//kn
      f27S    =  (D27C.f[d0M0])[ks   ];
      f27T    =  (D27C.f[d00P])[kzero];//kt
      f27B    =  (D27C.f[d00M])[kb   ];
      f27NE   =  (D27C.f[dPP0])[kzero];//kne
      f27SW   =  (D27C.f[dMM0])[ksw  ];
      f27SE   =  (D27C.f[dPM0])[ks   ];//kse
      f27NW   =  (D27C.f[dMP0])[kw   ];//knw
      f27TE   =  (D27C.f[dP0P])[kzero];//kte
      f27BW   =  (D27C.f[dM0M])[kbw  ];
      f27BE   =  (D27C.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27C.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27C.f[d0PP])[kzero];//ktn
      f27BS   =  (D27C.f[d0MM])[kbs  ];
      f27BN   =  (D27C.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27C.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27C.f[d000])[kzero];//kzero
      f27TNE   = (D27C.f[dPPP])[kzero];//ktne
      f27TSW   = (D27C.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27C.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27C.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27C.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27C.f[dMMM])[kbsw ];
      f27BSE   = (D27C.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27C.f[dMPM])[kbw  ];//kbnw

      Conc_C_NEB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW; 

      vx1  =f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  =f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  =f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_CNEB = (Conc_C_NEB * vx1 - Mx) * (c3o1*omegaD_C);
      real Diff_Conc_Y_CNEB = (Conc_C_NEB * vx2 - My) * (c3o1*omegaD_C);
      real Diff_Conc_Z_CNEB = (Conc_C_NEB * vx3 - Mz) * (c3o1*omegaD_C);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //linear
      //real Diff_Conc_XX = zero;
      //real Diff_Conc_YY = zero;
      //real Diff_Conc_ZZ = zero;
      //quadratic
      real Diff_Conc_XX = ((Diff_Conc_X_CNEB + Diff_Conc_X_CSEB + Diff_Conc_X_CNET + Diff_Conc_X_CSET) - (Diff_Conc_X_CNWB + Diff_Conc_X_CSWB + Diff_Conc_X_CNWT + Diff_Conc_X_CSWT)) * c1o4;
      real Diff_Conc_YY = ((Diff_Conc_Y_CNEB + Diff_Conc_Y_CNWB + Diff_Conc_Y_CNET + Diff_Conc_Y_CNWT) - (Diff_Conc_Y_CSEB + Diff_Conc_Y_CSWB + Diff_Conc_Y_CSET + Diff_Conc_Y_CSWT)) * c1o4;
      real Diff_Conc_ZZ = ((Diff_Conc_Z_CSET + Diff_Conc_Z_CSWT + Diff_Conc_Z_CNET + Diff_Conc_Z_CNWT) - (Diff_Conc_Z_CSEB + Diff_Conc_Z_CSWB + Diff_Conc_Z_CNEB + Diff_Conc_Z_CNWB)) * c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      real dx = c1o4 * (Conc_C_NEB + Conc_C_NET - Conc_C_NWB - Conc_C_NWT + Conc_C_SEB + Conc_C_SET - Conc_C_SWB - Conc_C_SWT);
      real dy = c1o4 * (Conc_C_NEB + Conc_C_NET + Conc_C_NWB + Conc_C_NWT - Conc_C_SEB - Conc_C_SET - Conc_C_SWB - Conc_C_SWT);
      real dz = c1o4 * (-Conc_C_NEB + Conc_C_NET - Conc_C_NWB + Conc_C_NWT - Conc_C_SEB + Conc_C_SET - Conc_C_SWB + Conc_C_SWT);
      real dxx = Diff_Conc_XX * c1o2;
      real dyy = Diff_Conc_YY * c1o2;
      real dzz = Diff_Conc_ZZ * c1o2;
      real dxy = c1o2 * (Conc_C_NEB + Conc_C_NET - Conc_C_NWB - Conc_C_NWT - Conc_C_SEB - Conc_C_SET + Conc_C_SWB + Conc_C_SWT);
      real dyz = c1o2 * (-Conc_C_NEB + Conc_C_NET - Conc_C_NWB + Conc_C_NWT + Conc_C_SEB - Conc_C_SET + Conc_C_SWB - Conc_C_SWT);
      real dxz = c1o2 * (-Conc_C_NEB + Conc_C_NET + Conc_C_NWB - Conc_C_NWT - Conc_C_SEB + Conc_C_SET + Conc_C_SWB - Conc_C_SWT);
      real dxyz = -Conc_C_NEB + Conc_C_NET + Conc_C_NWB - Conc_C_NWT + Conc_C_SEB - Conc_C_SET - Conc_C_SWB + Conc_C_SWT;
      real d0 = c1o8 * (-c2o1 * dxx - c2o1 * dyy - c2o1 * dzz + Conc_C_NEB + Conc_C_NET + Conc_C_NWB + Conc_C_NWT + Conc_C_SEB + Conc_C_SET + Conc_C_SWB + Conc_C_SWT);

     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // offset for refinement into the wall
      //
      // X------X
      // |      | x---x    
      // |   ---+-+-> |    ----> off-vector
      // |      | x---x 
      // X------X   
      //            
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      d0 = d0 + xoff * dx + yoff * dy + zoff * dz + xoff_sq * dxx + yoff_sq * dyy + zoff_sq * dzz + xoff*yoff*dxy + xoff*zoff*dxz + yoff*zoff*dyz;
      dx = dx + c2o1 * xoff * dxx + yoff * dxy + zoff * dxz;
      dy = dy + c2o1 * yoff * dyy + xoff * dxy + zoff * dyz;
      dz = dz + c2o1 * zoff * dzz + xoff * dxz + yoff * dyz;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position SWB -0.25, -0.25, -0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x = -c1o4;
      y = -c1o4;
      z = -c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 0
      k0zero= posFSWB[k];
      k0w   = neighborFX[k0zero];
      k0s   = neighborFY[k0zero];
      k0b   = neighborFZ[k0zero];
      k0sw  = neighborFY[k0w];
      k0bw  = neighborFZ[k0w];
      k0bs  = neighborFZ[k0s];
      k0bsw = neighborFZ[k0sw];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kzero= k0zero;
      kw   = k0w;   
      ks   = k0s;   
      kb   = k0b;   
      ksw  = k0sw;  
      kbw  = k0bw;  
      kbs  = k0bs;  
      kbsw = k0bsw; 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      real Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      real Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      real Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      real Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position SWT -0.25, -0.25, 0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x = -c1o4;
      y = -c1o4;
      z =  c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborFZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborFZ[kbw];  
      kbs  = neighborFZ[kbs];  
      kbsw = neighborFZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position SET 0.25, -0.25, 0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x =  c1o4;
      y = -c1o4;
      z =  c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborFX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborFX[ksw];  
      kbw  = neighborFX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborFX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position SEB 0.25, -0.25, -0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x =  c1o4;
      y = -c1o4;
      z = -c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborFX[k0w];   
      ks   = k0sw;   
      ksw  = neighborFX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position NWB -0.25, 0.25, -0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x = -c1o4;
      y =  c1o4;
      z = -c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 0
      k0zero= k0s;
      k0w   = k0sw;
      k0s   = neighborFY[k0s];
      k0b   = k0bs;
      k0sw  = neighborFY[k0sw];
      k0bw  = k0bsw;
      k0bs  = neighborFY[k0bs];
      k0bsw = neighborFY[k0bsw];
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= k0zero;
      kw   = k0w;   
      ks   = k0s;   
      kb   = k0b;   
      ksw  = k0sw;  
      kbw  = k0bw;  
      kbs  = k0bs;  
      kbsw = k0bsw; 
      ////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position NWT -0.25, 0.25, 0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x = -c1o4;
      y =  c1o4;
      z =  c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborFZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborFZ[kbw];  
      kbs  = neighborFZ[kbs];  
      kbsw = neighborFZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position NET 0.25, 0.25, 0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x = c1o4;
      y = c1o4;
      z = c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborFX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborFX[ksw];  
      kbw  = neighborFX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborFX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));








      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position NEB 0.25, 0.25, -0.25
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      x =  c1o4;
      y =  c1o4;
      z = -c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborFX[k0w];   
      ks   = k0sw;   
      ksw  = neighborFX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vx1=feF[kzero]+fneF[kzero]+fseF[ks]+fteF[kzero]+fbeF[kb]-fwF[kw]-fnwF[kw]-fswF[ksw]-ftwF[kw]-fbwF[kbw]+ftneF[kzero]-ftswF[ksw]+ftseF[ks]-ftnwF[kw]+fbneF[kb]-fbswF[kbsw]+fbseF[kbs]-fbnwF[kbw];
      vx2=fnF[kzero]+fneF[kzero]+fnwF[kw]+ftnF[kzero]+fbnF[kb]-fsF[ks]-fseF[ks]-fswF[ksw]-ftsF[ks]-fbsF[kbs]+ftneF[kzero]-ftswF[ksw]-ftseF[ks]+ftnwF[kw]+fbneF[kb]-fbswF[kbsw]-fbseF[kbs]+fbnwF[kbw];
      vx3=ftF[kzero]+fteF[kzero]+ftwF[kw]+ftnF[kzero]+ftsF[ks]-fbF[kb]-fbeF[kb]-fbwF[kbw]-fbnF[kb]-fbsF[kbs]+ftneF[kzero]+ftswF[ksw]+ftseF[ks]+ftnwF[kw]-fbneF[kb]-fbswF[kbsw]-fbseF[kbs]-fbnwF[kbw];

      Conc_F = d0 + dx*x + dy*y + dz*z + dxx*x*x + dyy*y*y + dzz*z*z + dxy*x*y +  dxz*x*z + dyz*y*z + dxyz*x*y*z;

      Diff_Conc_X = dx + x * dxx + y * dxy + z * dxz + y * z * dxyz;
      Diff_Conc_Y = dy + y * dyy + x * dxy + z * dyz + x * z * dxyz;
      Diff_Conc_Z = dz + z * dzz + x * dxz + y * dyz + x * y * dxyz;

      Mx = Conc_F*vx1-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_X;
      My = Conc_F*vx2-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Y;
      Mz = Conc_F*vx3-(c1o1)/(c3o1*omegaD_F)*c1o2*Diff_Conc_Z;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27F.f[d000])[kzero] =   c8o27* Conc_F*(c1o1-cu_sq);
      (D27F.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_F*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27F.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_F*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27F.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_F*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27F.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_F*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27F.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_F*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27F.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_F*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27F.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_F*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27F.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_F*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27F.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_F*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27F.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_F*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27F.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_F*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27F.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_F*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27F.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_F*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27F.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_F*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27F.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_F*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27F.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_F*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27F.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_F*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27F.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_F*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27F.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27F.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27F.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_F*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27F.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_F*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27F.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27F.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27F.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_F*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27F.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_F*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));
   }
}
