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

__global__ void scaleFineToCoarseAdvectionDiffusion_Device(
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
    uint* posC, 
    uint* posFSWB, 
    uint kFC, 
    real nu,
    real diffusivity_coarse,
    ICellNeigh offFC)
{
   real *feF, *fwF, *fnF, *fsF, *ftF, *fbF, 
       *fneF, *fswF, *fseF, *fnwF, *fteF, *fbwF, 
       *fbeF, *ftwF, *ftnF, *fbsF, *fbnF, *ftsF,
       *ftneF, *ftswF, *ftseF, *ftnwF, 
       *fbneF, *fbswF, *fbseF, *fbnwF;

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

   real *feC, *fwC, *fnC, *fsC, *ftC, *fbC, 
       *fneC, *fswC, *fseC, *fnwC, *fteC, *fbwC, 
       *fbeC, *ftwC, *ftnC, *fbsC, *fbnC, *ftsC, 
       *ftneC, *ftswC, *ftseC, *ftnwC, 
       *fbneC, *fbswC, *fbseC, *fbnwC;

   if (isEvenTimestep)
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
   if (isEvenTimestep)
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

   ////////////////////////////////////////////////////////////////////////////////
   real vx1,vx2,vx3,cu_sq;
   real f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_TNE,f_TSW,f_TSE,f_TNW,f_BNE,f_BSW,f_BSE,f_BNW;

   real f27E,f27W,f27N,f27S,f27T,f27B,f27NE,f27SW,f27SE,f27NW,f27TE,f27BW,f27BE,f27TW,f27TN,f27BS,f27BN,f27TS,f27ZERO,f27TNE,f27TSW,f27TSE,f27TNW,f27BNE,f27BSW,f27BSE,f27BNW;
   real Mx,My,Mz; 
   real Conc_F_SWB, Conc_F_SWT, Conc_F_SET, Conc_F_SEB, Conc_F_NWB, Conc_F_NWT, Conc_F_NET, Conc_F_NEB;

   real omegaD_C = c2o1 / (c6o1 * diffusivity_coarse + c1o1);
   real omegaD_F = c2o1 / (c6o1 * diffusivity_coarse*c2o1 + c1o1);

   real xoff,    yoff,    zoff;
   real xoff_sq, yoff_sq, zoff_sq;

   if(k<kFC){
      //////////////////////////////////////////////////////////////////////////
      xoff    = offFC.x[k];
      yoff    = offFC.y[k];
      zoff    = offFC.z[k];
      xoff_sq = xoff * xoff;
      yoff_sq = yoff * yoff;
      zoff_sq = zoff * zoff;
      //////////////////////////////////////////////////////////////////////////
      //SWB//
      //////////////////////////////////////////////////////////////////////////
      //index 0
      unsigned int k0zero= posFSWB[k];
      unsigned int k0w   = neighborFX[k0zero];
      unsigned int k0s   = neighborFY[k0zero];
      unsigned int k0b   = neighborFZ[k0zero];
      unsigned int k0sw  = neighborFY[k0w];
      unsigned int k0bw  = neighborFZ[k0w];
      unsigned int k0bs  = neighborFZ[k0s];
      unsigned int k0bsw = neighborFZ[k0sw];
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
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_SWB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FSWB = (Conc_F_SWB * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FSWB = (Conc_F_SWB * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FSWB = (Conc_F_SWB * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //SWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborFZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborFZ[kbw];  
      kbs  = neighborFZ[kbs];  
      kbsw = neighborFZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_SWT = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FSWT = (Conc_F_SWT * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FSWT = (Conc_F_SWT * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FSWT = (Conc_F_SWT * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //SET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborFX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborFX[ksw];  
      kbw  = neighborFX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborFX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_SET = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FSET = (Conc_F_SET * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FSET = (Conc_F_SET * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FSET = (Conc_F_SET * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //SEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborFX[k0w];   
      ks   = k0sw;   
      ksw  = neighborFX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_SEB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FSEB = (Conc_F_SEB * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FSEB = (Conc_F_SEB * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FSEB = (Conc_F_SEB * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //NWB//
      //////////////////////////////////////////////////////////////////////////
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
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_NWB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FNWB = (Conc_F_NWB * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FNWB = (Conc_F_NWB * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FNWB = (Conc_F_NWB * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //NWT//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kb;
      kw   = kbw;   
      ks   = kbs;   
      kb   = neighborFZ[kb];   
      ksw  = kbsw;  
      kbw  = neighborFZ[kbw];  
      kbs  = neighborFZ[kbs];  
      kbsw = neighborFZ[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_NWT = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FNWT = (Conc_F_NWT * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FNWT = (Conc_F_NWT * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FNWT = (Conc_F_NWT * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //NET//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kzero= kw;
      kw   = neighborFX[kw];   
      ks   = ksw;   
      kb   = kbw;   
      ksw  = neighborFX[ksw];  
      kbw  = neighborFX[kbw];  
      kbs  = kbsw;  
      kbsw = neighborFX[kbsw]; 
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_NET = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FNET = (Conc_F_NET * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FNET = (Conc_F_NET * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FNET = (Conc_F_NET * vx3 - Mz) * (c3o1*omegaD_F);




      //////////////////////////////////////////////////////////////////////////
      //NEB//
      //////////////////////////////////////////////////////////////////////////
      //index 
      kb   = kzero;   
      kbw  = kw;  
      kbs  = ks;  
      kbsw = ksw; 
      kzero= k0w;
      kw   = neighborFX[k0w];   
      ks   = k0sw;   
      ksw  = neighborFX[k0sw];  
      ////////////////////////////////////////////////////////////////////////////////
      f_E    = feF[kzero];
      f_W    = fwF[kw];
      f_N    = fnF[kzero];
      f_S    = fsF[ks];
      f_T    = ftF[kzero];
      f_B    = fbF[kb];
      f_NE   = fneF[kzero];
      f_SW   = fswF[ksw];
      f_SE   = fseF[ks];
      f_NW   = fnwF[kw];
      f_TE   = fteF[kzero];
      f_BW   = fbwF[kbw];
      f_BE   = fbeF[kb];
      f_TW   = ftwF[kw];
      f_TN   = ftnF[kzero];
      f_BS   = fbsF[kbs];
      f_BN   = fbnF[kb];
      f_TS   = ftsF[ks];
      f_TNE  = ftneF[kzero];
      f_TSW  = ftswF[ksw];
      f_TSE  = ftseF[ks];
      f_TNW  = ftnwF[kw];
      f_BNE  = fbneF[kb];
      f_BSW  = fbswF[kbsw];
      f_BSE  = fbseF[kbs];
      f_BNW  = fbnwF[kbw];
      //////////////////////////////////////////////////////////////////////////////////
      f27E    =  (D27F.f[dP00])[kzero];//ke
      f27W    =  (D27F.f[dM00])[kw   ];
      f27N    =  (D27F.f[d0P0])[kzero];//kn
      f27S    =  (D27F.f[d0M0])[ks   ];
      f27T    =  (D27F.f[d00P])[kzero];//kt
      f27B    =  (D27F.f[d00M])[kb   ];
      f27NE   =  (D27F.f[dPP0])[kzero];//kne
      f27SW   =  (D27F.f[dMM0])[ksw  ];
      f27SE   =  (D27F.f[dPM0])[ks   ];//kse
      f27NW   =  (D27F.f[dMP0])[kw   ];//knw
      f27TE   =  (D27F.f[dP0P])[kzero];//kte
      f27BW   =  (D27F.f[dM0M])[kbw  ];
      f27BE   =  (D27F.f[dP0M])[kb   ];//kbe
      f27TW   =  (D27F.f[dM0P])[kw   ];//ktw
      f27TN   =  (D27F.f[d0PP])[kzero];//ktn
      f27BS   =  (D27F.f[d0MM])[kbs  ];
      f27BN   =  (D27F.f[d0PM])[kb   ];//kbn
      f27TS   =  (D27F.f[d0MP])[ks   ];//kts
      f27ZERO =  (D27F.f[d000])[kzero];//kzero
      f27TNE   = (D27F.f[dPPP])[kzero];//ktne
      f27TSW   = (D27F.f[dMMP])[ksw  ];//ktsw
      f27TSE   = (D27F.f[dPMP])[ks   ];//ktse
      f27TNW   = (D27F.f[dMPP])[kw   ];//ktnw
      f27BNE   = (D27F.f[dPPM])[kb   ];//kbne
      f27BSW   = (D27F.f[dMMM])[kbsw ];
      f27BSE   = (D27F.f[dPMM])[kbs  ];//kbse
      f27BNW   = (D27F.f[dMPM])[kbw  ];//kbnw

      Conc_F_NEB = f27E + f27W + f27N + f27S + f27T + f27B + f27NE + f27SW + f27SE + f27NW + 
                   f27TE + f27BW + f27BE + f27TW + f27TN + f27BS + f27BN + f27TS + f27ZERO + 
                   f27TNE + f27TSW + f27TSE + f27TNW + f27BNE + f27BSW + f27BSE + f27BNW;

      vx1  = f_E+f_NE+f_SE+f_TE+f_BE-f_W-f_NW-f_SW-f_TW-f_BW+f_TNE-f_TSW+f_TSE-f_TNW+f_BNE-f_BSW+f_BSE-f_BNW;
      vx2  = f_N+f_NE+f_NW+f_TN+f_BN-f_S-f_SE-f_SW-f_TS-f_BS+f_TNE-f_TSW-f_TSE+f_TNW+f_BNE-f_BSW-f_BSE+f_BNW;
      vx3  = f_T+f_TE+f_TW+f_TN+f_TS-f_B-f_BE-f_BW-f_BN-f_BS+f_TNE+f_TSW+f_TSE+f_TNW-f_BNE-f_BSW-f_BSE-f_BNW;
      Mx   =f27E+f27NE+f27SE+f27TE+f27BE-f27W-f27NW-f27SW-f27TW-f27BW+f27TNE-f27TSW+f27TSE-f27TNW+f27BNE-f27BSW+f27BSE-f27BNW;
      My   =f27N+f27NE+f27NW+f27TN+f27BN-f27S-f27SE-f27SW-f27TS-f27BS+f27TNE-f27TSW-f27TSE+f27TNW+f27BNE-f27BSW-f27BSE+f27BNW;
      Mz   =f27T+f27TE+f27TW+f27TN+f27TS-f27B-f27BE-f27BW-f27BN-f27BS+f27TNE+f27TSW+f27TSE+f27TNW-f27BNE-f27BSW-f27BSE-f27BNW;

      real Diff_Conc_X_FNEB = (Conc_F_NEB * vx1 - Mx) * (c3o1*omegaD_F);
      real Diff_Conc_Y_FNEB = (Conc_F_NEB * vx2 - My) * (c3o1*omegaD_F);
      real Diff_Conc_Z_FNEB = (Conc_F_NEB * vx3 - Mz) * (c3o1*omegaD_F);




      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //linear
      //real Diff_Conc_XX = zero;
      //real Diff_Conc_YY = zero;
      //real Diff_Conc_ZZ = zero;
      //quadratic
      real Diff_Conc_XX = ((Diff_Conc_X_FNEB + Diff_Conc_X_FSEB + Diff_Conc_X_FNET + Diff_Conc_X_FSET) - (Diff_Conc_X_FNWB + Diff_Conc_X_FSWB + Diff_Conc_X_FNWT + Diff_Conc_X_FSWT)) * c1o4;
      real Diff_Conc_YY = ((Diff_Conc_Y_FNEB + Diff_Conc_Y_FNWB + Diff_Conc_Y_FNET + Diff_Conc_Y_FNWT) - (Diff_Conc_Y_FSEB + Diff_Conc_Y_FSWB + Diff_Conc_Y_FSET + Diff_Conc_Y_FSWT)) * c1o4;
      real Diff_Conc_ZZ = ((Diff_Conc_Z_FSET + Diff_Conc_Z_FSWT + Diff_Conc_Z_FNET + Diff_Conc_Z_FNWT) - (Diff_Conc_Z_FSEB + Diff_Conc_Z_FSWB + Diff_Conc_Z_FNEB + Diff_Conc_Z_FNWB)) * c1o4;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      real dx = c1o4 * (Conc_F_NEB + Conc_F_NET - Conc_F_NWB - Conc_F_NWT + Conc_F_SEB + Conc_F_SET - Conc_F_SWB - Conc_F_SWT);
      real dy = c1o4 * (Conc_F_NEB + Conc_F_NET + Conc_F_NWB + Conc_F_NWT - Conc_F_SEB - Conc_F_SET - Conc_F_SWB - Conc_F_SWT);
      real dz = c1o4 * (-Conc_F_NEB + Conc_F_NET - Conc_F_NWB + Conc_F_NWT - Conc_F_SEB + Conc_F_SET - Conc_F_SWB + Conc_F_SWT);
      real dxx = Diff_Conc_XX * c1o2;
      real dyy = Diff_Conc_YY * c1o2;
      real dzz = Diff_Conc_ZZ * c1o2;
      real dxy = c1o2 * (Conc_F_NEB + Conc_F_NET - Conc_F_NWB - Conc_F_NWT - Conc_F_SEB - Conc_F_SET + Conc_F_SWB + Conc_F_SWT);
      real dyz = c1o2 * (-Conc_F_NEB + Conc_F_NET - Conc_F_NWB + Conc_F_NWT + Conc_F_SEB - Conc_F_SET + Conc_F_SWB - Conc_F_SWT);
      real dxz = c1o2 * (-Conc_F_NEB + Conc_F_NET + Conc_F_NWB - Conc_F_NWT - Conc_F_SEB + Conc_F_SET + Conc_F_SWB - Conc_F_SWT);

      real d0 = c1o8 * (-c2o1 * dxx - c2o1 * dyy - c2o1 * dzz + Conc_F_NEB + Conc_F_NET + Conc_F_NWB + Conc_F_NWT + Conc_F_SEB + Conc_F_SET + Conc_F_SWB + Conc_F_SWT);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // offset for refinement into the wall
      //
      // X------X
      // |      |
      // |   ---+-->X     ----> off-vector
      // |      |
      // X------X
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      d0 = d0 + xoff * dx + yoff * dy + zoff * dz + xoff_sq * dxx + yoff_sq * dyy + zoff_sq * dzz + xoff*yoff*dxy + xoff*zoff*dxz + yoff*zoff*dyz;
      dx = dx + c2o1 * xoff * dxx + yoff * dxy + zoff * dxz;
      dy = dy + c2o1 * yoff * dyy + xoff * dxy + zoff * dyz;
      dz = dz + c2o1 * zoff * dzz + xoff * dxz + yoff * dyz;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Position C 0.5, 0.5, 0.5
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //x = zero;
      //y = zero;
      //z = zero;

      //////////////////////////////////////////////////////////////////////////
      //index 0
      kzero= posC[k];
      kw   = neighborCX[kzero];
      ks   = neighborCY[kzero];
      kb   = neighborCZ[kzero];
      ksw  = neighborCY[kw];
      kbw  = neighborCZ[kw];
      kbs  = neighborCZ[ks];
      kbsw = neighborCZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      vx1=feC[kzero]+fneC[kzero]+fseC[ks]+fteC[kzero]+fbeC[kb]-fwC[kw]-fnwC[kw]-fswC[ksw]-ftwC[kw]-fbwC[kbw]+ftneC[kzero]-ftswC[ksw]+ftseC[ks]-ftnwC[kw]+fbneC[kb]-fbswC[kbsw]+fbseC[kbs]-fbnwC[kbw];
      vx2=fnC[kzero]+fneC[kzero]+fnwC[kw]+ftnC[kzero]+fbnC[kb]-fsC[ks]-fseC[ks]-fswC[ksw]-ftsC[ks]-fbsC[kbs]+ftneC[kzero]-ftswC[ksw]-ftseC[ks]+ftnwC[kw]+fbneC[kb]-fbswC[kbsw]-fbseC[kbs]+fbnwC[kbw];
      vx3=ftC[kzero]+fteC[kzero]+ftwC[kw]+ftnC[kzero]+ftsC[ks]-fbC[kb]-fbeC[kb]-fbwC[kbw]-fbnC[kb]-fbsC[kbs]+ftneC[kzero]+ftswC[ksw]+ftseC[ks]+ftnwC[kw]-fbneC[kb]-fbswC[kbsw]-fbseC[kbs]-fbnwC[kbw];

      real Conc_C = d0;

      real Diff_Conc_X_C = dx;
      real Diff_Conc_Y_C = dy;
      real Diff_Conc_Z_C = dz;

      Mx = Conc_C*vx1-(c1o1)/(c3o1*omegaD_C)*c2o1*Diff_Conc_X_C;
      My = Conc_C*vx2-(c1o1)/(c3o1*omegaD_C)*c2o1*Diff_Conc_Y_C;
      Mz = Conc_C*vx3-(c1o1)/(c3o1*omegaD_C)*c2o1*Diff_Conc_Z_C;

      cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      (D27C.f[d000])[kzero] =   c8o27* Conc_C*(c1o1-cu_sq);
      (D27C.f[dP00])[kzero] =   c2o27* (c3o1*( Mx        )+Conc_C*(c1o1+c9o2*( vx1        )*( vx1        )-cu_sq));
      (D27C.f[dM00])[kw   ] =   c2o27* (c3o1*(-Mx        )+Conc_C*(c1o1+c9o2*(-vx1        )*(-vx1        )-cu_sq));
      (D27C.f[d0P0])[kzero] =   c2o27* (c3o1*(     My    )+Conc_C*(c1o1+c9o2*(     vx2    )*(     vx2    )-cu_sq));
      (D27C.f[d0M0])[ks   ] =   c2o27* (c3o1*(    -My    )+Conc_C*(c1o1+c9o2*(    -vx2    )*(    -vx2    )-cu_sq));
      (D27C.f[d00P])[kzero] =   c2o27* (c3o1*(         Mz)+Conc_C*(c1o1+c9o2*(         vx3)*(         vx3)-cu_sq));
      (D27C.f[d00M])[kb   ] =   c2o27* (c3o1*(        -Mz)+Conc_C*(c1o1+c9o2*(        -vx3)*(        -vx3)-cu_sq));
      (D27C.f[dPP0])[kzero] =   c1o54* (c3o1*( Mx +My    )+Conc_C*(c1o1+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq));
      (D27C.f[dMM0])[ksw  ] =   c1o54* (c3o1*(-Mx -My    )+Conc_C*(c1o1+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq));
      (D27C.f[dPM0])[ks   ] =   c1o54* (c3o1*( Mx -My    )+Conc_C*(c1o1+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq));
      (D27C.f[dMP0])[kw   ] =   c1o54* (c3o1*(-Mx +My    )+Conc_C*(c1o1+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq));
      (D27C.f[dP0P])[kzero] =   c1o54* (c3o1*( Mx     +Mz)+Conc_C*(c1o1+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq));
      (D27C.f[dM0M])[kbw  ] =   c1o54* (c3o1*(-Mx     -Mz)+Conc_C*(c1o1+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq));
      (D27C.f[dP0M])[kb   ] =   c1o54* (c3o1*( Mx     -Mz)+Conc_C*(c1o1+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq));
      (D27C.f[dM0P])[kw   ] =   c1o54* (c3o1*(-Mx     +Mz)+Conc_C*(c1o1+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq));
      (D27C.f[d0PP])[kzero] =   c1o54* (c3o1*(     My +Mz)+Conc_C*(c1o1+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq));
      (D27C.f[d0MM])[kbs  ] =   c1o54* (c3o1*(    -My -Mz)+Conc_C*(c1o1+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq));
      (D27C.f[d0PM])[kb   ] =   c1o54* (c3o1*(     My -Mz)+Conc_C*(c1o1+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq));
      (D27C.f[d0MP])[ks   ] =   c1o54* (c3o1*(    -My +Mz)+Conc_C*(c1o1+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq));
      (D27C.f[dPPP])[kzero] =   c1o216*(c3o1*( Mx +My +Mz)+Conc_C*(c1o1+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq));
      (D27C.f[dMMM])[kbsw ] =   c1o216*(c3o1*(-Mx -My -Mz)+Conc_C*(c1o1+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq));
      (D27C.f[dPPM])[kb   ] =   c1o216*(c3o1*( Mx +My -Mz)+Conc_C*(c1o1+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq));
      (D27C.f[dMMP])[ksw  ] =   c1o216*(c3o1*(-Mx -My +Mz)+Conc_C*(c1o1+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq));
      (D27C.f[dPMP])[ks   ] =   c1o216*(c3o1*( Mx -My +Mz)+Conc_C*(c1o1+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq));
      (D27C.f[dMPM])[kbw  ] =   c1o216*(c3o1*(-Mx +My -Mz)+Conc_C*(c1o1+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq));
      (D27C.f[dPMM])[kbs  ] =   c1o216*(c3o1*( Mx -My -Mz)+Conc_C*(c1o1+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq));
      (D27C.f[dMPP])[kw   ] =   c1o216*(c3o1*(-Mx +My +Mz)+Conc_C*(c1o1+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq));

   }
}

