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
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "lbm/MacroscopicQuantities.h"
#include "Utilities/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void PressureNonEquilibriumCompressible_Device(
    real* rhoBC,
    real* distributions,
    int* bcNodeIndices,
    int* bcNeighborIndices,
    int numberOfBCnodes,
    real omega1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    size_t direction)
{
   ////////////////////////////////////////////////////////////////////////////////
   //! The pressure boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   ////////////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(nodeIndex < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local pressure
      //!
      real rhoBClocal = rhoBC[nodeIndex];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int KQK  = bcNodeIndices[nodeIndex];
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
      //! - Set neighbor indices (necessary for indirect addressing) for neighboring node
      //!
      unsigned int K1QK  = bcNeighborIndices[nodeIndex];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions for neighboring node
      //!
      real f1_W    = (dist.f[dP00])[k1e   ];
      real f1_E    = (dist.f[dM00])[k1w   ];
      real f1_S    = (dist.f[d0P0])[k1n   ];
      real f1_N    = (dist.f[d0M0])[k1s   ];
      real f1_B    = (dist.f[d00P])[k1t   ];
      real f1_T    = (dist.f[d00M])[k1b   ];
      real f1_SW   = (dist.f[dPP0])[k1ne  ];
      real f1_NE   = (dist.f[dMM0])[k1sw  ];
      real f1_NW   = (dist.f[dPM0])[k1se  ];
      real f1_SE   = (dist.f[dMP0])[k1nw  ];
      real f1_BW   = (dist.f[dP0P])[k1te  ];
      real f1_TE   = (dist.f[dM0M])[k1bw  ];
      real f1_TW   = (dist.f[dP0M])[k1be  ];
      real f1_BE   = (dist.f[dM0P])[k1tw  ];
      real f1_BS   = (dist.f[d0PP])[k1tn  ];
      real f1_TN   = (dist.f[d0MM])[k1bs  ];
      real f1_TS   = (dist.f[d0PM])[k1bn  ];
      real f1_BN   = (dist.f[d0MP])[k1ts  ];
      real f1_ZERO = (dist.f[d000])[k1zero];
      real f1_BSW  = (dist.f[dPPP])[k1tne ];
      real f1_BNE  = (dist.f[dMMP])[k1tsw ];
      real f1_BNW  = (dist.f[dPMP])[k1tse ];
      real f1_BSE  = (dist.f[dMPP])[k1tnw ];
      real f1_TSW  = (dist.f[dPPM])[k1bne ];
      real f1_TNE  = (dist.f[dMMM])[k1bsw ];
      real f1_TNW  = (dist.f[dPMM])[k1bse ];
      real f1_TSE  = (dist.f[dMPM])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities (for neighboring node)
      //!
      real drho1 = f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                   f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW +
                   f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((dist.f[d000])[kzero]);

      real vx1  = (((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                   (f1_E - f1_W)) / (c1o1 + drho1);

      real vx2  = ((-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                   (f1_N - f1_S)) / (c1o1 + drho1);

      real vx3  = (((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                   (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                   (f1_T - f1_B)) / (c1o1 + drho1);

      real cusq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      //! subtract the equilibrium (eq) to obtain the non-equilibrium (neq) (for neighboring node)
      //!
      f1_ZERO  -= c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      ////////////////////////////////////////////////////////////////////////////////
      //! redefine drho1 with rhoBClocal
      //!
      drho1 = rhoBClocal;

      ////////////////////////////////////////////////////////////////////////////////
      //! add the equilibrium (eq), which is calculated with rhoBClocal (for neighboring node)
      //!
      f1_ZERO  += c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes (only for the relevant directions)
      //!

      switch (direction)
      {
         case dM00:
            (dist.f[dP00])[ke   ] = f1_W   ;
            (dist.f[d0P0])[kn   ] = f1_S   ;
            (dist.f[d0M0])[ks   ] = f1_N   ;
            (dist.f[d00P])[kt   ] = f1_B   ;
            (dist.f[d00M])[kb   ] = f1_T   ;
            (dist.f[dPP0])[kne  ] = f1_SW  ;
            (dist.f[dPM0])[kse  ] = f1_NW  ;
            (dist.f[dP0P])[kte  ] = f1_BW  ;
            (dist.f[dP0M])[kbe  ] = f1_TW  ;
            (dist.f[d0PP])[ktn  ] = f1_BS  ;
            (dist.f[d0MM])[kbs  ] = f1_TN  ;
            (dist.f[d0PM])[kbn  ] = f1_TS  ;
            (dist.f[d0MP])[kts  ] = f1_BN  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dPPP])[ktne ] = f1_BSW ;
            (dist.f[dPMP])[ktse ] = f1_BNW ;
            (dist.f[dPPM])[kbne ] = f1_TSW ;
            (dist.f[dPMM])[kbse ] = f1_TNW ;
            break;
         case dP00:
            (dist.f[dM00])[kw   ] = f1_E   ;
            (dist.f[d0P0])[kn   ] = f1_S   ;
            (dist.f[d0M0])[ks   ] = f1_N   ;
            (dist.f[d00P])[kt   ] = f1_B   ;
            (dist.f[d00M])[kb   ] = f1_T   ;
            (dist.f[dMM0])[ksw  ] = f1_NE  ;
            (dist.f[dMP0])[knw  ] = f1_SE  ;
            (dist.f[dM0M])[kbw  ] = f1_TE  ;
            (dist.f[dM0P])[ktw  ] = f1_BE  ;
            (dist.f[d0PP])[ktn  ] = f1_BS  ;
            (dist.f[d0MM])[kbs  ] = f1_TN  ;
            (dist.f[d0PM])[kbn  ] = f1_TS  ;
            (dist.f[d0MP])[kts  ] = f1_BN  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dMMP])[ktsw ] = f1_BNE ;
            (dist.f[dMPP])[ktnw ] = f1_BSE ;
            (dist.f[dMMM])[kbsw ] = f1_TNE ;
            (dist.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         case d0M0:
            (dist.f[dP00])[ke   ] = f1_W   ;
            (dist.f[dM00])[kw   ] = f1_E   ;
            (dist.f[d0P0])[kn   ] = f1_S   ;
            (dist.f[d00P])[kt   ] = f1_B   ;
            (dist.f[d00M])[kb   ] = f1_T   ;
            (dist.f[dPP0])[kne  ] = f1_SW  ;
            (dist.f[dMP0])[knw  ] = f1_SE  ;
            (dist.f[dP0P])[kte  ] = f1_BW  ;
            (dist.f[dM0M])[kbw  ] = f1_TE  ;
            (dist.f[dP0M])[kbe  ] = f1_TW  ;
            (dist.f[dM0P])[ktw  ] = f1_BE  ;
            (dist.f[d0PP])[ktn  ] = f1_BS  ;
            (dist.f[d0PM])[kbn  ] = f1_TS  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dPPP])[ktne ] = f1_BSW ;
            (dist.f[dMPP])[ktnw ] = f1_BSE ;
            (dist.f[dPPM])[kbne ] = f1_TSW ;
            (dist.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         case d0P0:
            (dist.f[dP00])[ke   ] = f1_W   ;
            (dist.f[dM00])[kw   ] = f1_E   ;
            (dist.f[d0M0])[ks   ] = f1_N   ;
            (dist.f[d00P])[kt   ] = f1_B   ;
            (dist.f[d00M])[kb   ] = f1_T   ;
            (dist.f[dMM0])[ksw  ] = f1_NE  ;
            (dist.f[dPM0])[kse  ] = f1_NW  ;
            (dist.f[dP0P])[kte  ] = f1_BW  ;
            (dist.f[dM0M])[kbw  ] = f1_TE  ;
            (dist.f[dP0M])[kbe  ] = f1_TW  ;
            (dist.f[dM0P])[ktw  ] = f1_BE  ;
            (dist.f[d0MM])[kbs  ] = f1_TN  ;
            (dist.f[d0MP])[kts  ] = f1_BN  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dMMP])[ktsw ] = f1_BNE ;
            (dist.f[dPMP])[ktse ] = f1_BNW ;
            (dist.f[dMMM])[kbsw ] = f1_TNE ;
            (dist.f[dPMM])[kbse ] = f1_TNW ;
            break;
         case d00M:
            (dist.f[dP00])[ke   ] = f1_W   ;
            (dist.f[dM00])[kw   ] = f1_E   ;
            (dist.f[d0P0])[kn   ] = f1_S   ;
            (dist.f[d0M0])[ks   ] = f1_N   ;
            (dist.f[d00P])[kt   ] = f1_B   ;
            (dist.f[dPP0])[kne  ] = f1_SW  ;
            (dist.f[dMM0])[ksw  ] = f1_NE  ;
            (dist.f[dPM0])[kse  ] = f1_NW  ;
            (dist.f[dMP0])[knw  ] = f1_SE  ;
            (dist.f[dP0P])[kte  ] = f1_BW  ;
            (dist.f[dM0P])[ktw  ] = f1_BE  ;
            (dist.f[d0PP])[ktn  ] = f1_BS  ;
            (dist.f[d0MP])[kts  ] = f1_BN  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dPPP])[ktne ] = f1_BSW ;
            (dist.f[dMMP])[ktsw ] = f1_BNE ;
            (dist.f[dPMP])[ktse ] = f1_BNW ;
            (dist.f[dMPP])[ktnw ] = f1_BSE ;
            break;
         case d00P:
            (dist.f[dP00])[ke   ] = f1_W   ;
            (dist.f[dM00])[kw   ] = f1_E   ;
            (dist.f[d0P0])[kn   ] = f1_S   ;
            (dist.f[d0M0])[ks   ] = f1_N   ;
            (dist.f[d00M])[kb   ] = f1_T   ;
            (dist.f[dPP0])[kne  ] = f1_SW  ;
            (dist.f[dMM0])[ksw  ] = f1_NE  ;
            (dist.f[dPM0])[kse  ] = f1_NW  ;
            (dist.f[dMP0])[knw  ] = f1_SE  ;
            (dist.f[dM0M])[kbw  ] = f1_TE  ;
            (dist.f[dP0M])[kbe  ] = f1_TW  ;
            (dist.f[d0MM])[kbs  ] = f1_TN  ;
            (dist.f[d0PM])[kbn  ] = f1_TS  ;
            (dist.f[d000])[kzero] = f1_ZERO;
            (dist.f[dPPM])[kbne ] = f1_TSW ;
            (dist.f[dMMM])[kbsw ] = f1_TNE ;
            (dist.f[dPMM])[kbse ] = f1_TNW ;
            (dist.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         default:
            break; 
      }
   }
}

