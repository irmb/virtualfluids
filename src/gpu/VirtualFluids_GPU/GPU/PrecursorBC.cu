#include "LBM/LB.h" 
#include <lbm/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/MacroscopicQuantities.h>

#include "VirtualFluids_GPU/Kernel/Utilities/DistributionHelper.cuh"
#include "VirtualFluids_GPU/GPU/KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

extern "C" __global__ void QPrecursorDeviceCompZeroPress( 	int* subgridDistanceIndices,
															int numberOfBCnodes,
                                                            int sizeQ,
                                                            real omega,
															real* distributions,
                                                            real* subgridDistances,
                                                            uint* neighborX, 
                                                            uint* neighborY, 
                                                            uint* neighborZ,
															uint* neighborsNT, 
															uint* neighborsNB,
															uint* neighborsST,
															uint* neighborsSB,
															real* weightsNT, 
															real* weightsNB,
															real* weightsST,
															real* weightsSB,
															real* vxLast, 
															real* vyLast, 
															real* vzLast,
															real* vxCurrent,
															real* vyCurrent,
															real* vzCurrent,
                                                            real velocityX,
                                                            real velocityY,
                                                            real velocityZ,
															real tRatio,
                                                            real velocityRatio,
															unsigned long long numberOfLBnodes,
															bool isEvenTimestep)
{
    const unsigned k = vf::gpu::getNodeIndex();

    if(k>=numberOfBCnodes) return;

    ////////////////////////////////////////////////////////////////////////////////
    // interpolation of velocity
    real vxLastInterpd, vyLastInterpd, vzLastInterpd; 
    real vxNextInterpd, vyNextInterpd, vzNextInterpd; 

    uint kNT = neighborsNT[k];
    real dNT = weightsNT[k];

    if(dNT < 1e6)
    {
        uint kNB = neighborsNB[k];
        uint kST = neighborsST[k];
        uint kSB = neighborsSB[k];

        real dNB = weightsNB[k];
        real dST = weightsST[k];
        real dSB = weightsSB[k];

        real invWeightSum = 1.f/(dNT+dNB+dST+dSB);

        vxLastInterpd = (vxLast[kNT]*dNT + vxLast[kNB]*dNB + vxLast[kST]*dST + vxLast[kSB]*dSB)*invWeightSum;
        vyLastInterpd = (vyLast[kNT]*dNT + vyLast[kNB]*dNB + vyLast[kST]*dST + vyLast[kSB]*dSB)*invWeightSum;
        vzLastInterpd = (vzLast[kNT]*dNT + vzLast[kNB]*dNB + vzLast[kST]*dST + vzLast[kSB]*dSB)*invWeightSum;

        vxNextInterpd = (vxCurrent[kNT]*dNT + vxCurrent[kNB]*dNB + vxCurrent[kST]*dST + vxCurrent[kSB]*dSB)*invWeightSum;
        vyNextInterpd = (vyCurrent[kNT]*dNT + vyCurrent[kNB]*dNB + vyCurrent[kST]*dST + vyCurrent[kSB]*dSB)*invWeightSum;
        vzNextInterpd = (vzCurrent[kNT]*dNT + vzCurrent[kNB]*dNB + vzCurrent[kST]*dST + vzCurrent[kSB]*dSB)*invWeightSum;
    }
    else
    {
        vxLastInterpd = vxLast[kNT];
        vyLastInterpd = vyLast[kNT];
        vzLastInterpd = vzLast[kNT];

        vxNextInterpd = vxCurrent[kNT];
        vyNextInterpd = vyCurrent[kNT];
        vzNextInterpd = vzCurrent[kNT];
    }

    // if(k==16300) printf("%f %f %f\n", vxLastInterpd, vyLastInterpd, vzLastInterpd);
    real VeloX = (velocityX + (1.f-tRatio)*vxLastInterpd + tRatio*vxNextInterpd)/velocityRatio;
    real VeloY = (velocityY + (1.f-tRatio)*vyLastInterpd + tRatio*vyNextInterpd)/velocityRatio; 
    real VeloZ = (velocityZ + (1.f-tRatio)*vzLastInterpd + tRatio*vzNextInterpd)/velocityRatio;
    // From here on just a copy of QVelDeviceCompZeroPress
    ////////////////////////////////////////////////////////////////////////////////

    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

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
    
    SubgridDistances27 subgridD;
    getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
    
    ////////////////////////////////////////////////////////////////////////////////
      real drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                     f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                     f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[REST])[kzero]); 

      real vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                      ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                      (f_E - f_W)) / (c1o1 + drho); 
         

      real vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                       ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                       (f_N - f_S)) / (c1o1 + drho); 

      real vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                       (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                       (f_T - f_B)) / (c1o1 + drho); 

    
    if(k==16383 || k==0) printf("k %d kQ %d drho = %f u %f v %f w %f\n",k, KQK, drho, vx1, vx2, vx3);
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);
    //////////////////////////////////////////////////////////////////////////

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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////