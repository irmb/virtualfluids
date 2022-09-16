#include "LBM/LB.h" 
#include <lbm/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/MacroscopicQuantities.h>

#include "VirtualFluids_GPU/Kernel/Utilities/DistributionHelper.cuh"
#include "VirtualFluids_GPU/GPU/KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

__global__ void QPrecursorDeviceCompZeroPress( 	int* subgridDistanceIndices,
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
    real f_W    = (dist.f[DIR_P00   ])[ke   ];
    real f_E    = (dist.f[DIR_M00   ])[kw   ];
    real f_S    = (dist.f[DIR_0P0   ])[kn   ];
    real f_N    = (dist.f[DIR_0M0   ])[ks   ];
    real f_B    = (dist.f[DIR_00P   ])[kt   ];
    real f_T    = (dist.f[DIR_00M   ])[kb   ];
    real f_SW   = (dist.f[DIR_PP0  ])[kne  ];
    real f_NE   = (dist.f[DIR_MM0  ])[ksw  ];
    real f_NW   = (dist.f[DIR_PM0  ])[kse  ];
    real f_SE   = (dist.f[DIR_MP0  ])[knw  ];
    real f_BW   = (dist.f[DIR_P0P  ])[kte  ];
    real f_TE   = (dist.f[DIR_M0M  ])[kbw  ];
    real f_TW   = (dist.f[DIR_P0M  ])[kbe  ];
    real f_BE   = (dist.f[DIR_M0P  ])[ktw  ];
    real f_BS   = (dist.f[DIR_0PP  ])[ktn  ];
    real f_TN   = (dist.f[DIR_0MM  ])[kbs  ];
    real f_TS   = (dist.f[DIR_0PM  ])[kbn  ];
    real f_BN   = (dist.f[DIR_0MP  ])[kts  ];
    real f_BSW  = (dist.f[DIR_PPP ])[ktne ];
    real f_BNE  = (dist.f[DIR_MMP ])[ktsw ];
    real f_BNW  = (dist.f[DIR_PMP ])[ktse ];
    real f_BSE  = (dist.f[DIR_MPP ])[ktnw ];
    real f_TSW  = (dist.f[DIR_PPM ])[kbne ];
    real f_TNE  = (dist.f[DIR_MMM ])[kbsw ];
    real f_TNW  = (dist.f[DIR_PMM ])[kbse ];
    real f_TSE  = (dist.f[DIR_MPM ])[kbnw ];
    
    SubgridDistances27 subgridD;
    getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
    
    ////////////////////////////////////////////////////////////////////////////////
      real drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                     f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                     f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[DIR_000])[kzero]); 

      real vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                      ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                      (f_E - f_W)) / (c1o1 + drho); 
         

      real vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                       ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                       (f_N - f_S)) / (c1o1 + drho); 

      real vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                       (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                       (f_T - f_B)) / (c1o1 + drho); 

    
    // if(k==16383 || k==0) printf("k %d kQ %d drho = %f u %f v %f w %f\n",k, KQK, drho, vx1, vx2, vx3);
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);
    //////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    //! - Update distributions with subgrid distance (q) between zero and one
    real feq, q, velocityLB, velocityBC;
    q = (subgridD.q[DIR_P00])[k];
    if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
    {
        velocityLB = vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloX;
        (dist.f[DIR_M00])[kw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_E, f_W, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_M00])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloX;
        (dist.f[DIR_P00])[ke] = getInterpolatedDistributionForVeloWithPressureBC(q, f_W, f_E, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_0P0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloY;
        (dist.f[DIR_0M0])[DIR_0M0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_N, f_S, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_0M0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloY;
        (dist.f[DIR_0P0])[kn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_S, f_N, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_00P])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloZ;
        (dist.f[DIR_00M])[kb] = getInterpolatedDistributionForVeloWithPressureBC(q, f_T, f_B, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_00M])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloZ;
        (dist.f[DIR_00P])[kt] = getInterpolatedDistributionForVeloWithPressureBC(q, f_B, f_T, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_PP0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloY;
        (dist.f[DIR_MM0])[ksw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NE, f_SW, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_MM0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloY;
        (dist.f[DIR_PP0])[kne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SW, f_NE, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_PM0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloY;
        (dist.f[DIR_MP0])[knw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SE, f_NW, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_MP0])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloY;
        (dist.f[DIR_PM0])[kse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NW, f_SE, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_P0P])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloZ;
        (dist.f[DIR_M0M])[kbw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TE, f_BW, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_M0M])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloZ;
        (dist.f[DIR_P0P])[kte] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BW, f_TE, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_P0M])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloZ;
        (dist.f[DIR_M0P])[ktw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BE, f_TW, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_M0P])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloZ;
        (dist.f[DIR_P0M])[kbe] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TW, f_BE, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0PP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY + VeloZ;
        (dist.f[DIR_0MM])[kbs] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TN, f_BS, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0MM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY - VeloZ;
        (dist.f[DIR_0PP])[ktn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BS, f_TN, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0PM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY - VeloZ;
        (dist.f[DIR_0MP])[kts] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BN, f_TS, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0MP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY + VeloZ;
        (dist.f[DIR_0PM])[kbn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TS, f_BN, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_PPP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY + VeloZ;
        (dist.f[DIR_MMM])[kbsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNE, f_BSW, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MMM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY - VeloZ;
        (dist.f[DIR_PPP])[ktne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSW, f_TNE, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PPM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY - VeloZ;
        (dist.f[DIR_MMP])[ktsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNE, f_TSW, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MMP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY + VeloZ;
        (dist.f[DIR_PPM])[kbne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSW, f_BNE, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PMP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY + VeloZ;
        (dist.f[DIR_MPM])[kbnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSE, f_BNW, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MPM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY - VeloZ;
        (dist.f[DIR_PMP])[ktse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNW, f_TSE, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PMM])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY - VeloZ;
        (dist.f[DIR_MPP])[ktnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSE, f_TNW, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MPP])[k];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY + VeloZ;
        (dist.f[DIR_PMM])[kbse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNW, f_BSE, feq, omega, drho, velocityBC, c1o216);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


__global__ void PrecursorDeviceEQ27( 	int* subgridDistanceIndices,
                                                int numberOfBCnodes,
                                                real omega,
                                                real* distributions,
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // based on BGK Plus Comp
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real f_W    = (dist.f[DIR_P00])[ke   ];
    real f_E    = (dist.f[DIR_M00])[kw   ];
    real f_S    = (dist.f[DIR_0P0])[kn   ];
    real f_N    = (dist.f[DIR_0M0])[ks   ];
    real f_B    = (dist.f[DIR_00P])[kt   ];
    real f_T    = (dist.f[DIR_00M])[kb   ];
    real f_SW   = (dist.f[DIR_PP0])[kne  ];
    real f_NE   = (dist.f[DIR_MM0])[ksw  ];
    real f_NW   = (dist.f[DIR_PM0])[kse  ];
    real f_SE   = (dist.f[DIR_MP0])[knw  ];
    real f_BW   = (dist.f[DIR_P0P])[kte  ];
    real f_TE   = (dist.f[DIR_M0M])[kbw  ];
    real f_TW   = (dist.f[DIR_P0M])[kbe  ];
    real f_BE   = (dist.f[DIR_M0P])[ktw  ];
    real f_BS   = (dist.f[DIR_0PP])[ktn  ];
    real f_TN   = (dist.f[DIR_0MM])[kbs  ];
    real f_TS   = (dist.f[DIR_0PM])[kbn  ];
    real f_BN   = (dist.f[DIR_0MP])[kts  ];
    real f_ZERO = (dist.f[DIR_000])[kzero];
    real f_BSW  = (dist.f[DIR_PPP])[ktne ];
    real f_BNE  = (dist.f[DIR_MMP])[ktsw ];
    real f_BNW  = (dist.f[DIR_PMP])[ktse ];
    real f_BSE  = (dist.f[DIR_MPP])[ktnw ];
    real f_TSW  = (dist.f[DIR_PPM])[kbne ];
    real f_TNE  = (dist.f[DIR_MMM])[kbsw ];
    real f_TNW  = (dist.f[DIR_PMM])[kbse ];
    real f_TSE  = (dist.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set macroscopic quantities
      //!
      real drho = c0o1;

      real vx1  = VeloX;          

      real vx2  = VeloY; 

      real vx3  = VeloZ; 

      real cusq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      f_ZERO  = c8o27*  (drho-(drho+c1o1)*cusq);
      f_E     = c2o27*  (drho+(drho+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f_W     = c2o27*  (drho+(drho+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f_N     = c2o27*  (drho+(drho+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f_S     = c2o27*  (drho+(drho+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f_T     = c2o27*  (drho+(drho+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f_B     = c2o27*  (drho+(drho+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f_NE    = c1o54*  (drho+(drho+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f_SW    = c1o54*  (drho+(drho+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f_SE    =  c1o54* (drho+(drho+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f_NW    =  c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f_TE    =  c1o54* (drho+(drho+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f_BW    =  c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f_BE    =  c1o54* (drho+(drho+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f_TW    =  c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f_TN    =  c1o54* (drho+(drho+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f_BS    =  c1o54* (drho+(drho+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f_BN    =  c1o54* (drho+(drho+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f_TS    =  c1o54* (drho+(drho+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f_TNE   =  c1o216*(drho+(drho+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f_BSW   =  c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f_BNE   =  c1o216*(drho+(drho+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f_TSW   =  c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f_TSE   =  c1o216*(drho+(drho+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f_BNW   =  c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f_BSE   =  c1o216*(drho+(drho+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f_TNW   =  c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes
      //!
      (dist.f[DIR_P00   ])[ke  ] = f_W   ;
      (dist.f[DIR_M00   ])[kw  ] = f_E   ;
      (dist.f[DIR_0P0   ])[kn  ] = f_S   ;
      (dist.f[DIR_0M0   ])[ks  ] = f_N   ;
      (dist.f[DIR_00P   ])[kt  ] = f_B   ;
      (dist.f[DIR_00M   ])[kb  ] = f_T   ;
      (dist.f[DIR_PP0  ])[kne  ] = f_SW  ;
      (dist.f[DIR_MM0  ])[ksw  ] = f_NE  ;
      (dist.f[DIR_PM0  ])[kse  ] = f_NW  ;
      (dist.f[DIR_MP0  ])[knw  ] = f_SE  ;
      (dist.f[DIR_P0P  ])[kte  ] = f_BW  ;
      (dist.f[DIR_M0M  ])[kbw  ] = f_TE  ;
      (dist.f[DIR_P0M  ])[kbe  ] = f_TW  ;
      (dist.f[DIR_M0P  ])[ktw  ] = f_BE  ;
      (dist.f[DIR_0PP  ])[ktn  ] = f_BS  ;
      (dist.f[DIR_0MM  ])[kbs  ] = f_TN  ;
      (dist.f[DIR_0PM  ])[kbn  ] = f_TS  ;
      (dist.f[DIR_0MP  ])[kts  ] = f_BN  ;
      (dist.f[DIR_000])[kzero] = f_ZERO;
      (dist.f[DIR_PPP ])[ktne ] = f_BSW ;
      (dist.f[DIR_MMP ])[ktsw ] = f_BNE ;
      (dist.f[DIR_PMP ])[ktse ] = f_BNW ;
      (dist.f[DIR_MPP ])[ktnw ] = f_BSE ;
      (dist.f[DIR_PPM ])[kbne ] = f_TSW ;
      (dist.f[DIR_MMM ])[kbsw ] = f_TNE ;
      (dist.f[DIR_PMM ])[kbse ] = f_TNW ;
      (dist.f[DIR_MPM ])[kbnw ] = f_TSE ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////