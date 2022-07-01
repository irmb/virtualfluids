#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>
#include <lbm/MacroscopicQuantities.h>

#include "VirtualFluids_GPU/Kernel/Utilities/DistributionHelper.cuh"

using namespace vf::lbm::constant;

extern "C" __global__ void QPrecursorDeviceCompZeroPress( 	int* k_Q,
															int numberOfBCNodes,
                                                            int sizeQ,
                                                            real om1,
															real* DD,
                                                            real* QQ,
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
															unsigned long long size_Mat,
															bool isEvenTimestep)
{
    const unsigned k = vf::gpu::getNodeIndex();

    if(k>=numberOfBCNodes) return;

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
    DistributionReferences27 Q = vf::gpu::getDistributionReferences27(QQ, sizeQ, true);
    
    vf::gpu::DistributionWrapper distWrapper = vf::gpu::DistributionWrapper(DD, size_Mat, !isEvenTimestep, k_Q[k], neighborX, neighborY, neighborZ);
    real (&f)[27] = distWrapper.distribution.f;

    ////////////////////////////////////////////////////////////////////////////////
    real drho = vf::lbm::getDensity(f);
    real vx1 = vf::lbm::getCompressibleVelocityX1(f, drho);
    real vx2 = vf::lbm::getCompressibleVelocityX2(f, drho);
    real vx3 = vf::lbm::getCompressibleVelocityX3(f, drho);

    if(k==16383 || k==0) printf("k %d kQ %d drho = %f u %f v %f w %f\n",k, k_Q[k], drho, vx1, vx2, vx3);
    real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

    // have to make a copy of the distributions so we can write to them later
    real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
        f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

    f_W    = f[vf::lbm::dir::W   ];
    f_E    = f[vf::lbm::dir::E   ];
    f_S    = f[vf::lbm::dir::S   ];
    f_N    = f[vf::lbm::dir::N   ];
    f_B    = f[vf::lbm::dir::B   ];
    f_T    = f[vf::lbm::dir::T   ];
    f_SW   = f[vf::lbm::dir::SW  ];
    f_NE   = f[vf::lbm::dir::NE  ];
    f_NW   = f[vf::lbm::dir::NW  ];
    f_SE   = f[vf::lbm::dir::SE  ];
    f_BW   = f[vf::lbm::dir::BW  ];
    f_TE   = f[vf::lbm::dir::TE  ];
    f_TW   = f[vf::lbm::dir::TW  ];
    f_BE   = f[vf::lbm::dir::BE  ];
    f_BS   = f[vf::lbm::dir::BS  ];
    f_TN   = f[vf::lbm::dir::TN  ];
    f_TS   = f[vf::lbm::dir::TS  ];
    f_BN   = f[vf::lbm::dir::BN  ];
    f_BSW  = f[vf::lbm::dir::BSW ];
    f_BNE  = f[vf::lbm::dir::BNE ];
    f_BNW  = f[vf::lbm::dir::BNW ];
    f_BSE  = f[vf::lbm::dir::BSE ];
    f_TSW  = f[vf::lbm::dir::TSW ];
    f_TNE  = f[vf::lbm::dir::TNE ];
    f_TNW  = f[vf::lbm::dir::TNW ];
    f_TSE  = f[vf::lbm::dir::TSE ];
    ////////////////////////////////////////////////////////////////////////////////
    real feq, q;

    //////////////////////////////////////////////////////////////////////////
    q = Q.f[dirE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::W]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::E]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirN][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::S]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirS][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::N]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirT][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::B]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirB][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::T]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
    }

    q = Q.f[dirNE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::SW]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirSW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::NE]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirSE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::NW]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirNW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::SE]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirTE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BW]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirBW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TE]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirBE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TW]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirTW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BE]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirTN][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BS]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirBS][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TN]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirBN][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TS]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirTS][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BN]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
    }

    q = Q.f[dirTNE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BSW]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirBSW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TNE]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirBNE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TSW]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirTSW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BNE]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirTSE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BNW]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirBNW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TSE]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirBSE][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::TNW]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
    }

    q = Q.f[dirTNW][k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
        f[vf::lbm::dir::BSE]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
    }
    //////////////////////////////////////////////////////////////////////////
    distWrapper.write();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////