#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>
#include "VirtualFluids_GPU/Kernel/Utilities/DistributionHelper.cuh"

using namespace vf::lbm::constant;

extern "C" __global__ void QPrecursorDeviceCompZeroPress( 	int* k_Q,
															int kQ,
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
															real tRatio,
                                                            real velocityRatio,
															unsigned long long size_Mat,
															bool evenOrOdd)
{
    const unsigned  x = threadIdx.x;  // Globaler x-Index 
    const unsigned  y = blockIdx.x;   // Globaler y-Index 
    const unsigned  z = blockIdx.y;   // Globaler z-Index 

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx*(ny*z + y) + x;

    if(k>kQ) return;

    DistributionReferences27 D = vf::gpu::getDistributionReferences27(DD, size_Mat, evenOrOdd);
    
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
    // if(k==100)
        // printf("last u %f v %f next u %f v %f\n", vxLastInterpd, vyLastInterpd, vxNextInterpd, vyNextInterpd);
    real VeloX = ((1.f-tRatio)*vxLastInterpd + tRatio*vxNextInterpd)/velocityRatio;
    real VeloY = ((1.f-tRatio)*vyLastInterpd + tRatio*vyNextInterpd)/velocityRatio; 
    real VeloZ = ((1.f-tRatio)*vzLastInterpd + tRatio*vzNextInterpd)/velocityRatio;
    // From here on just a copy of the velocity boundary condition
    ////////////////////////////////////////////////////////////////////////////////
    real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
    real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
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
    real vx1, vx2, vx3, drho, feq, q;
    drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

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
    D = vf::gpu::getDistributionReferences27(DD, size_Mat, !evenOrOdd);

    q = q_dirE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirW])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirW])[kw]=zero;
    }

    q = q_dirW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirE])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirE])[ke]=zero;
    }

    q = q_dirN[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirS])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirS])[ks]=zero;
    }

    q = q_dirS[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirN])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirN])[kn]=zero;
    }

    q = q_dirT[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirB])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirB])[kb]=one;
    }

    q = q_dirB[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirT])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
        //(D.f[dirT])[kt]=zero;
    }

    q = q_dirNE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirSW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirSW])[ksw]=zero;
    }

    q = q_dirSW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirNE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirNE])[kne]=zero;
    }

    q = q_dirSE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirNW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirNW])[knw]=zero;
    }

    q = q_dirNW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
        (D.f[dirSE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirSE])[kse]=zero;
    }

    q = q_dirTE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirBW])[kbw]=zero;
    }

    q = q_dirBW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirTE])[kte]=zero;
    }

    q = q_dirBE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirTW])[ktw]=zero;
    }

    q = q_dirTW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirBE])[kbe]=zero;
    }

    q = q_dirTN[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirBS])[kbs]=zero;
    }

    q = q_dirBS[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirTN])[ktn]=zero;
    }

    q = q_dirBN[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirTS])[kts]=zero;
    }

    q = q_dirTS[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
        //(D.f[dirBN])[kbn]=zero;
    }

    q = q_dirTNE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirBSW])[kbsw]=zero;
    }

    q = q_dirBSW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirTNE])[ktne]=zero;
    }

    q = q_dirBNE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirTSW])[ktsw]=zero;
    }

    q = q_dirTSW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirBNE])[kbne]=zero;
    }

    q = q_dirTSE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirBNW])[kbnw]=zero;
    }

    q = q_dirBNW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirTSE])[ktse]=zero;
    }

    q = q_dirBSE[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirTNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirTNW])[ktnw]=zero;
    }

    q = q_dirTNW[k];
    if (q>=c0o1 && q<=c1o1)
    {
        feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
        (D.f[dirBSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
        //(D.f[dirBSE])[kbse]=zero;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////