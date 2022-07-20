#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

namespace vf::lbm::dir
{

static constexpr int STARTDIR = 0;
static constexpr int ENDDIR   = 26;

// used in the CPU and the GPU version
static constexpr int REST = 0;
static constexpr int E    = 1;
static constexpr int W    = 2;
static constexpr int N    = 3;
static constexpr int S    = 4;
static constexpr int T    = 5;
static constexpr int B    = 6;
static constexpr int NE   = 7;
static constexpr int SW   = 8;
static constexpr int SE   = 9;
static constexpr int NW   = 10;
static constexpr int TE   = 11;
static constexpr int BW   = 12;
static constexpr int BE   = 13;
static constexpr int TW   = 14;
static constexpr int TN   = 15;
static constexpr int BS   = 16;
static constexpr int BN   = 17;
static constexpr int TS   = 18;
static constexpr int TNE  = 19;
static constexpr int TNW  = 20;
static constexpr int TSE  = 21;
static constexpr int TSW  = 22;
static constexpr int BNE  = 23;
static constexpr int BNW  = 24;
static constexpr int BSE  = 25;
static constexpr int BSW  = 26;

// used in the CPU version
// static constexpr int INV_P00 = DIR_M00;
// static constexpr int INV_M00 = DIR_P00;
// static constexpr int INV_0P0 = DIR_0M0;
// static constexpr int INV_0M0 = DIR_0P0;
// static constexpr int INV_00P = DIR_00M;
// static constexpr int INV_00M = DIR_00P;
// static constexpr int INV_PP0 = DIR_MM0;
// static constexpr int INV_MM0 = DIR_PP0;
// static constexpr int INV_PM0 = DIR_MP0;
// static constexpr int INV_MP0 = DIR_PM0;
// static constexpr int INV_P0P = DIR_M0M;
// static constexpr int INV_M0M = DIR_P0P;
// static constexpr int INV_P0M = DIR_M0P;
// static constexpr int INV_M0P = DIR_P0M;
// static constexpr int INV_0PP = DIR_0MM;
// static constexpr int INV_0MM = DIR_0PP;
// static constexpr int INV_0PM = DIR_0MP;
// static constexpr int INV_0MP = DIR_0PM;
// static constexpr int INV_PPP = DIR_MMM;
// static constexpr int INV_MPP = DIR_PMM;
// static constexpr int INV_PMP = DIR_MPM;
// static constexpr int INV_MMP = DIR_PPM;
// static constexpr int INV_PPM = DIR_MMP;
// static constexpr int INV_MPM = DIR_PMP;
// static constexpr int INV_PMM = DIR_MPP;
// static constexpr int INV_MMM = DIR_PPP;

// static constexpr int SGD_P00 = 0;
// static constexpr int SGD_M00 = 1;
// static constexpr int SGD_0P0 = 2;
// static constexpr int SGD_0M0 = 3;
// static constexpr int SGD_00P = 4;
// static constexpr int SGD_00M = 5;
// static constexpr int SGD_PP0 = 6;
// static constexpr int SGD_MM0 = 7;
// static constexpr int SGD_PM0 = 8;
// static constexpr int SGD_MP0 = 9;
// static constexpr int SGD_P0P = 10;
// static constexpr int SGD_M0M = 11;
// static constexpr int SGD_P0M = 12;
// static constexpr int SGD_M0P = 13;
// static constexpr int SGD_0PP = 14;
// static constexpr int SGD_0MM = 15;
// static constexpr int SGD_0PM = 16;
// static constexpr int SGD_0MP = 17;
// static constexpr int SGD_PPP = 18;
// static constexpr int SGD_MPP = 19;
// static constexpr int SGD_PMP = 20;
// static constexpr int SGD_MMP = 21;
// static constexpr int SGD_PPM = 22;
// static constexpr int SGD_MPM = 23;
// static constexpr int SGD_PMM = 24;
// static constexpr int SGD_MMM = 25;


// DEPRECATED
static constexpr int ZZZ = REST;
static constexpr int PZZ = E;
static constexpr int MZZ = W;
static constexpr int ZPZ = N;
static constexpr int ZMZ = S;
static constexpr int ZZP = T;
static constexpr int ZZM = B;
static constexpr int PPZ = NE;
static constexpr int MMZ = SW;
static constexpr int PMZ = SE;
static constexpr int MPZ = NW;
static constexpr int PZP = TE;
static constexpr int MZM = BW;
static constexpr int PZM = BE;
static constexpr int MZP = TW;
static constexpr int ZPP = TN;
static constexpr int ZMM = BS;
static constexpr int ZPM = BN;
static constexpr int ZMP = TS;
static constexpr int PPP = TNE;
static constexpr int MPP = TNW;
static constexpr int PMP = TSE;
static constexpr int MMP = TSW;
static constexpr int PPM = BNE;
static constexpr int MPM = BNW;
static constexpr int PMM = BSE;
static constexpr int MMM = BSW;
}
#endif
