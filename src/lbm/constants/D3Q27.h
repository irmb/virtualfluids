#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

namespace vf
{
namespace lbm
{
namespace dir
{

static constexpr int E    = 0;
static constexpr int W    = 1;
static constexpr int N    = 2;
static constexpr int S    = 3;
static constexpr int T    = 4;
static constexpr int B    = 5;
static constexpr int NE   = 6;
static constexpr int SW   = 7;
static constexpr int SE   = 8;
static constexpr int NW   = 9;
static constexpr int TE   = 10;
static constexpr int BW   = 11;
static constexpr int BE   = 12;
static constexpr int TW   = 13;
static constexpr int TN   = 14;
static constexpr int BS   = 15;
static constexpr int BN   = 16;
static constexpr int TS   = 17;
static constexpr int TNE  = 18;
static constexpr int TNW  = 19;
static constexpr int TSE  = 20;
static constexpr int TSW  = 21;
static constexpr int BNE  = 22;
static constexpr int BNW  = 23;
static constexpr int BSE  = 24;
static constexpr int BSW  = 25;
static constexpr int REST = 26;

static constexpr int PZZ = 0;
static constexpr int MZZ = 1;
static constexpr int ZPZ = 2;
static constexpr int ZMZ = 3;
static constexpr int ZZP = 4;
static constexpr int ZZM = 5;
static constexpr int PPZ = 6;
static constexpr int MMZ = 7;
static constexpr int PMZ = 8;
static constexpr int MPZ = 9;
static constexpr int PZP = 10;
static constexpr int MZM = 11;
static constexpr int PZM = 12;
static constexpr int MZP = 13;
static constexpr int ZPP = 14;
static constexpr int ZMM = 15;
static constexpr int ZPM = 16;
static constexpr int ZMP = 17;
static constexpr int PPP = 18;
static constexpr int MPP = 19;
static constexpr int PMP = 20;
static constexpr int MMP = 21;
static constexpr int PPM = 22;
static constexpr int MPM = 23;
static constexpr int PMM = 24;
static constexpr int MMM = 25;
static constexpr int ZZZ = 26;

static constexpr int DIR_000 = 0;
static constexpr int DIR_P00 = 1;
static constexpr int DIR_M00 = 2;
static constexpr int DIR_0P0 = 3;
static constexpr int DIR_0M0 = 4;
static constexpr int DIR_00P = 5;
static constexpr int DIR_00M = 6;
static constexpr int DIR_PP0 = 7;
static constexpr int DIR_MM0 = 8;
static constexpr int DIR_PM0 = 9;
static constexpr int DIR_MP0 = 10;
static constexpr int DIR_P0P = 11;
static constexpr int DIR_M0M = 12;
static constexpr int DIR_P0M = 13;
static constexpr int DIR_M0P = 14;
static constexpr int DIR_0PP = 15;
static constexpr int DIR_0MM = 16;
static constexpr int DIR_0PM = 17;
static constexpr int DIR_0MP = 18;
static constexpr int DIR_PPP = 19;
static constexpr int DIR_MPP = 20;
static constexpr int DIR_PMP = 21;
static constexpr int DIR_MMP = 22;
static constexpr int DIR_PPM = 23;
static constexpr int DIR_MPM = 24;
static constexpr int DIR_PMM = 25;
static constexpr int DIR_MMM = 26;

static constexpr int INV_P00 = DIR_M00;
static constexpr int INV_M00 = DIR_P00;
static constexpr int INV_0P0 = DIR_0M0;
static constexpr int INV_0M0 = DIR_0P0;
static constexpr int INV_00P = DIR_00M;
static constexpr int INV_00M = DIR_00P;
static constexpr int INV_PP0 = DIR_MM0;
static constexpr int INV_MM0 = DIR_PP0;
static constexpr int INV_PM0 = DIR_MP0;
static constexpr int INV_MP0 = DIR_PM0;
static constexpr int INV_P0P = DIR_M0M;
static constexpr int INV_M0M = DIR_P0P;
static constexpr int INV_P0M = DIR_M0P;
static constexpr int INV_M0P = DIR_P0M;
static constexpr int INV_0PP = DIR_0MM;
static constexpr int INV_0MM = DIR_0PP;
static constexpr int INV_0PM = DIR_0MP;
static constexpr int INV_0MP = DIR_0PM;
static constexpr int INV_PPP = DIR_MMM;
static constexpr int INV_MPP = DIR_PMM;
static constexpr int INV_PMP = DIR_MPM;
static constexpr int INV_MMP = DIR_PPM;
static constexpr int INV_PPM = DIR_MMP;
static constexpr int INV_MPM = DIR_PMP;
static constexpr int INV_PMM = DIR_MPP;
static constexpr int INV_MMM = DIR_PPP;

static constexpr int SGD_P00 = 0;
static constexpr int SGD_M00 = 1;
static constexpr int SGD_0P0 = 2;
static constexpr int SGD_0M0 = 3;
static constexpr int SGD_00P = 4;
static constexpr int SGD_00M = 5;
static constexpr int SGD_PP0 = 6;
static constexpr int SGD_MM0 = 7;
static constexpr int SGD_PM0 = 8;
static constexpr int SGD_MP0 = 9;
static constexpr int SGD_P0P = 10;
static constexpr int SGD_M0M = 11;
static constexpr int SGD_P0M = 12;
static constexpr int SGD_M0P = 13;
static constexpr int SGD_0PP = 14;
static constexpr int SGD_0MM = 15;
static constexpr int SGD_0PM = 16;
static constexpr int SGD_0MP = 17;
static constexpr int SGD_PPP = 18;
static constexpr int SGD_MPP = 19;
static constexpr int SGD_PMP = 20;
static constexpr int SGD_MMP = 21;
static constexpr int SGD_PPM = 22;
static constexpr int SGD_MPM = 23;
static constexpr int SGD_PMM = 24;
static constexpr int SGD_MMM = 25;

}
}
}

#endif
