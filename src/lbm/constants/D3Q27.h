#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

namespace vf::lbm::dir
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

}
#endif
