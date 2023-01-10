#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

#include <map>
#include "basics/Core/DataTypes.h"

namespace vf::lbm::dir
{

static constexpr size_t STARTDIR = 0;
static constexpr size_t ENDDIR = 26;

// used in the CPU and the GPU version
static constexpr size_t DIR_000 = 0;
static constexpr size_t DIR_P00 = 1;
static constexpr size_t DIR_M00 = 2;
static constexpr size_t DIR_0P0 = 3;
static constexpr size_t DIR_0M0 = 4;
static constexpr size_t DIR_00P = 5;
static constexpr size_t DIR_00M = 6;
static constexpr size_t DIR_PP0 = 7;
static constexpr size_t DIR_MM0 = 8;
static constexpr size_t DIR_PM0 = 9;
static constexpr size_t DIR_MP0 = 10;
static constexpr size_t DIR_P0P = 11;
static constexpr size_t DIR_M0M = 12;
static constexpr size_t DIR_P0M = 13;
static constexpr size_t DIR_M0P = 14;
static constexpr size_t DIR_0PP = 15;
static constexpr size_t DIR_0MM = 16;
static constexpr size_t DIR_0PM = 17;
static constexpr size_t DIR_0MP = 18;
static constexpr size_t DIR_PPP = 19;
static constexpr size_t DIR_MPP = 20;
static constexpr size_t DIR_PMP = 21;
static constexpr size_t DIR_MMP = 22;
static constexpr size_t DIR_PPM = 23;
static constexpr size_t DIR_MPM = 24;
static constexpr size_t DIR_PMM = 25;
static constexpr size_t DIR_MMM = 26;

static constexpr size_t INV_P00 = DIR_M00;
static constexpr size_t INV_M00 = DIR_P00;
static constexpr size_t INV_0P0 = DIR_0M0;
static constexpr size_t INV_0M0 = DIR_0P0;
static constexpr size_t INV_00P = DIR_00M;
static constexpr size_t INV_00M = DIR_00P;
static constexpr size_t INV_PP0 = DIR_MM0;
static constexpr size_t INV_MM0 = DIR_PP0;
static constexpr size_t INV_PM0 = DIR_MP0;
static constexpr size_t INV_MP0 = DIR_PM0;
static constexpr size_t INV_P0P = DIR_M0M;
static constexpr size_t INV_M0M = DIR_P0P;
static constexpr size_t INV_P0M = DIR_M0P;
static constexpr size_t INV_M0P = DIR_P0M;
static constexpr size_t INV_0PP = DIR_0MM;
static constexpr size_t INV_0MM = DIR_0PP;
static constexpr size_t INV_0PM = DIR_0MP;
static constexpr size_t INV_0MP = DIR_0PM;
static constexpr size_t INV_PPP = DIR_MMM;
static constexpr size_t INV_MPP = DIR_PMM;
static constexpr size_t INV_PMP = DIR_MPM;
static constexpr size_t INV_MMP = DIR_PPM;
static constexpr size_t INV_PPM = DIR_MMP;
static constexpr size_t INV_MPM = DIR_PMP;
static constexpr size_t INV_PMM = DIR_MPP;
static constexpr size_t INV_MMM = DIR_PPP;

static constexpr size_t SGD_P00 = 0;
static constexpr size_t SGD_M00 = 1;
static constexpr size_t SGD_0P0 = 2;
static constexpr size_t SGD_0M0 = 3;
static constexpr size_t SGD_00P = 4;
static constexpr size_t SGD_00M = 5;
static constexpr size_t SGD_PP0 = 6;
static constexpr size_t SGD_MM0 = 7;
static constexpr size_t SGD_PM0 = 8;
static constexpr size_t SGD_MP0 = 9;
static constexpr size_t SGD_P0P = 10;
static constexpr size_t SGD_M0M = 11;
static constexpr size_t SGD_P0M = 12;
static constexpr size_t SGD_M0P = 13;
static constexpr size_t SGD_0PP = 14;
static constexpr size_t SGD_0MM = 15;
static constexpr size_t SGD_0PM = 16;
static constexpr size_t SGD_0MP = 17;
static constexpr size_t SGD_PPP = 18;
static constexpr size_t SGD_MPP = 19;
static constexpr size_t SGD_PMP = 20;
static constexpr size_t SGD_MMP = 21;
static constexpr size_t SGD_PPM = 22;
static constexpr size_t SGD_MPM = 23;
static constexpr size_t SGD_PMM = 24;
static constexpr size_t SGD_MMM = 25;

struct countersForPointerChasing{
    uint counterInverse;
    uint counterX;
    uint counterY;
    uint counterZ;
};

const std::map<const size_t, const countersForPointerChasing> mapForPointerChasing = 
{
    {DIR_000, countersForPointerChasing{0, 0, 0, 0}},
    {DIR_P00, countersForPointerChasing{0, 1, 0, 0}},
    {DIR_M00, countersForPointerChasing{1, 0, 1, 1}},
    {DIR_0P0, countersForPointerChasing{0, 0, 1, 0}},
    {DIR_0M0, countersForPointerChasing{1, 1, 0, 1}},
    {DIR_00P, countersForPointerChasing{0, 0, 0, 1}},
    {DIR_00M, countersForPointerChasing{1, 1, 1, 0}},

    {DIR_PP0, countersForPointerChasing{0, 1, 1, 0}},
    {DIR_MM0, countersForPointerChasing{1, 0, 0, 1}},
    {DIR_PM0, countersForPointerChasing{1, 2, 0, 1}},
    {DIR_MP0, countersForPointerChasing{1, 0, 2, 1}},
    {DIR_P0P, countersForPointerChasing{0, 1, 0, 1}},
    {DIR_M0M, countersForPointerChasing{1, 0, 1, 0}},
    {DIR_P0M, countersForPointerChasing{1, 2, 1, 0}},
    {DIR_M0P, countersForPointerChasing{1, 0, 1, 2}},
    {DIR_0PP, countersForPointerChasing{0, 0, 1, 1}},
    {DIR_0MM, countersForPointerChasing{1, 1, 0, 0}},
    {DIR_0PM, countersForPointerChasing{1, 1, 2, 0}},
    {DIR_0MP, countersForPointerChasing{1, 1, 0, 2}},

    {DIR_PPP, countersForPointerChasing{0, 1, 1, 1}},
    {DIR_MPP, countersForPointerChasing{1, 0, 2, 2}},
    {DIR_PMP, countersForPointerChasing{1, 2, 0, 2}},
    {DIR_MMP, countersForPointerChasing{1, 0, 0, 2}},
    {DIR_PPM, countersForPointerChasing{1, 2, 2, 0}},
    {DIR_MPM, countersForPointerChasing{1, 0, 2, 0}},
    {DIR_PMM, countersForPointerChasing{1, 2, 0, 0}},
    {DIR_MMM, countersForPointerChasing{1, 0, 0, 0}}
};



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
static constexpr int ZZZ = DIR_000;
static constexpr int PZZ = DIR_P00;
static constexpr int MZZ = DIR_M00;
static constexpr int ZPZ = DIR_0P0;
static constexpr int ZMZ = DIR_0M0;
static constexpr int ZZP = DIR_00P;
static constexpr int ZZM = DIR_00M;
static constexpr int PPZ = DIR_PP0;
static constexpr int MMZ = DIR_MM0;
static constexpr int PMZ = DIR_PM0;
static constexpr int MPZ = DIR_MP0;
static constexpr int PZP = DIR_P0P;
static constexpr int MZM = DIR_M0M;
static constexpr int PZM = DIR_P0M;
static constexpr int MZP = DIR_M0P;
static constexpr int ZPP = DIR_0PP;
static constexpr int ZMM = DIR_0MM;
static constexpr int ZPM = DIR_0PM;
static constexpr int ZMP = DIR_0MP;
static constexpr int PPP = DIR_PPP;
static constexpr int MPP = DIR_MPP;
static constexpr int PMP = DIR_PMP;
static constexpr int MMP = DIR_MMP;
static constexpr int PPM = DIR_PPM;
static constexpr int MPM = DIR_MPM;
static constexpr int PMM = DIR_PMM;
static constexpr int MMM = DIR_MMM;
}
#endif
