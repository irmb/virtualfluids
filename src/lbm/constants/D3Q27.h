#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

#include <map>

#include <basics/DataTypes.h>

namespace vf::lbm::dir
{

static constexpr size_t STARTDIR = 0;
static constexpr size_t ENDDIR = 26;

// used in the CPU and the GPU version
static constexpr size_t d000 = 0;
static constexpr size_t dP00 = 1;
static constexpr size_t dM00 = 2;
static constexpr size_t d0P0 = 3;
static constexpr size_t d0M0 = 4;
static constexpr size_t d00P = 5;
static constexpr size_t d00M = 6;
static constexpr size_t dPP0 = 7;
static constexpr size_t dMM0 = 8;
static constexpr size_t dPM0 = 9;
static constexpr size_t dMP0 = 10;
static constexpr size_t dP0P = 11;
static constexpr size_t dM0M = 12;
static constexpr size_t dP0M = 13;
static constexpr size_t dM0P = 14;
static constexpr size_t d0PP = 15;
static constexpr size_t d0MM = 16;
static constexpr size_t d0PM = 17;
static constexpr size_t d0MP = 18;
static constexpr size_t dPPP = 19;
static constexpr size_t dMPP = 20;
static constexpr size_t dPMP = 21;
static constexpr size_t dMMP = 22;
static constexpr size_t dPPM = 23;
static constexpr size_t dMPM = 24;
static constexpr size_t dPMM = 25;
static constexpr size_t dMMM = 26;

static constexpr size_t iP00 = dM00;
static constexpr size_t iM00 = dP00;
static constexpr size_t i0P0 = d0M0;
static constexpr size_t i0M0 = d0P0;
static constexpr size_t i00P = d00M;
static constexpr size_t i00M = d00P;
static constexpr size_t iPP0 = dMM0;
static constexpr size_t iMM0 = dPP0;
static constexpr size_t iPM0 = dMP0;
static constexpr size_t iMP0 = dPM0;
static constexpr size_t iP0P = dM0M;
static constexpr size_t iM0M = dP0P;
static constexpr size_t iP0M = dM0P;
static constexpr size_t iM0P = dP0M;
static constexpr size_t i0PP = d0MM;
static constexpr size_t i0MM = d0PP;
static constexpr size_t i0PM = d0MP;
static constexpr size_t i0MP = d0PM;
static constexpr size_t iPPP = dMMM;
static constexpr size_t iMPP = dPMM;
static constexpr size_t iPMP = dMPM;
static constexpr size_t iMMP = dPPM;
static constexpr size_t iPPM = dMMP;
static constexpr size_t iMPM = dPMP;
static constexpr size_t iPMM = dMPP;
static constexpr size_t iMMM = dPPP;

struct countersForPointerChasing
{
    uint counterInverse;
    uint counterX;
    uint counterY;
    uint counterZ;
};

const std::map<const size_t, const countersForPointerChasing> mapForPointerChasing = 
{
    {d000, countersForPointerChasing{0, 0, 0, 0}},
    {dP00, countersForPointerChasing{0, 1, 0, 0}},
    {dM00, countersForPointerChasing{1, 0, 1, 1}},
    {d0P0, countersForPointerChasing{0, 0, 1, 0}},
    {d0M0, countersForPointerChasing{1, 1, 0, 1}},
    {d00P, countersForPointerChasing{0, 0, 0, 1}},
    {d00M, countersForPointerChasing{1, 1, 1, 0}},

    {dPP0, countersForPointerChasing{0, 1, 1, 0}},
    {dMM0, countersForPointerChasing{1, 0, 0, 1}},
    {dPM0, countersForPointerChasing{1, 2, 0, 1}},
    {dMP0, countersForPointerChasing{1, 0, 2, 1}},
    {dP0P, countersForPointerChasing{0, 1, 0, 1}},
    {dM0M, countersForPointerChasing{1, 0, 1, 0}},
    {dP0M, countersForPointerChasing{1, 2, 1, 0}},
    {dM0P, countersForPointerChasing{1, 0, 1, 2}},
    {d0PP, countersForPointerChasing{0, 0, 1, 1}},
    {d0MM, countersForPointerChasing{1, 1, 0, 0}},
    {d0PM, countersForPointerChasing{1, 1, 2, 0}},
    {d0MP, countersForPointerChasing{1, 1, 0, 2}},

    {dPPP, countersForPointerChasing{0, 1, 1, 1}},
    {dMPP, countersForPointerChasing{1, 0, 2, 2}},
    {dPMP, countersForPointerChasing{1, 2, 0, 2}},
    {dMMP, countersForPointerChasing{1, 0, 0, 2}},
    {dPPM, countersForPointerChasing{1, 2, 2, 0}},
    {dMPM, countersForPointerChasing{1, 0, 2, 0}},
    {dPMM, countersForPointerChasing{1, 2, 0, 0}},
    {dMMM, countersForPointerChasing{1, 0, 0, 0}}
};

}
#endif
