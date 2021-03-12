#include "D3Q27System.h"

#include "lbm/MacroscopicQuantities.h"

namespace D3Q27System
{
using namespace UbMath;

// index            0   1   2   3   4   5  6   7   8   9   10  11  12  13  14  15  16  17  18  20  21  22  23  24  25  26
// f:               E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
const int DX1[] = { 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
const int DX2[] = { 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1 };
const int DX3[] = { 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1 };

const double WEIGTH[] = { c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  c1o54,  c1o54,  c1o54,
                          c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,
                          c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c8o27 };

const int INVDIR[] = { INV_E,   INV_W,   INV_N,   INV_S,   INV_T,   INV_B,   INV_NE,  INV_SW, INV_SE,
                       INV_NW,  INV_TE,  INV_BW,  INV_BE,  INV_TW,  INV_TN,  INV_BS,  INV_BN, INV_TS,
                       INV_TNE, INV_TNW, INV_TSE, INV_TSW, INV_BNE, INV_BNW, INV_BSE, INV_BSW };




LBMReal getDensity(const LBMReal *const &f /*[27]*/)
{
    return VF::LBM::getDensity(f);
}

LBMReal getIncompVelocityX1(const LBMReal *const &f /*[27]*/)
{
    return VF::LBM::getIncompressibleVelocityX1(f);
}

LBMReal getIncompVelocityX2(const LBMReal *const &f /*[27]*/)
{
    return VF::LBM::getIncompressibleVelocityX2(f);
}

LBMReal getIncompVelocityX3(const LBMReal *const &f /*[27]*/)
{
    return VF::LBM::getIncompressibleVelocityX3(f);
}





} // namespace D3Q27System