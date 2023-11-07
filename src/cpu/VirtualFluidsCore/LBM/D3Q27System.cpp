#include "D3Q27System.h"

#include "lbm/MacroscopicQuantities.h"

namespace D3Q27System
{
//using namespace UbMath;
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

// index            0   1   2   3   4   5  6   7   8   9   10  11  12  13  14  15  16  17   18  19  20  21  22  23  24  25
// f:               E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
//const int DX1[] = { 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
//const int DX2[] = { 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1 };
//const int DX3[] = { 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1 };

// index            0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18   19  20  21  22  23  24  25  26
// f:             REST, E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
const int DX1[] = { 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   1, -1,  1, -1,  1, -1,  1, -1 };
const int DX2[] = { 0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   1,  1, -1, -1,  1,  1, -1, -1 };
const int DX3[] = { 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   1,  1,  1,  1, -1, -1, -1, -1 };

const real WEIGTH[] = { c8o27,  
                          c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  
                          c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,  c1o54,
                          c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216 };

const int INVDIR[] = { d000, INV_P00, INV_M00, INV_0P0, INV_0M0, INV_00P, INV_00M, 
                                INV_PP0, INV_MM0, INV_PM0, INV_MP0, INV_P0P, INV_M0M, INV_P0M, INV_M0P, INV_0PP, INV_0MM, INV_0PM, INV_0MP,
                                INV_PPP, INV_MPP, INV_PMP, INV_MMP, INV_PPM, INV_MPM, INV_PMM, INV_MMM };

// index             0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18
// direction:        E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE
// BSW
//const int EX1[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
//const int EX2[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1 };
//const int EX3[] = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1 };

//////////////////////////////////////////////////////////////////////////



real getDensity(const real *const &f /*[27]*/)
{
    return vf::lbm::getDensity(f);
}

real getIncompVelocityX1(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX1(f);
}

real getIncompVelocityX2(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX2(f);
}

real getIncompVelocityX3(const real *const &f /*[27]*/)
{
    return vf::lbm::getIncompressibleVelocityX3(f);
}





} // namespace D3Q27System