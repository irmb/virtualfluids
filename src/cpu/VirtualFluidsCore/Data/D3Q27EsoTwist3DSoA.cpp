#include "D3Q27EsoTwist3DSoA.h"
#include "EsoTwistD3Q27System.h"
#include <D3Q27System.h>

D3Q27EsoTwist3DSoA::D3Q27EsoTwist3DSoA() = default;
//////////////////////////////////////////////////////////////////////////
D3Q27EsoTwist3DSoA::D3Q27EsoTwist3DSoA(const size_t &nx1, const size_t &nx2, const size_t &nx3, real value)
{
    this->NX1 = nx1;
    this->NX2 = nx2;
    this->NX3 = nx3;

    d.E = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.W = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.N = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.S = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.T = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.B = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.NE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.SW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.SE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.NW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TN = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BS = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BN = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TS = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TNE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TNW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TSE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.TSW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BNE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BNW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BSE = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.BSW = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(
        new CbArray3D<real, IndexerX3X2X1>(nx1 + 1, nx2 + 1, nx3 + 1, value));
    d.REST =
        CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(nx1, nx2, nx3, value));
}
//////////////////////////////////////////////////////////////////////////
D3Q27EsoTwist3DSoA::~D3Q27EsoTwist3DSoA() = default;
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::swap()
{
    std::swap(d.E, d.W);
    std::swap(d.N, d.S);
    std::swap(d.T, d.B);
    std::swap(d.NE, d.SW);
    std::swap(d.NW, d.SE);
    std::swap(d.TE, d.BW);
    std::swap(d.TW, d.BE);
    std::swap(d.TN, d.BS);
    std::swap(d.TS, d.BN);
    std::swap(d.TNE, d.BSW);
    std::swap(d.TNW, d.BSE);
    std::swap(d.TSE, d.BNW);
    std::swap(d.TSW, d.BNE);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::getDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    size_t x1p = x1 + 1;
    size_t x2p = x2 + 1;
    size_t x3p = x3 + 1;

    f[DIR_P00]   = (*d.E)(x1, x2, x3);
    f[DIR_0P0]   = (*d.N)(x1, x2, x3);
    f[DIR_00P]   = (*d.T)(x1, x2, x3);
    f[DIR_PP0]  = (*d.NE)(x1, x2, x3);
    f[DIR_MP0]  = (*d.NW)(x1p, x2, x3);
    f[DIR_P0P]  = (*d.TE)(x1, x2, x3);
    f[DIR_M0P]  = (*d.TW)(x1p, x2, x3);
    f[DIR_0PP]  = (*d.TN)(x1, x2, x3);
    f[DIR_0MP]  = (*d.TS)(x1, x2p, x3);
    f[DIR_PPP] = (*d.TNE)(x1, x2, x3);
    f[DIR_MPP] = (*d.TNW)(x1p, x2, x3);
    f[DIR_PMP] = (*d.TSE)(x1, x2p, x3);
    f[DIR_MMP] = (*d.TSW)(x1p, x2p, x3);

    f[DIR_M00]   = (*d.W)(x1p, x2, x3);
    f[DIR_0M0]   = (*d.S)(x1, x2p, x3);
    f[DIR_00M]   = (*d.B)(x1, x2, x3p);
    f[DIR_MM0]  = (*d.SW)(x1p, x2p, x3);
    f[DIR_PM0]  = (*d.SE)(x1, x2p, x3);
    f[DIR_M0M]  = (*d.BW)(x1p, x2, x3p);
    f[DIR_P0M]  = (*d.BE)(x1, x2, x3p);
    f[DIR_0MM]  = (*d.BS)(x1, x2p, x3p);
    f[DIR_0PM]  = (*d.BN)(x1, x2, x3p);
    f[DIR_MMM] = (*d.BSW)(x1p, x2p, x3p);
    f[DIR_PMM] = (*d.BSE)(x1, x2p, x3p);
    f[DIR_MPM] = (*d.BNW)(x1p, x2, x3p);
    f[DIR_PPM] = (*d.BNE)(x1, x2, x3p);

    f[DIR_000] = (*d.REST)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    size_t x1p = x1 + 1;
    size_t x2p = x2 + 1;
    size_t x3p = x3 + 1;

    (*d.E)(x1, x2, x3)     = f[INV_P00];
    (*d.N)(x1, x2, x3)     = f[INV_0P0];
    (*d.T)(x1, x2, x3)     = f[INV_00P];
    (*d.NE)(x1, x2, x3)    = f[INV_PP0];
    (*d.NW)(x1p, x2, x3)   = f[INV_MP0];
    (*d.TE)(x1, x2, x3)    = f[INV_P0P];
    (*d.TW)(x1p, x2, x3)   = f[INV_M0P];
    (*d.TN)(x1, x2, x3)    = f[INV_0PP];
    (*d.TS)(x1, x2p, x3)   = f[INV_0MP];
    (*d.TNE)(x1, x2, x3)   = f[INV_PPP];
    (*d.TNW)(x1p, x2, x3)  = f[INV_MPP];
    (*d.TSE)(x1, x2p, x3)  = f[INV_PMP];
    (*d.TSW)(x1p, x2p, x3) = f[INV_MMP];

    (*d.W)(x1p, x2, x3)     = f[INV_M00];
    (*d.S)(x1, x2p, x3)     = f[INV_0M0];
    (*d.B)(x1, x2, x3p)     = f[INV_00M];
    (*d.SW)(x1p, x2p, x3)   = f[INV_MM0];
    (*d.SE)(x1, x2p, x3)    = f[INV_PM0];
    (*d.BW)(x1p, x2, x3p)   = f[INV_M0M];
    (*d.BE)(x1, x2, x3p)    = f[INV_P0M];
    (*d.BS)(x1, x2p, x3p)   = f[INV_0MM];
    (*d.BN)(x1, x2, x3p)    = f[INV_0PM];
    (*d.BSW)(x1p, x2p, x3p) = f[INV_MMM];
    (*d.BSE)(x1, x2p, x3p)  = f[INV_PMM];
    (*d.BNW)(x1p, x2, x3p)  = f[INV_MPM];
    (*d.BNE)(x1, x2, x3p)   = f[INV_PPM];

    (*d.REST)(x1, x2, x3) = f[DIR_000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::getDistributionInv(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    f[INV_P00]   = (*d.E)(x1, x2, x3);
    f[INV_0P0]   = (*d.N)(x1, x2, x3);
    f[INV_00P]   = (*d.T)(x1, x2, x3);
    f[INV_PP0]  = (*d.NE)(x1, x2, x3);
    f[INV_MP0]  = (*d.NW)(x1 + 1, x2, x3);
    f[INV_P0P]  = (*d.TE)(x1, x2, x3);
    f[INV_M0P]  = (*d.TW)(x1 + 1, x2, x3);
    f[INV_0PP]  = (*d.TN)(x1, x2, x3);
    f[INV_0MP]  = (*d.TS)(x1, x2 + 1, x3);
    f[INV_PPP] = (*d.TNE)(x1, x2, x3);
    f[INV_MPP] = (*d.TNW)(x1 + 1, x2, x3);
    f[INV_PMP] = (*d.TSE)(x1, x2 + 1, x3);
    f[INV_MMP] = (*d.TSW)(x1 + 1, x2 + 1, x3);

    f[INV_M00]   = (*d.W)(x1 + 1, x2, x3);
    f[INV_0M0]   = (*d.S)(x1, x2 + 1, x3);
    f[INV_00M]   = (*d.B)(x1, x2, x3 + 1);
    f[INV_MM0]  = (*d.SW)(x1 + 1, x2 + 1, x3);
    f[INV_PM0]  = (*d.SE)(x1, x2 + 1, x3);
    f[INV_M0M]  = (*d.BW)(x1 + 1, x2, x3 + 1);
    f[INV_P0M]  = (*d.BE)(x1, x2, x3 + 1);
    f[INV_0MM]  = (*d.BS)(x1, x2 + 1, x3 + 1);
    f[INV_0PM]  = (*d.BN)(x1, x2, x3 + 1);
    f[INV_MMM] = (*d.BSW)(x1 + 1, x2 + 1, x3 + 1);
    f[INV_PMM] = (*d.BSE)(x1, x2 + 1, x3 + 1);
    f[INV_MPM] = (*d.BNW)(x1 + 1, x2, x3 + 1);
    f[INV_PPM] = (*d.BNE)(x1, x2, x3 + 1);

    f[DIR_000] = (*d.REST)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistributionInv(const real *const f, size_t x1, size_t x2, size_t x3)
{
    //(*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::DIR_P00];
    //(*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::DIR_0P0];
    //(*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::DIR_00P];
    //(*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::DIR_PP0];
    //(*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::DIR_MP0];
    //(*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::DIR_P0P];
    //(*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::DIR_M0P];
    //(*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::DIR_0PP];
    //(*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::DIR_0MP];
    //(*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::DIR_PPP];
    //(*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::DIR_MPP];
    //(*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::DIR_PMP];
    //(*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::DIR_MMP];

    //(*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::DIR_M00 ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::DIR_0M0 ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::DIR_00M ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::DIR_MM0];
    //(*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::DIR_PM0];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::DIR_M0M];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_P0M];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::DIR_0MM];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_0PM];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::DIR_DIR_PPM];

    //(*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                     unsigned long int direction)
{
    // bool directionFlag = false;
    // if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::DIR_P00]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::DIR_M00]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::DIR_0M0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::DIR_0P0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::DIR_00M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::DIR_00P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::DIR_MM0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::DIR_PP0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::DIR_MP0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::DIR_PM0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::DIR_M0M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::DIR_P0P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_M0P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::DIR_P0M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::DIR_0MM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::DIR_0PP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_0MP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::DIR_0PM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::BSW]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::DIR_PPP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::BSE]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::DIR_MPP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::BNW]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::DIR_PMP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
    //   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::DIR_DIR_PPM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::DIR_MMP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
    //   (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST]; directionFlag=true;
    //#ifdef _DEBUG
    //   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction)
{
    // switch (direction)
    //{
    // case D3Q27System::DIR_P00 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
    //   break;
    // case D3Q27System::DIR_M00 :
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_0M0 :
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_0P0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
    //   break;
    // case D3Q27System::DIR_00M :
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_00P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
    //   break;
    // case D3Q27System::DIR_MM0 :
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_PP0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::DIR_MP0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::DIR_PM0 :
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_M0M :
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_P0P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_M0P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_P0M :
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_0MM :
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_0PP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_0MP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_0PM :
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::BSW :
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_PPP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::BSE :
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_MPP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::BNW :
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::DIR_PMP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
    //   break;
    // case D3Q27System::DIR_DIR_PPM :
    //   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f;
    //   break;
    // case D3Q27System::DIR_MMP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f;
    //   break;
    // case D3Q27System::REST :
    //   (*this->zeroDistributions)(x1,x2,x3) = f;
    //   break;
    // default:
    //   UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //}
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistributionInvForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                        unsigned long int direction)
{
    //   bool directionFlag = false;
    //   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
    //      (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::DIR_P00]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::DIR_M00]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::DIR_0M0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
    //      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::DIR_0P0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::DIR_00M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
    //      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::DIR_00P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::DIR_MM0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
    //      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::DIR_PP0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
    //      (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::DIR_MP0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::DIR_PM0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::DIR_M0M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
    //      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::DIR_P0P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
    //      (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::DIR_M0P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_P0M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::DIR_0MM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
    //      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::DIR_0PP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
    //      (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::DIR_0MP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::DIR_0PM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
    //      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::DIR_PPP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
    //      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::DIR_MPP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
    //      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::DIR_PMP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1)= f[D3Q27System::DIR_DIR_PPM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
    //      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::DIR_MMP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
    //      (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST]; directionFlag=true;
    //#ifdef _DEBUG
    //   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setDistributionInvForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                        unsigned long int direction)
{
    // switch (direction)
    //{
    // case D3Q27System::DIR_P00 :
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_M00 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
    //   break;
    // case D3Q27System::DIR_0M0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
    //   break;
    // case D3Q27System::DIR_0P0 :
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_00M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
    //   break;
    // case D3Q27System::DIR_00P :
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_MM0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::DIR_PP0 :
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_MP0 :
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_PM0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::DIR_M0M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_P0P :
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_M0P :
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_P0M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_0MM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
    //   break;
    // case D3Q27System::DIR_0PP :
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::DIR_0MP :
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::DIR_0PM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::BSW :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::DIR_PPP :
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::BSE :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::DIR_MPP :
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::BNW :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
    //   break;
    // case D3Q27System::DIR_PMP :
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::DIR_DIR_PPM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f;
    //   break;
    // case D3Q27System::DIR_MMP :
    //   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f;
    //   break;
    // case D3Q27System::REST :
    //   (*this->zeroDistributions)(x1,x2,x3) = f;
    //   break;
    // default:
    //   UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //}
}
//////////////////////////////////////////////////////////////////////////
real D3Q27EsoTwist3DSoA::getDistributionInvForDirection(size_t /*x1*/, size_t /*x2*/, size_t /*x3*/,
                                                           int /*direction*/)
{
    // switch (direction)
    //{
    // case D3Q27System::DIR_P00 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    );
    // case D3Q27System::DIR_M00 :
    //   return (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3);
    // case D3Q27System::DIR_0M0 :
    //   return (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3);
    // case D3Q27System::DIR_0P0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    );
    // case D3Q27System::DIR_00M :
    //   return (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3);
    // case D3Q27System::DIR_00P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  );
    // case D3Q27System::DIR_MM0 :
    //   return (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3);
    // case D3Q27System::DIR_PP0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   );
    // case D3Q27System::DIR_MP0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   );
    // case D3Q27System::DIR_PM0 :
    //   return (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3);
    // case D3Q27System::DIR_M0M :
    //   return (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3);
    // case D3Q27System::DIR_P0P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 );
    // case D3Q27System::DIR_M0P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 );
    // case D3Q27System::DIR_P0M :
    //   return (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3);
    // case D3Q27System::DIR_0MM :
    //   return (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3);
    // case D3Q27System::DIR_0PP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 );
    // case D3Q27System::DIR_0MP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 );
    // case D3Q27System::DIR_0PM :
    //   return (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3);
    // case D3Q27System::BSW :
    //   return (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3);
    // case D3Q27System::DIR_PPP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
    // case D3Q27System::BSE :
    //   return (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3);
    // case D3Q27System::DIR_MPP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1);
    // case D3Q27System::BNW :
    //   return (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3);
    // case D3Q27System::DIR_PMP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1);
    // case D3Q27System::DIR_DIR_PPM :
    //   return (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
    // case D3Q27System::DIR_MMP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1);
    // case D3Q27System::REST :
    //   return (*this->zeroDistributions)(x1,x2,x3);
    // default:
    //   UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //}
    return 0;
}
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSoA::getNX1() const { return NX1; }
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSoA::getNX2() const { return NX2; }
//////////////////////////////////////////////////////////////////////////
size_t D3Q27EsoTwist3DSoA::getNX3() const { return NX3; }
//////////////////////////////////////////////////////////////////////////
Distributions D3Q27EsoTwist3DSoA::getDistributions() { return d; }
//////////////////////////////////////////////////////////////////////////
