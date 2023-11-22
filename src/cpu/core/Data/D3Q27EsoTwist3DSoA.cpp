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
void D3Q27EsoTwist3DSoA::getPreCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    size_t x1p = x1 + 1;
    size_t x2p = x2 + 1;
    size_t x3p = x3 + 1;

    f[dP00]   = (*d.E)(x1, x2, x3);
    f[d0P0]   = (*d.N)(x1, x2, x3);
    f[d00P]   = (*d.T)(x1, x2, x3);
    f[dPP0]  = (*d.NE)(x1, x2, x3);
    f[dMP0]  = (*d.NW)(x1p, x2, x3);
    f[dP0P]  = (*d.TE)(x1, x2, x3);
    f[dM0P]  = (*d.TW)(x1p, x2, x3);
    f[d0PP]  = (*d.TN)(x1, x2, x3);
    f[d0MP]  = (*d.TS)(x1, x2p, x3);
    f[dPPP] = (*d.TNE)(x1, x2, x3);
    f[dMPP] = (*d.TNW)(x1p, x2, x3);
    f[dPMP] = (*d.TSE)(x1, x2p, x3);
    f[dMMP] = (*d.TSW)(x1p, x2p, x3);

    f[dM00]   = (*d.W)(x1p, x2, x3);
    f[d0M0]   = (*d.S)(x1, x2p, x3);
    f[d00M]   = (*d.B)(x1, x2, x3p);
    f[dMM0]  = (*d.SW)(x1p, x2p, x3);
    f[dPM0]  = (*d.SE)(x1, x2p, x3);
    f[dM0M]  = (*d.BW)(x1p, x2, x3p);
    f[dP0M]  = (*d.BE)(x1, x2, x3p);
    f[d0MM]  = (*d.BS)(x1, x2p, x3p);
    f[d0PM]  = (*d.BN)(x1, x2, x3p);
    f[dMMM] = (*d.BSW)(x1p, x2p, x3p);
    f[dPMM] = (*d.BSE)(x1, x2p, x3p);
    f[dMPM] = (*d.BNW)(x1p, x2, x3p);
    f[dPPM] = (*d.BNE)(x1, x2, x3p);

    f[d000] = (*d.REST)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setPostCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    size_t x1p = x1 + 1;
    size_t x2p = x2 + 1;
    size_t x3p = x3 + 1;

    (*d.E)(x1, x2, x3)     = f[iP00];
    (*d.N)(x1, x2, x3)     = f[i0P0];
    (*d.T)(x1, x2, x3)     = f[i00P];
    (*d.NE)(x1, x2, x3)    = f[iPP0];
    (*d.NW)(x1p, x2, x3)   = f[iMP0];
    (*d.TE)(x1, x2, x3)    = f[iP0P];
    (*d.TW)(x1p, x2, x3)   = f[iM0P];
    (*d.TN)(x1, x2, x3)    = f[i0PP];
    (*d.TS)(x1, x2p, x3)   = f[i0MP];
    (*d.TNE)(x1, x2, x3)   = f[iPPP];
    (*d.TNW)(x1p, x2, x3)  = f[iMPP];
    (*d.TSE)(x1, x2p, x3)  = f[iPMP];
    (*d.TSW)(x1p, x2p, x3) = f[iMMP];

    (*d.W)(x1p, x2, x3)     = f[iM00];
    (*d.S)(x1, x2p, x3)     = f[i0M0];
    (*d.B)(x1, x2, x3p)     = f[i00M];
    (*d.SW)(x1p, x2p, x3)   = f[iMM0];
    (*d.SE)(x1, x2p, x3)    = f[iPM0];
    (*d.BW)(x1p, x2, x3p)   = f[iM0M];
    (*d.BE)(x1, x2, x3p)    = f[iP0M];
    (*d.BS)(x1, x2p, x3p)   = f[i0MM];
    (*d.BN)(x1, x2, x3p)    = f[i0PM];
    (*d.BSW)(x1p, x2p, x3p) = f[iMMM];
    (*d.BSE)(x1, x2p, x3p)  = f[iPMM];
    (*d.BNW)(x1p, x2, x3p)  = f[iMPM];
    (*d.BNE)(x1, x2, x3p)   = f[iPPM];

    (*d.REST)(x1, x2, x3) = f[d000];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::getPostCollisionDistribution(real *const f, size_t x1, size_t x2, size_t x3)
{
    using namespace vf::lbm::dir;

    f[iP00]   = (*d.E)(x1, x2, x3);
    f[i0P0]   = (*d.N)(x1, x2, x3);
    f[i00P]   = (*d.T)(x1, x2, x3);
    f[iPP0]  = (*d.NE)(x1, x2, x3);
    f[iMP0]  = (*d.NW)(x1 + 1, x2, x3);
    f[iP0P]  = (*d.TE)(x1, x2, x3);
    f[iM0P]  = (*d.TW)(x1 + 1, x2, x3);
    f[i0PP]  = (*d.TN)(x1, x2, x3);
    f[i0MP]  = (*d.TS)(x1, x2 + 1, x3);
    f[iPPP] = (*d.TNE)(x1, x2, x3);
    f[iMPP] = (*d.TNW)(x1 + 1, x2, x3);
    f[iPMP] = (*d.TSE)(x1, x2 + 1, x3);
    f[iMMP] = (*d.TSW)(x1 + 1, x2 + 1, x3);

    f[iM00]   = (*d.W)(x1 + 1, x2, x3);
    f[i0M0]   = (*d.S)(x1, x2 + 1, x3);
    f[i00M]   = (*d.B)(x1, x2, x3 + 1);
    f[iMM0]  = (*d.SW)(x1 + 1, x2 + 1, x3);
    f[iPM0]  = (*d.SE)(x1, x2 + 1, x3);
    f[iM0M]  = (*d.BW)(x1 + 1, x2, x3 + 1);
    f[iP0M]  = (*d.BE)(x1, x2, x3 + 1);
    f[i0MM]  = (*d.BS)(x1, x2 + 1, x3 + 1);
    f[i0PM]  = (*d.BN)(x1, x2, x3 + 1);
    f[iMMM] = (*d.BSW)(x1 + 1, x2 + 1, x3 + 1);
    f[iPMM] = (*d.BSE)(x1, x2 + 1, x3 + 1);
    f[iMPM] = (*d.BNW)(x1 + 1, x2, x3 + 1);
    f[iPPM] = (*d.BNE)(x1, x2, x3 + 1);

    f[d000] = (*d.REST)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setPreCollisionDistribution(const real *const f, size_t x1, size_t x2, size_t x3)
{
    //(*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::dP00];
    //(*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::d0P0];
    //(*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::d00P];
    //(*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::dPP0];
    //(*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::dMP0];
    //(*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::dP0P];
    //(*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::dM0P];
    //(*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::d0PP];
    //(*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::d0MP];
    //(*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::dPPP];
    //(*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::dMPP];
    //(*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::dPMP];
    //(*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::dMMP];

    //(*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::dM00 ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::d0M0 ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::d00M ];
    //(*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::dMM0];
    //(*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::dPM0];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::dM0M];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::dP0M];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::d0MM];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::d0PM];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
    //(*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::DIR_DIR_PPM];

    //(*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setPostCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                     unsigned long int direction)
{
    // bool directionFlag = false;
    // if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::dP00]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::dM00]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::d0M0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::d0P0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::d00M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::d00P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::dMM0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::dPP0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::dMP0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::dPM0]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::dM0M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::dP0P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::dM0P]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::dP0M]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::d0MM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::d0PP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::d0MP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::d0PM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::BSW]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::dPPP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::BSE]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::dMPP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::BNW]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::dPMP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
    //   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::DIR_DIR_PPM]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::dMMP]; directionFlag=true;
    // if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
    //   (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST]; directionFlag=true;
    //#ifdef _DEBUG
    //   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setPostCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3, int direction)
{
    // switch (direction)
    //{
    // case D3Q27System::dP00 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
    //   break;
    // case D3Q27System::dM00 :
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d0M0 :
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d0P0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
    //   break;
    // case D3Q27System::d00M :
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d00P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
    //   break;
    // case D3Q27System::dMM0 :
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dPP0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::dMP0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::dPM0 :
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::dM0M :
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dP0P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::dM0P :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::dP0M :
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::d0MM :
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d0PP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
    //   break;
    // case D3Q27System::d0MP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::d0PM :
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::BSW :
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dPPP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::BSE :
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::dMPP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::BNW :
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::dPMP :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
    //   break;
    // case D3Q27System::DIR_DIR_PPM :
    //   (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f;
    //   break;
    // case D3Q27System::dMMP :
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
void D3Q27EsoTwist3DSoA::setPreCollisionDistributionForDirection(const real *const f, size_t x1, size_t x2, size_t x3,
                                                        unsigned long int direction)
{
    //   bool directionFlag = false;
    //   if ((direction & EsoTwistD3Q27System::etE) == EsoTwistD3Q27System::etE)
    //      (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::dP00]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etW) == EsoTwistD3Q27System::etW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::dM00]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etS) == EsoTwistD3Q27System::etS)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::d0M0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etN) == EsoTwistD3Q27System::etN)
    //      (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::d0P0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etB) == EsoTwistD3Q27System::etB)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::d00M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etT) == EsoTwistD3Q27System::etT)
    //      (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::d00P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etSW) == EsoTwistD3Q27System::etSW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::dMM0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etNE) == EsoTwistD3Q27System::etNE)
    //      (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::dPP0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etNW) == EsoTwistD3Q27System::etNW)
    //      (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::dMP0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etSE) == EsoTwistD3Q27System::etSE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::dPM0]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBW) == EsoTwistD3Q27System::etBW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::dM0M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTE) == EsoTwistD3Q27System::etTE)
    //      (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::dP0P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTW) == EsoTwistD3Q27System::etTW)
    //      (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::dM0P]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBE) == EsoTwistD3Q27System::etBE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::dP0M]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBS) == EsoTwistD3Q27System::etBS)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::d0MM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTN) == EsoTwistD3Q27System::etTN)
    //      (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::d0PP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTS) == EsoTwistD3Q27System::etTS)
    //      (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::d0MP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBN) == EsoTwistD3Q27System::etBN)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::d0PM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBSW) == EsoTwistD3Q27System::etBSW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::BSW];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTNE) == EsoTwistD3Q27System::etTNE)
    //      (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::dPPP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBSE) == EsoTwistD3Q27System::etBSE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::BSE];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTNW) == EsoTwistD3Q27System::etTNW)
    //      (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::dMPP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBNW) == EsoTwistD3Q27System::etBNW)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::BNW];
    //      directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTSE) == EsoTwistD3Q27System::etTSE)
    //      (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::dPMP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etBNE) == EsoTwistD3Q27System::etBNE)
    //      (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1)= f[D3Q27System::DIR_DIR_PPM]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::etTSW) == EsoTwistD3Q27System::etTSW)
    //      (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::dMMP]; directionFlag=true;
    //   if ((direction & EsoTwistD3Q27System::REST) == EsoTwistD3Q27System::REST)
    //      (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::REST]; directionFlag=true;
    //#ifdef _DEBUG
    //   if(!directionFlag)UB_THROW( UbException(UB_EXARGS, "Direction didn't find") );
    //#endif //DEBUG
}
//////////////////////////////////////////////////////////////////////////
void D3Q27EsoTwist3DSoA::setPreCollisionDistributionForDirection(real f, size_t x1, size_t x2, size_t x3,
                                                        unsigned long int direction)
{
    // switch (direction)
    //{
    // case D3Q27System::dP00 :
    //   (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dM00 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f;
    //   break;
    // case D3Q27System::d0M0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f;
    //   break;
    // case D3Q27System::d0P0 :
    //   (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d00M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f;
    //   break;
    // case D3Q27System::d00P :
    //   (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dMM0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::dPP0 :
    //   (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dMP0 :
    //   (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::dPM0 :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f;
    //   break;
    // case D3Q27System::dM0M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::dP0P :
    //   (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::dM0P :
    //   (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::dP0M :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::d0MM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f;
    //   break;
    // case D3Q27System::d0PP :
    //   (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::d0MP :
    //   (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::d0PM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f;
    //   break;
    // case D3Q27System::BSW :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::dPPP :
    //   (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f;
    //   break;
    // case D3Q27System::BSE :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f;
    //   break;
    // case D3Q27System::dMPP :
    //   (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f;
    //   break;
    // case D3Q27System::BNW :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f;
    //   break;
    // case D3Q27System::dPMP :
    //   (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f;
    //   break;
    // case D3Q27System::DIR_DIR_PPM :
    //   (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f;
    //   break;
    // case D3Q27System::dMMP :
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
real D3Q27EsoTwist3DSoA::getPostCollisionDistributionForDirection(size_t /*x1*/, size_t /*x2*/, size_t /*x3*/,
                                                           int /*direction*/)
{
    // switch (direction)
    //{
    // case D3Q27System::dP00 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    );
    // case D3Q27System::dM00 :
    //   return (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3);
    // case D3Q27System::d0M0 :
    //   return (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3);
    // case D3Q27System::d0P0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    );
    // case D3Q27System::d00M :
    //   return (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3);
    // case D3Q27System::d00P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  );
    // case D3Q27System::dMM0 :
    //   return (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3);
    // case D3Q27System::dPP0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   );
    // case D3Q27System::dMP0 :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   );
    // case D3Q27System::dPM0 :
    //   return (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3);
    // case D3Q27System::dM0M :
    //   return (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3);
    // case D3Q27System::dP0P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 );
    // case D3Q27System::dM0P :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 );
    // case D3Q27System::dP0M :
    //   return (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3);
    // case D3Q27System::d0MM :
    //   return (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3);
    // case D3Q27System::d0PP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 );
    // case D3Q27System::d0MP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 );
    // case D3Q27System::d0PM :
    //   return (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3);
    // case D3Q27System::BSW :
    //   return (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3);
    // case D3Q27System::dPPP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
    // case D3Q27System::BSE :
    //   return (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3);
    // case D3Q27System::dMPP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1);
    // case D3Q27System::BNW :
    //   return (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3);
    // case D3Q27System::dPMP :
    //   return (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1);
    // case D3Q27System::DIR_DIR_PPM :
    //   return (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
    // case D3Q27System::dMMP :
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
