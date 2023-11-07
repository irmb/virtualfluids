#include "BGKLBMKernel.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "D3Q27EsoTwist3DSoA.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "Block3D.h"
#include "basics/constants/NumericConstants.h"

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
BGKLBMKernel::BGKLBMKernel() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
BGKLBMKernel::~BGKLBMKernel(void) = default;
//////////////////////////////////////////////////////////////////////////
void BGKLBMKernel::initDataSet()
{
    SPtr<DistributionArray3D> d(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> BGKLBMKernel::clone()
{
    SPtr<LBMKernel> kernel(new BGKLBMKernel());
    kernel->setNX(nx);
    std::dynamic_pointer_cast<BGKLBMKernel>(kernel)->initDataSet();
    kernel->setCollisionFactor(this->collFactor);
    kernel->setBCSet(bcSet->clone(kernel));
    kernel->setWithForcing(withForcing);
    kernel->setForcingX1(muForcingX1);
    kernel->setForcingX2(muForcingX2);
    kernel->setForcingX3(muForcingX3);
    kernel->setIndex(ix1, ix2, ix3);
    kernel->setDeltaT(deltaT);
    kernel->setBlock(block.lock());
    return kernel;
}
//////////////////////////////////////////////////////////////////////////
void BGKLBMKernel::calculate(int step)
{
    using namespace D3Q27System;
 //   using namespace UbMath;
   using namespace vf::basics::constant;
   using namespace vf::lbm::dir;

    // initializing of forcing stuff
    if (withForcing) {
        muForcingX1.DefineVar("x1", &muX1);
        muForcingX1.DefineVar("x2", &muX2);
        muForcingX1.DefineVar("x3", &muX3);
        muForcingX2.DefineVar("x1", &muX1);
        muForcingX2.DefineVar("x2", &muX2);
        muForcingX2.DefineVar("x3", &muX3);
        muForcingX3.DefineVar("x1", &muX1);
        muForcingX3.DefineVar("x2", &muX2);
        muForcingX3.DefineVar("x3", &muX3);
        forcingX1 = 0;
        forcingX2 = 0;
        forcingX3 = 0;
    }
    /////////////////////////////////////

    localDistributions =
        std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
    nonLocalDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())
                                ->getNonLocalDistributions();
    zeroDistributions =
        std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

    SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();
    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    real drho, vx1, vx2, vx3;
    const int bcArrayMaxX1 = (int)bcArray->getNX1();
    const int bcArrayMaxX2 = (int)bcArray->getNX2();
    const int bcArrayMaxX3 = (int)bcArray->getNX3();

    int minX1 = ghostLayerWidth;
    int minX2 = ghostLayerWidth;
    int minX3 = ghostLayerWidth;
    int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
    int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
    int maxX3 = bcArrayMaxX3 - ghostLayerWidth;

    for (int x3 = minX3; x3 < maxX3; x3++) {
        for (int x2 = minX2; x2 < maxX2; x2++) {
            for (int x1 = minX1; x1 < maxX1; x1++) {
                if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                    int x1p = x1 + 1;
                    int x2p = x2 + 1;
                    int x3p = x3 + 1;
                    //////////////////////////////////////////////////////////////////////////
                    // read distribution
                    ////////////////////////////////////////////////////////////////////////////
                    f[d000] = (*this->zeroDistributions)(x1, x2, x3);

                    f[dP00]   = (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
                    f[d0P0]   = (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
                    f[d00P]   = (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
                    f[dPP0]  = (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
                    f[dMP0]  = (*this->localDistributions)(D3Q27System::ET_NW, x1p, x2, x3);
                    f[dP0P]  = (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
                    f[dM0P]  = (*this->localDistributions)(D3Q27System::ET_TW, x1p, x2, x3);
                    f[d0PP]  = (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
                    f[d0MP]  = (*this->localDistributions)(D3Q27System::ET_TS, x1, x2p, x3);
                    f[dPPP] = (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
                    f[dMPP] = (*this->localDistributions)(D3Q27System::ET_TNW, x1p, x2, x3);
                    f[dPMP] = (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2p, x3);
                    f[dMMP] = (*this->localDistributions)(D3Q27System::ET_TSW, x1p, x2p, x3);

                    f[dM00]   = (*this->nonLocalDistributions)(D3Q27System::ET_W, x1p, x2, x3);
                    f[d0M0]   = (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2p, x3);
                    f[d00M]   = (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3p);
                    f[dMM0]  = (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1p, x2p, x3);
                    f[dPM0]  = (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2p, x3);
                    f[dM0M]  = (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1p, x2, x3p);
                    f[dP0M]  = (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3p);
                    f[d0MM]  = (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2p, x3p);
                    f[d0PM]  = (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3p);
                    f[dMMM] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                    f[dPMM] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2p, x3p);
                    f[dMPM] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1p, x2, x3p);
                    f[dPPM] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3p);
                    //////////////////////////////////////////////////////////////////////////

                    drho = f[d000] + f[dP00] + f[dM00] + f[d0P0] + f[d0M0] + f[d00P] + f[d00M] + f[dPP0] + f[dMM0] + f[dPM0] + f[dMP0] + f[dP0P] +
                           f[dM0M] + f[dP0M] + f[dM0P] + f[d0PP] + f[d0MM] + f[d0PM] + f[d0MP] + f[dPPP] + f[dMMP] + f[dPMP] + f[dMPP] +
                           f[dPPM] + f[dMMM] + f[dPMM] + f[dMPM];

                    vx1 = f[dP00] - f[dM00] + f[dPP0] - f[dMM0] + f[dPM0] - f[dMP0] + f[dP0P] - f[dM0M] + f[dP0M] - f[dM0P] + f[dPPP] -
                          f[dMMP] + f[dPMP] - f[dMPP] + f[dPPM] - f[dMMM] + f[dPMM] - f[dMPM];

                    vx2 = f[d0P0] - f[d0M0] + f[dPP0] - f[dMM0] - f[dPM0] + f[dMP0] + f[d0PP] - f[d0MM] + f[d0PM] - f[d0MP] + f[dPPP] -
                          f[dMMP] - f[dPMP] + f[dMPP] + f[dPPM] - f[dMMM] - f[dPMM] + f[dMPM];

                    vx3 = f[d00P] - f[d00M] + f[dP0P] - f[dM0M] - f[dP0M] + f[dM0P] + f[d0PP] - f[d0MM] - f[d0PM] + f[d0MP] + f[dPPP] +
                          f[dMMP] + f[dPMP] + f[dMPP] - f[dPPM] - f[dMMM] - f[dPMM] - f[dMPM];

                    real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

                    feq[d000] = c8o27 * (drho - cu_sq);
                    feq[dP00]    = c2o27 * (drho + c3o1 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq);
                    feq[dM00]    = c2o27 * (drho + c3o1 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq);
                    feq[d0P0]    = c2o27 * (drho + c3o1 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq);
                    feq[d0M0]    = c2o27 * (drho + c3o1 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq);
                    feq[d00P]    = c2o27 * (drho + c3o1 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq);
                    feq[d00M]    = c2o27 * (drho + c3o1 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq);
                    feq[dPP0]   = c1o54 * (drho + c3o1 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
                    feq[dMM0]   = c1o54 * (drho + c3o1 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
                    feq[dPM0]   = c1o54 * (drho + c3o1 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
                    feq[dMP0]   = c1o54 * (drho + c3o1 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
                    feq[dP0P]   = c1o54 * (drho + c3o1 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
                    feq[dM0M]   = c1o54 * (drho + c3o1 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
                    feq[dP0M]   = c1o54 * (drho + c3o1 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
                    feq[dM0P]   = c1o54 * (drho + c3o1 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
                    feq[d0PP]   = c1o54 * (drho + c3o1 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
                    feq[d0MM]   = c1o54 * (drho + c3o1 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
                    feq[d0PM]   = c1o54 * (drho + c3o1 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
                    feq[d0MP]   = c1o54 * (drho + c3o1 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
                    feq[dPPP]  = c1o216 *
                               (drho + c3o1 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
                    feq[dMMM] = c1o216 * (drho + c3o1 * (-vx1 - vx2 - vx3) +
                                         c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
                    feq[dPPM] = c1o216 *
                               (drho + c3o1 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
                    feq[dMMP] = c1o216 * (drho + c3o1 * (-vx1 - vx2 + vx3) +
                                         c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
                    feq[dPMP] = c1o216 *
                               (drho + c3o1 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
                    feq[dMPM] = c1o216 * (drho + c3o1 * (-vx1 + vx2 - vx3) +
                                         c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
                    feq[dPMM] = c1o216 *
                               (drho + c3o1 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
                    feq[dMPP] = c1o216 * (drho + c3o1 * (-vx1 + vx2 + vx3) +
                                         c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);

                    // Relaxation
                    f[d000] += (feq[d000] - f[d000]) * collFactor;
                    f[dP00] += (feq[dP00] - f[dP00]) * collFactor;
                    f[dM00] += (feq[dM00] - f[dM00]) * collFactor;
                    f[d0P0] += (feq[d0P0] - f[d0P0]) * collFactor;
                    f[d0M0] += (feq[d0M0] - f[d0M0]) * collFactor;
                    f[d00P] += (feq[d00P] - f[d00P]) * collFactor;
                    f[d00M] += (feq[d00M] - f[d00M]) * collFactor;
                    f[dPP0] += (feq[dPP0] - f[dPP0]) * collFactor;
                    f[dMM0] += (feq[dMM0] - f[dMM0]) * collFactor;
                    f[dPM0] += (feq[dPM0] - f[dPM0]) * collFactor;
                    f[dMP0] += (feq[dMP0] - f[dMP0]) * collFactor;
                    f[dP0P] += (feq[dP0P] - f[dP0P]) * collFactor;
                    f[dM0M] += (feq[dM0M] - f[dM0M]) * collFactor;
                    f[dP0M] += (feq[dP0M] - f[dP0M]) * collFactor;
                    f[dM0P] += (feq[dM0P] - f[dM0P]) * collFactor;
                    f[d0PP] += (feq[d0PP] - f[d0PP]) * collFactor;
                    f[d0MM] += (feq[d0MM] - f[d0MM]) * collFactor;
                    f[d0PM] += (feq[d0PM] - f[d0PM]) * collFactor;
                    f[d0MP] += (feq[d0MP] - f[d0MP]) * collFactor;

                    f[dPPP] += (feq[dPPP] - f[dPPP]) * collFactor;
                    f[dMMM] += (feq[dMMM] - f[dMMM]) * collFactor;
                    f[dPPM] += (feq[dPPM] - f[dPPM]) * collFactor;
                    f[dMMP] += (feq[dMMP] - f[dMMP]) * collFactor;
                    f[dPMP] += (feq[dPMP] - f[dPMP]) * collFactor;
                    f[dMPM] += (feq[dMPM] - f[dMPM]) * collFactor;
                    f[dPMM] += (feq[dPMM] - f[dPMM]) * collFactor;
                    f[dMPP] += (feq[dMPP] - f[dMPP]) * collFactor;

                    //////////////////////////////////////////////////////////////////////////
                    // forcing
                    if (withForcing) {
                        muX1 = x1 + ix1 * bcArrayMaxX1;
                        muX2 = x2 + ix2 * bcArrayMaxX2;
                        muX3 = x3 + ix3 * bcArrayMaxX3;

                        forcingX1 = muForcingX1.Eval();
                        forcingX2 = muForcingX2.Eval();
                        forcingX3 = muForcingX3.Eval();

                        f[d000] += c0o1;
                        f[dP00] += c3o1 * c2o27 * (forcingX1);
                        f[dM00] += c3o1 * c2o27 * (-forcingX1);
                        f[d0P0] += c3o1 * c2o27 * (forcingX2);
                        f[d0M0] += c3o1 * c2o27 * (-forcingX2);
                        f[d00P] += c3o1 * c2o27 * (forcingX3);
                        f[d00M] += c3o1 * c2o27 * (-forcingX3);
                        f[dPP0] += c3o1 * c1o54 * (forcingX1 + forcingX2);
                        f[dMM0] += c3o1 * c1o54 * (-forcingX1 - forcingX2);
                        f[dPM0] += c3o1 * c1o54 * (forcingX1 - forcingX2);
                        f[dMP0] += c3o1 * c1o54 * (-forcingX1 + forcingX2);
                        f[dP0P] += c3o1 * c1o54 * (forcingX1 + forcingX3);
                        f[dM0M] += c3o1 * c1o54 * (-forcingX1 - forcingX3);
                        f[dP0M] += c3o1 * c1o54 * (forcingX1 - forcingX3);
                        f[dM0P] += c3o1 * c1o54 * (-forcingX1 + forcingX3);
                        f[d0PP] += c3o1 * c1o54 * (forcingX2 + forcingX3);
                        f[d0MM] += c3o1 * c1o54 * (-forcingX2 - forcingX3);
                        f[d0PM] += c3o1 * c1o54 * (forcingX2 - forcingX3);
                        f[d0MP] += c3o1 * c1o54 * (-forcingX2 + forcingX3);
                        f[dPPP] += c3o1 * c1o216 * (forcingX1 + forcingX2 + forcingX3);
                        f[dMMM] += c3o1 * c1o216 * (-forcingX1 - forcingX2 - forcingX3);
                        f[dPPM] += c3o1 * c1o216 * (forcingX1 + forcingX2 - forcingX3);
                        f[dMMP] += c3o1 * c1o216 * (-forcingX1 - forcingX2 + forcingX3);
                        f[dPMP] += c3o1 * c1o216 * (forcingX1 - forcingX2 + forcingX3);
                        f[dMPM] += c3o1 * c1o216 * (-forcingX1 + forcingX2 - forcingX3);
                        f[dPMM] += c3o1 * c1o216 * (forcingX1 - forcingX2 - forcingX3);
                        f[dMPP] += c3o1 * c1o216 * (-forcingX1 + forcingX2 + forcingX3);
                    }
                    //////////////////////////////////////////////////////////////////////////
#ifdef PROOF_CORRECTNESS
                    real rho_post = f[d000] + f[dP00] + f[dM00] + f[d0P0] + f[d0M0] + f[d00P] + f[d00M] + f[dPP0] + f[dMM0] + f[dPM0] +
                                       f[dMP0] + f[dP0P] + f[dM0M] + f[dP0M] + f[dM0P] + f[d0PP] + f[d0MM] + f[d0PM] + f[d0MP] + f[dPPP] +
                                       f[dMMP] + f[dPMP] + f[dMPP] + f[dPPM] + f[dMMM] + f[dPMM] + f[dMPM];
                    real dif = drho - rho_post;
#ifdef SINGLEPRECISION
                    if (dif > 10.0E-7 || dif < -10.0E-7)
#else
                    if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
                    {
                      UB_THROW(UbException(UB_EXARGS, "rho="+UbSystem::toString(drho)+", rho_post="+UbSystem::toString(rho_post)
                         +" dif="+UbSystem::toString(dif)
                         +" rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)
                         +" in " + block.lock()->toString()+" step = "+UbSystem::toString(step)));
                    }
#endif
                    //////////////////////////////////////////////////////////////////////////
                    // write distribution
                    //////////////////////////////////////////////////////////////////////////
                    (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3)     = f[iP00];
                    (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3)     = f[INV_0P0];
                    (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3)     = f[INV_00P];
                    (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3)    = f[INV_PP0];
                    (*this->localDistributions)(D3Q27System::ET_NW, x1p, x2, x3)   = f[INV_MP0];
                    (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3)    = f[INV_P0P];
                    (*this->localDistributions)(D3Q27System::ET_TW, x1p, x2, x3)   = f[INV_M0P];
                    (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3)    = f[INV_0PP];
                    (*this->localDistributions)(D3Q27System::ET_TS, x1, x2p, x3)   = f[INV_0MP];
                    (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3)   = f[INV_PPP];
                    (*this->localDistributions)(D3Q27System::ET_TNW, x1p, x2, x3)  = f[INV_MPP];
                    (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2p, x3)  = f[INV_PMP];
                    (*this->localDistributions)(D3Q27System::ET_TSW, x1p, x2p, x3) = f[INV_MMP];

                    (*this->nonLocalDistributions)(D3Q27System::ET_W, x1p, x2, x3)     = f[INV_M00];
                    (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2p, x3)     = f[INV_0M0];
                    (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3p)     = f[INV_00M];
                    (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1p, x2p, x3)   = f[INV_MM0];
                    (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2p, x3)    = f[INV_PM0];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1p, x2, x3p)   = f[INV_M0M];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3p)    = f[INV_P0M];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2p, x3p)   = f[INV_0MM];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3p)    = f[INV_0PM];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1p, x2p, x3p) = f[INV_MMM];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2p, x3p)  = f[INV_PMM];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1p, x2, x3p)  = f[INV_MPM];
                    (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3p)   = f[INV_PPM];

                    (*this->zeroDistributions)(x1, x2, x3) = f[d000];
                    //////////////////////////////////////////////////////////////////////////
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
real BGKLBMKernel::getCalculationTime() { return vf::basics::constant::c0o1; }
