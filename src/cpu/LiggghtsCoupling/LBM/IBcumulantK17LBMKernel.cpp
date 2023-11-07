//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file IBcumulantK17LBMKernel.cpp
//! \ingroup LiggghtsCoupling
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================
#include "IBcumulantK17LBMKernel.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <cmath>
#include "DataSet3D.h"
#include "LBMKernel.h"
#include "Block3D.h"
#include "BCArray3D.h"

#define PROOF_CORRECTNESS

//using namespace UbMath;

//////////////////////////////////////////////////////////////////////////
IBcumulantK17LBMKernel::IBcumulantK17LBMKernel()
{
    this->compressible = true;
    forcingX1 = 0.0;
    forcingX2 = 0.0;
    forcingX3 = 0.0;
}
//////////////////////////////////////////////////////////////////////////
void IBcumulantK17LBMKernel::initDataSet()
{
    SPtr<DistributionArray3D> d(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);

    particleData = std::make_shared<CbArray3D<SPtr<IBdynamicsParticleData>, IndexerX3X2X1>>(nx[0] + 2, nx[1] + 2, nx[2] + 2);

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;
    int maxX1 = nx[0]+2;
    int maxX2 = nx[1]+2;
    int maxX3 = nx[2]+2;

    LBMReal omega = collFactor;

    for (int x3 = minX3; x3 < maxX3; x3++) {
        for (int x2 = minX2; x2 < maxX2; x2++) {
            for (int x1 = minX1; x1 < maxX1; x1++) {
                (*particleData)(x1, x2, x3) = std::make_shared<IBdynamicsParticleData>();
            }
        }
    }

}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> IBcumulantK17LBMKernel::clone()
{
    SPtr<LBMKernel> kernel(new IBcumulantK17LBMKernel());
    kernel->setNX(nx);
    std::dynamic_pointer_cast<IBcumulantK17LBMKernel>(kernel)->initDataSet();
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
void IBcumulantK17LBMKernel::calculate(int step)
{
    //////////////////////////////////////////////////////////////////////////
    //! Cumulant K17 Kernel is based on
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //! and
    //! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
    //!
    //! The cumulant kernel is executed in the following steps
    //!
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from thredIdx, blockIdx, blockDim and gridDim.
    //!

    using namespace std;
    using namespace vf::basics::constant;

    //initializing of forcing stuff
    if (withForcing)
    {
        muForcingX1.DefineVar("x1", &muX1); muForcingX1.DefineVar("x2", &muX2); muForcingX1.DefineVar("x3", &muX3);
        muForcingX2.DefineVar("x1", &muX1); muForcingX2.DefineVar("x2", &muX2); muForcingX2.DefineVar("x3", &muX3);
        muForcingX3.DefineVar("x1", &muX1); muForcingX3.DefineVar("x2", &muX2); muForcingX3.DefineVar("x3", &muX3);

        muDeltaT = deltaT;

        muForcingX1.DefineVar("dt", &muDeltaT);
        muForcingX2.DefineVar("dt", &muDeltaT);
        muForcingX3.DefineVar("dt", &muDeltaT);

        muNu = (1.0 / 3.0) * (1.0 / collFactor - 1.0 / 2.0);

        muForcingX1.DefineVar("nu", &muNu);
        muForcingX2.DefineVar("nu", &muNu);
        muForcingX3.DefineVar("nu", &muNu);
    }
    /////////////////////////////////////

    localDistributionsF = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
    nonLocalDistributionsF = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
    restDistributionsF = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

    SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

    const int bcArrayMaxX1 = (int)bcArray->getNX1();
    const int bcArrayMaxX2 = (int)bcArray->getNX2();
    const int bcArrayMaxX3 = (int)bcArray->getNX3();

    int minX1 = ghostLayerWidth;
    int minX2 = ghostLayerWidth;
    int minX3 = ghostLayerWidth;
    int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
    int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
    int maxX3 = bcArrayMaxX3 - ghostLayerWidth;

    LBMReal omega = collFactor;



    for (int x3 = minX3; x3 < maxX3; x3++)
    {
        for (int x2 = minX2; x2 < maxX2; x2++)
        {
            for (int x1 = minX1; x1 < maxX1; x1++)
            {
                if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                    int x1p = x1 + 1;
                    int x2p = x2 + 1;
                    int x3p = x3 + 1;
                    //////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////
                    //! - Read distributions: style of reading and writing the distributions from/to stored arrays
                    //! dependent on timestep is based on the esoteric twist algorithm <a
                    //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
                    //! DOI:10.3390/computation5020019 ]</b></a>
                    //!
                    ////////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////

                    // E   N  T
                    // c   c  c
                    //////////
                    // W   S  B
                    // a   a  a

                    // Rest is b

                    // mfxyz
                    // a - negative
                    // b - null
                    // c - positive

                    // a b c
                    //-1 0 1

                    LBMReal mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
                    LBMReal mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
                    LBMReal mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
                    LBMReal mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
                    LBMReal mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
                    LBMReal mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
                    LBMReal mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
                    LBMReal mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
                    LBMReal mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
                    LBMReal mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
                    LBMReal mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
                    LBMReal mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
                    LBMReal mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);

                    LBMReal mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
                    LBMReal mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
                    LBMReal mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
                    LBMReal mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
                    LBMReal mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
                    LBMReal mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
                    LBMReal mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
                    LBMReal mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
                    LBMReal mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
                    LBMReal mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                    LBMReal mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
                    LBMReal mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
                    LBMReal mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

                    LBMReal mfbbb = (*this->restDistributionsF)(x1, x2, x3);

                    LBMReal f[D3Q27System::ENDF + 1];
                    LBMReal fEq[D3Q27System::ENDF + 1];
                    LBMReal fEqSolid[D3Q27System::ENDF + 1];
                    LBMReal fPre[D3Q27System::ENDF + 1];

                    f[vf::lbm::dir::d000] = mfbbb;

                    f[vf::lbm::dir::dP00] = mfcbb;
                    f[vf::lbm::dir::DIR_0P0] = mfbcb;
                    f[vf::lbm::dir::DIR_00P] = mfbbc;
                    f[vf::lbm::dir::DIR_PP0] = mfccb;
                    f[vf::lbm::dir::DIR_MP0] = mfacb;
                    f[vf::lbm::dir::DIR_P0P] = mfcbc;
                    f[vf::lbm::dir::DIR_M0P] = mfabc;
                    f[vf::lbm::dir::DIR_0PP] = mfbcc;
                    f[vf::lbm::dir::DIR_0MP] = mfbac;
                    f[vf::lbm::dir::DIR_PPP] = mfccc;
                    f[vf::lbm::dir::DIR_MPP] = mfacc;
                    f[vf::lbm::dir::DIR_PMP] = mfcac;
                    f[vf::lbm::dir::DIR_MMP] = mfaac;

                    f[vf::lbm::dir::dM00] = mfabb;
                    f[vf::lbm::dir::DIR_0M0] = mfbab;
                    f[vf::lbm::dir::DIR_00M] = mfbba;
                    f[vf::lbm::dir::DIR_MM0] = mfaab;
                    f[vf::lbm::dir::DIR_PM0] = mfcab;
                    f[vf::lbm::dir::DIR_M0M] = mfaba;
                    f[vf::lbm::dir::DIR_P0M] = mfcba;
                    f[vf::lbm::dir::DIR_0MM] = mfbaa;
                    f[vf::lbm::dir::DIR_0PM] = mfbca;
                    f[vf::lbm::dir::DIR_MMM] = mfaaa;
                    f[vf::lbm::dir::DIR_PMM] = mfcaa;
                    f[vf::lbm::dir::DIR_MPM] = mfaca;
                    f[vf::lbm::dir::DIR_PPM] = mfcca;

                    if ((*particleData)(x1, x2, x3)->solidFraction > SOLFRAC_MIN) {
                        fPre[vf::lbm::dir::d000] = mfbbb;

                        fPre[vf::lbm::dir::dP00] = mfcbb;
                        fPre[vf::lbm::dir::DIR_0P0] = mfbcb;
                        fPre[vf::lbm::dir::DIR_00P] = mfbbc;
                        fPre[vf::lbm::dir::DIR_PP0] = mfccb;
                        fPre[vf::lbm::dir::DIR_MP0] = mfacb;
                        fPre[vf::lbm::dir::DIR_P0P] = mfcbc;
                        fPre[vf::lbm::dir::DIR_M0P] = mfabc;
                        fPre[vf::lbm::dir::DIR_0PP] = mfbcc;
                        fPre[vf::lbm::dir::DIR_0MP] = mfbac;
                        fPre[vf::lbm::dir::DIR_PPP] = mfccc;
                        fPre[vf::lbm::dir::DIR_MPP] = mfacc;
                        fPre[vf::lbm::dir::DIR_PMP] = mfcac;
                        fPre[vf::lbm::dir::DIR_MMP] = mfaac;
                          
                        fPre[vf::lbm::dir::dM00] = mfabb;
                        fPre[vf::lbm::dir::DIR_0M0] = mfbab;
                        fPre[vf::lbm::dir::DIR_00M] = mfbba;
                        fPre[vf::lbm::dir::DIR_MM0] = mfaab;
                        fPre[vf::lbm::dir::DIR_PM0] = mfcab;
                        fPre[vf::lbm::dir::DIR_M0M] = mfaba;
                        fPre[vf::lbm::dir::DIR_P0M] = mfcba;
                        fPre[vf::lbm::dir::DIR_0MM] = mfbaa;
                        fPre[vf::lbm::dir::DIR_0PM] = mfbca;
                        fPre[vf::lbm::dir::DIR_MMM] = mfaaa;
                        fPre[vf::lbm::dir::DIR_PMM] = mfcaa;
                        fPre[vf::lbm::dir::DIR_MPM] = mfaca;
                        fPre[vf::lbm::dir::DIR_PPM] = mfcca;
                    }

                    (*particleData)(x1, x2, x3)->hydrodynamicForce.fill(0.0);

                    if ((*particleData)(x1, x2, x3)->solidFraction <= SOLFRAC_MAX) {

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq.
                    //! (J1)-(J3) <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
                    //!
                    LBMReal drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                                    (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) +
                                     ((mfacb + mfcab) + (mfaab + mfccb))) +
                                    ((mfabb + mfcbb) + (mfbab + mfbcb)) + (mfbba + mfbbc)) +
                                   mfbbb;

                    LBMReal rho   = c1o1 + drho;
                    LBMReal OOrho = c1o1 / rho;
                    ////////////////////////////////////////////////////////////////////////////////////
                    LBMReal vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                                   (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                                   (mfcbb - mfabb)) /
                                  rho;
                    LBMReal vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                                   (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                                   (mfbcb - mfbab)) /
                                  rho;
                    LBMReal vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                                   (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                                   (mfbbc - mfbba)) /
                                  rho;
                    ////////////////////////////////////////////////////////////////////////////////////
                    // forcing
                    ///////////////////////////////////////////////////////////////////////////////////////////
                    if (withForcing) {
                        muX1 = static_cast<double>(x1 - 1 + ix1 * maxX1);
                        muX2 = static_cast<double>(x2 - 1 + ix2 * maxX2);
                        muX3 = static_cast<double>(x3 - 1 + ix3 * maxX3);

                        forcingX1 = muForcingX1.Eval();
                        forcingX2 = muForcingX2.Eval();
                        forcingX3 = muForcingX3.Eval();

                        ////////////////////////////////////////////////////////////////////////////////////
                        //! - Add half of the acceleration (body force) to the velocity as in Eq. (42)
                        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
                        //!
                        vvx += forcingX1 * deltaT * c1o2; // X
                        vvy += forcingX2 * deltaT * c1o2; // Y
                        vvz += forcingX3 * deltaT * c1o2; // Z
                    }
                    ////////////////////////////////////////////////////////////////////////////////////
                    // calculate the square of velocities for this lattice node
                    LBMReal vx2 = vvx * vvx;
                    LBMReal vy2 = vvy * vvy;
                    LBMReal vz2 = vvz * vvz;
                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Set relaxation limiters for third order cumulants to default value \f$ \lambda=0.001 \f$
                    //! according to section 6 in <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et
                    //! al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    LBMReal wadjust;
                    LBMReal qudricLimitP = c1o100;
                    LBMReal qudricLimitM = c1o100;
                    LBMReal qudricLimitD = c1o100;
                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Chimera transform from well conditioned distributions to central moments as defined in
                    //! Appendix J in <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (6)-(14) in <a
                    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Z - Dir
                    forwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
                    forwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9);
                    forwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
                    forwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9);
                    forwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9);
                    forwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9);
                    forwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
                    forwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9);
                    forwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);

                    ////////////////////////////////////////////////////////////////////////////////////
                    // Y - Dir
                    forwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6);
                    forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
                    forwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
                    forwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3);
                    forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
                    forwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9);
                    forwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6);
                    forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
                    forwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);

                    ////////////////////////////////////////////////////////////////////////////////////
                    // X - Dir
                    forwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
                    forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
                    forwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
                    forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
                    forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
                    forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
                    forwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
                    forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
                    forwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c9o1, c1o9);

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Setting relaxation rates for non-hydrodynamic cumulants (default values). Variable names and
                    //! equations according to <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al.
                    //! (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!  => [NAME IN PAPER]=[NAME IN CODE]=[DEFAULT VALUE].
                    //!  - Trace of second order cumulants \f$ C_{200}+C_{020}+C_{002} \f$ used to adjust bulk
                    //!  viscosity:\f$\omega_2=OxxPyyPzz=1.0 \f$.
                    //!  - Third order cumulants \f$ C_{120}+C_{102} \f$, \f$ C_{210}+C_{012} \f$, \f$ C_{201}+C_{021}
                    //!  \f$: \f$\omega_3=OxyyPxzz\f$ set according to Eq. (111) with simplifications assuming
                    //!  \f$\omega_2=1.0\f$.
                    //!  - Third order cumulants \f$ C_{120}-C_{102} \f$, \f$ C_{210}-C_{012} \f$, \f$ C_{201}-C_{021}
                    //!  \f$: \f$\omega_4 = OxyyMxzz\f$ set according to Eq. (112) with simplifications assuming
                    //!  \f$\omega_2 = 1.0\f$.
                    //!  - Third order cumulants \f$ C_{111} \f$: \f$\omega_5 = Oxyz\f$ set according to Eq. (113) with
                    //!  simplifications assuming \f$\omega_2 = 1.0\f$  (modify for different bulk viscosity).
                    //!  - Fourth order cumulants \f$ C_{220} \f$, \f$ C_{202} \f$, \f$ C_{022} \f$, \f$ C_{211} \f$,
                    //!  \f$ C_{121} \f$, \f$ C_{112} \f$: for simplification all set to the same default value \f$
                    //!  \omega_6=\omega_7=\omega_8=O4=1.0 \f$.
                    //!  - Fifth order cumulants \f$ C_{221}\f$, \f$C_{212}\f$, \f$C_{122}\f$: \f$\omega_9=O5=1.0\f$.
                    //!  - Sixth order cumulant \f$ C_{222}\f$: \f$\omega_{10}=O6=1.0\f$.
                    //!
                    ////////////////////////////////////////////////////////////
                    // 2.
                    LBMReal OxxPyyPzz = c1o1;
                    ////////////////////////////////////////////////////////////
                    // 3.
                    LBMReal OxyyPxzz =
                        c8o1 * (-c2o1 + omega) * (c1o1 + c2o1 * omega) / (-c8o1 - c14o1 * omega + c7o1 * omega * omega);
                    LBMReal OxyyMxzz =
                        c8o1 * (-c2o1 + omega) * (-c7o1 + c4o1 * omega) / (c56o1 - c50o1 * omega + c9o1 * omega * omega);
                    LBMReal Oxyz = c24o1 * (-c2o1 + omega) * (-c2o1 - c7o1 * omega + c3o1 * omega * omega) /
                                   (c48o1 + c152o1 * omega - c130o1 * omega * omega + c29o1 * omega * omega * omega);
                    ////////////////////////////////////////////////////////////
                    // 4.
                    LBMReal O4 = c1o1;
                    ////////////////////////////////////////////////////////////
                    // 5.
                    LBMReal O5 = c1o1;
                    ////////////////////////////////////////////////////////////
                    // 6.
                    LBMReal O6 = c1o1;

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - A and B: parameters for fourth order convergence of the diffusion term according to Eq. (114)
                    //! and (115) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a> with simplifications assuming \f$\omega_2 = 1.0\f$
                    //! (modify for different bulk viscosity).
                    //!
                    LBMReal A = (c4o1 + c2o1 * omega - c3o1 * omega * omega) / (c2o1 - c7o1 * omega + c5o1 * omega * omega);
                    LBMReal B = (c4o1 + c28o1 * omega - c14o1 * omega * omega) / (c6o1 - c21o1 * omega + c15o1 * omega * omega);

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Compute cumulants from central moments according to Eq. (20)-(23) in
                    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    ////////////////////////////////////////////////////////////
                    // 4.
                    LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) * OOrho;
                    LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) * OOrho;
                    LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) * OOrho;

                    LBMReal CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) * OOrho -
                                              c1o9 * (drho * OOrho));
                    LBMReal CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) * OOrho -
                                              c1o9 * (drho * OOrho));
                    LBMReal CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) * OOrho -
                                              c1o9 * (drho * OOrho));
                    ////////////////////////////////////////////////////////////
                    // 5.
                    LBMReal CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) +
                         c1o3 * (mfbca + mfbac)) * OOrho;
                    LBMReal CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) +
                         c1o3 * (mfcba + mfabc)) * OOrho;
                    LBMReal CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) +
                         c1o3 * (mfacb + mfcab)) * OOrho;
                    ////////////////////////////////////////////////////////////
                    // 6.
                    LBMReal CUMccc =
                        mfccc + ((-c4o1 * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                                  c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
                                  c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) *
                                     OOrho +
                                 (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                                  c2o1 * (mfcaa * mfaca * mfaac) + c16o1 * mfbba * mfbab * mfabb) *
                                     OOrho * OOrho -
                                 c1o3 * (mfacc + mfcac + mfcca) * OOrho - c1o9 * (mfcaa + mfaca + mfaac) * OOrho +
                                 (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                                  (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 * (mfaac + mfaca + mfcaa)) *
                                     OOrho * OOrho * c2o3 +
                                 c1o27 * ((drho * drho - drho) * OOrho * OOrho));

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Compute linear combinations of second and third order cumulants
                    //!
                    ////////////////////////////////////////////////////////////
                    // 2.
                    LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
                    LBMReal mxxMyy    = mfcaa - mfaca;
                    LBMReal mxxMzz    = mfcaa - mfaac;
                    ////////////////////////////////////////////////////////////
                    // 3.
                    LBMReal mxxyPyzz = mfcba + mfabc;
                    LBMReal mxxyMyzz = mfcba - mfabc;

                    LBMReal mxxzPyyz = mfcab + mfacb;
                    LBMReal mxxzMyyz = mfcab - mfacb;

                    LBMReal mxyyPxzz = mfbca + mfbac;
                    LBMReal mxyyMxzz = mfbca - mfbac;

                    ////////////////////////////////////////////////////////////////////////////////////
                    // incl. correction
                    ////////////////////////////////////////////////////////////
                    //! - Compute velocity  gradients from second order cumulants according to Eq. (27)-(32)
                    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a> Further explanations of the correction in viscosity in
                    //! Appendix H of <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> Note that the division by rho is omitted here as we
                    //! need rho times the gradients later.
                    //!
                    LBMReal Dxy  = -c3o1 * omega * mfbba;
                    LBMReal Dxz  = -c3o1 * omega * mfbab;
                    LBMReal Dyz  = -c3o1 * omega * mfabb;
                    LBMReal dxux = c1o2 * (-omega) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
                    LBMReal dyuy = dxux + omega * c3o2 * mxxMyy;
                    LBMReal dzuz = dxux + omega * c3o2 * mxxMzz;
                    ////////////////////////////////////////////////////////////
                    //! - Relaxation of second order cumulants with correction terms according to Eq. (33)-(35) in
                    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    mxxPyyPzz += OxxPyyPzz * (mfaaa - mxxPyyPzz) -
                                 c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
                    mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
                    mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////no correction
                    // mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);
                    // mxxMyy += -(-omega) * (-mxxMyy);
                    // mxxMzz += -(-omega) * (-mxxMzz);
                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    mfabb += omega * (-mfabb);
                    mfbab += omega * (-mfbab);
                    mfbba += omega * (-mfbba);

                    ////////////////////////////////////////////////////////////////////////////////////
                    // relax
                    //////////////////////////////////////////////////////////////////////////
                    // incl. limiter
                    //! - Relaxation of third order cumulants including limiter according to Eq. (116)-(123)
                    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    wadjust = Oxyz + (c1o1 - Oxyz) * abs(mfbbb) / (abs(mfbbb) + qudricLimitD);
                    mfbbb += wadjust * (-mfbbb);
                    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimitP);
                    mxxyPyzz += wadjust * (-mxxyPyzz);
                    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimitM);
                    mxxyMyzz += wadjust * (-mxxyMyzz);
                    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimitP);
                    mxxzPyyz += wadjust * (-mxxzPyyz);
                    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimitM);
                    mxxzMyyz += wadjust * (-mxxzMyyz);
                    wadjust = OxyyPxzz + (c1o1 - OxyyPxzz) * abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimitP);
                    mxyyPxzz += wadjust * (-mxyyPxzz);
                    wadjust = OxyyMxzz + (c1o1 - OxyyMxzz) * abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimitM);
                    mxyyMxzz += wadjust * (-mxyyMxzz);
                    //////////////////////////////////////////////////////////////////////////
                    // no limiter
                    // mfbbb += OxyyMxzz * (-mfbbb);
                    // mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
                    // mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
                    // mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
                    // mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
                    // mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
                    // mxyyMxzz += OxyyMxzz * (-mxyyMxzz);

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Compute inverse linear combinations of second and third order cumulants
                    //!
                    mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
                    mfaca = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz);
                    mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

                    mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
                    mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
                    mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
                    mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
                    mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
                    mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
                    //////////////////////////////////////////////////////////////////////////

                    //////////////////////////////////////////////////////////////////////////
                    // 4.
                    //  no limiter
                    //! - Relax fourth order cumulants to modified equilibrium for fourth order convergence of diffusion
                    //! according to Eq. (43)-(48) <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et
                    //! al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    CUMacc = -O4 * (c1o1 / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMacc);
                    CUMcac = -O4 * (c1o1 / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMcac);
                    CUMcca = -O4 * (c1o1 / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (c1o1 - O4) * (CUMcca);
                    CUMbbc = -O4 * (c1o1 / omega - c1o2) * Dxy * c1o3 * B + (c1o1 - O4) * (CUMbbc);
                    CUMbcb = -O4 * (c1o1 / omega - c1o2) * Dxz * c1o3 * B + (c1o1 - O4) * (CUMbcb);
                    CUMcbb = -O4 * (c1o1 / omega - c1o2) * Dyz * c1o3 * B + (c1o1 - O4) * (CUMcbb);

                    //////////////////////////////////////////////////////////////////////////
                    // 5.
                    CUMbcc += O5 * (-CUMbcc);
                    CUMcbc += O5 * (-CUMcbc);
                    CUMccb += O5 * (-CUMccb);

                    //////////////////////////////////////////////////////////////////////////
                    // 6.
                    CUMccc += O6 * (-CUMccc);

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Compute central moments from post collision cumulants according to Eq. (53)-(56) in
                    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!

                    //////////////////////////////////////////////////////////////////////////
                    // 4.
                    mfcbb = CUMcbb + c1o3 * ((c3o1 * mfcaa + c1o1) * mfabb + c6o1 * mfbba * mfbab) * OOrho;
                    mfbcb = CUMbcb + c1o3 * ((c3o1 * mfaca + c1o1) * mfbab + c6o1 * mfbba * mfabb) * OOrho;
                    mfbbc = CUMbbc + c1o3 * ((c3o1 * mfaac + c1o1) * mfbba + c6o1 * mfbab * mfabb) * OOrho;

                    mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) * c9o1 + c3o1 * (mfcaa + mfaca)) * OOrho -
                                      (drho * OOrho)) *
                                         c1o9;
                    mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) * c9o1 + c3o1 * (mfcaa + mfaac)) * OOrho -
                                      (drho * OOrho)) *
                                         c1o9;
                    mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) * c9o1 + c3o1 * (mfaac + mfaca)) * OOrho -
                                      (drho * OOrho)) *
                                         c1o9;

                    //////////////////////////////////////////////////////////////////////////
                    // 5.
                    mfbcc = CUMbcc + c1o3 *
                                         (c3o1 * (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb +
                                                c2o1 * (mfbab * mfacb + mfbba * mfabc)) +
                                          (mfbca + mfbac)) *
                                         OOrho;
                    mfcbc = CUMcbc + c1o3 *
                                         (c3o1 * (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb +
                                                c2o1 * (mfabb * mfcab + mfbba * mfbac)) +
                                          (mfcba + mfabc)) *
                                         OOrho;
                    mfccb = CUMccb + c1o3 *
                                         (c3o1 * (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb +
                                                c2o1 * (mfbab * mfbca + mfabb * mfcba)) +
                                          (mfacb + mfcab)) *
                                         OOrho;

                    //////////////////////////////////////////////////////////////////////////
                    // 6.
                    mfccc =
                        CUMccc - ((-c4o1 * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                                   c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
                                   c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) *
                                      OOrho +
                                  (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                                   c2o1 * (mfcaa * mfaca * mfaac) + c16o1 * mfbba * mfbab * mfabb) *
                                      OOrho * OOrho -
                                  c1o3 * (mfacc + mfcac + mfcca) * OOrho - c1o9 * (mfcaa + mfaca + mfaac) * OOrho +
                                  (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                                   (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 * (mfaac + mfaca + mfcaa)) *
                                      OOrho * OOrho * c2o3 +
                                  c1o27 * ((drho * drho - drho) * OOrho * OOrho));

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! -  Add acceleration (body force) to first order cumulants according to Eq. (85)-(87) in
                    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
                    //!
                    mfbaa = -mfbaa;
                    mfaba = -mfaba;
                    mfaab = -mfaab;
                    ////////////////////////////////////////////////////////////////////////////////////

                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Chimera transform from central moments to well conditioned distributions as defined in
                    //! Appendix J in <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
                    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
                    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017),
                    //! DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
                    //!
                    ////////////////////////////////////////////////////////////////////////////////////
                    // X - Dir
                    backwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
                    backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
                    backwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3);
                    backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
                    backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
                    backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
                    backwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3);
                    backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
                    backwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c9o1, c1o9);

                    ////////////////////////////////////////////////////////////////////////////////////
                    // Y - Dir
                    backwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6);
                    backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
                    backwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18);
                    backwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3);
                    backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
                    backwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9);
                    backwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6);
                    backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
                    backwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18);

                    ////////////////////////////////////////////////////////////////////////////////////
                    // Z - Dir
                    backwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36);
                    backwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9);
                    backwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36);
                    backwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9);
                    backwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9);
                    backwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9);
                    backwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36);
                    backwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9);
                    backwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36);
                    ////////////////////////////////////////////////////////////////////////////////////

                    //////////////////////////////////////////////////////////////////////////
                    // proof correctness
                    //////////////////////////////////////////////////////////////////////////
#ifdef PROOF_CORRECTNESS
                    LBMReal drho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
                                        (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
                                        (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) +
                                        (mfbba + mfbbc) + mfbbb;
                    LBMReal dif = drho - drho_post;
#ifdef SINGLEPRECISION
                    if (dif > 10.0E-7 || dif < -10.0E-7)
#else
                    if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
                    {
                        UB_THROW(UbException(
                            UB_EXARGS,
                            "rho=" + UbSystem::toString(drho) + ", rho_post=" + UbSystem::toString(drho_post) +
                                " dif=" + UbSystem::toString(dif) + " rho is not correct for node " +
                                UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3) +
                                " in " + block.lock()->toString() + " step = " + UbSystem::toString(step)));
                    }
#endif
                    ////////////////////////////////////////////////////////////////////////////////////
                    //! - Write distributions: style of reading and writing the distributions from/to stored arrays
                    //! dependent on timestep is based on the esoteric twist algorithm <a
                    //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
                    //! DOI:10.3390/computation5020019 ]</b></a>
                    //!
                    (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3)     = mfabb;
                    (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3)     = mfbab;
                    (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)     = mfbba;
                    (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3)    = mfaab;
                    (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3)   = mfcab;
                    (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3)    = mfaba;
                    (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3)   = mfcba;
                    (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3)    = mfbaa;
                    (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3)   = mfbca;
                    (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3)   = mfaaa;
                    (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3)  = mfcaa;
                    (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3)  = mfaca;
                    (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

                    (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)     = mfcbb;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3)     = mfbcb;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p)     = mfbbc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfccb;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3)    = mfacb;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfcbc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p)    = mfabc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfbcc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p)    = mfbac;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfacc;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfcac;
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfaac;

                    (*this->restDistributionsF)(x1, x2, x3) = mfbbb;
                    //////////////////////////////////////////////////////////////////////////
                    f[vf::lbm::dir::d000] = mfbbb;
                     
                    f[vf::lbm::dir::dP00]  = mfcbb;
                    f[vf::lbm::dir::DIR_0P0]  = mfbcb;
                    f[vf::lbm::dir::DIR_00P]  = mfbbc;
                    f[vf::lbm::dir::DIR_PP0]  = mfccb;
                    f[vf::lbm::dir::DIR_MP0]  = mfacb;
                    f[vf::lbm::dir::DIR_P0P]  = mfcbc;
                    f[vf::lbm::dir::DIR_M0P]  = mfabc;
                    f[vf::lbm::dir::DIR_0PP]  = mfbcc;
                    f[vf::lbm::dir::DIR_0MP]  = mfbac;
                    f[vf::lbm::dir::DIR_PPP]  = mfccc;
                    f[vf::lbm::dir::DIR_MPP]  = mfacc;
                    f[vf::lbm::dir::DIR_PMP]  = mfcac;
                    f[vf::lbm::dir::DIR_MMP]  = mfaac;
                                     
                    f[vf::lbm::dir::dM00]  = mfabb;
                    f[vf::lbm::dir::DIR_0M0]  = mfbab;
                    f[vf::lbm::dir::DIR_00M]  = mfbba;
                    f[vf::lbm::dir::DIR_MM0]  = mfaab;
                    f[vf::lbm::dir::DIR_PM0]  = mfcab;
                    f[vf::lbm::dir::DIR_M0M]  = mfaba;
                    f[vf::lbm::dir::DIR_P0M]  = mfcba;
                    f[vf::lbm::dir::DIR_0MM]  = mfbaa;
                    f[vf::lbm::dir::DIR_0PM]  = mfbca;
                    f[vf::lbm::dir::DIR_MMM]  = mfaaa;
                    f[vf::lbm::dir::DIR_PMM]  = mfcaa;
                    f[vf::lbm::dir::DIR_MPM]  = mfaca;
                    f[vf::lbm::dir::DIR_PPM]  = mfcca;
                }
                    if ((*particleData)(x1, x2, x3)->solidFraction >= SOLFRAC_MIN) {
                    //                    if ((*particleData)(x1, x2, x3)->solidFraction < SOLFRAC_MIN)
                    //                        continue;
                    //
                                        LBMReal vx1, vx2, vx3, drho;
                                        D3Q27System::calcCompMacroscopicValues(f, drho, vx1, vx2, vx3);
                                        D3Q27System::calcCompFeq(fEq, drho, vx1, vx2, vx3);
                    
                                        std::array<double, 3> uPart;
                                        uPart[0] = (*particleData)(x1, x2, x3)->uPart[0] * (1. + drho);
                                        uPart[1] = (*particleData)(x1, x2, x3)->uPart[1] * (1. + drho);
                                        uPart[2] = (*particleData)(x1, x2, x3)->uPart[2] * (1. + drho);
                    
                                        D3Q27System::calcCompFeq(fEqSolid, drho, uPart[0], uPart[1], uPart[2]);
                    
                                        if ((*particleData)(x1, x2, x3)->solidFraction > SOLFRAC_MAX) {
                                        double const bb0 = fEq[vf::lbm::dir::d000] - fEqSolid[vf::lbm::dir::d000];
                                        f[vf::lbm::dir::d000] = fPre[vf::lbm::dir::d000] + bb0;
                                            for (int iPop = D3Q27System::FSTARTDIR; iPop <= D3Q27System::FENDDIR; iPop++) {
                                                const int iOpp        = D3Q27System::INVDIR[iPop];
                                                double const bb       = ((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop] - fEqSolid[iPop]));
                                                double const bbOpp    = ((fPre[iPop] - fEq[iPop]) - (fPre[iOpp] - fEqSolid[iOpp]));
                    
                    
                                                f[iPop] = fPre[iPop] + bb;
                                                f[iOpp] = fPre[iOpp] + bbOpp;
                    
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[0] -= D3Q27System::DX1[iPop] * (bb - bbOpp);
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[1] -= D3Q27System::DX2[iPop] * (bb - bbOpp);
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[2] -= D3Q27System::DX3[iPop] * (bb - bbOpp);
                                            }
                                        } else { /* particleData.solidFraction < SOLFRAC_MAX */
                    //#ifdef LBDEM_USE_WEIGHING
                                            double const ooo = 1. / omega - 0.5;
                                            double const B   = (*particleData)(x1, x2, x3)->solidFraction * ooo / ((1. - (*particleData)(x1, x2, x3)->solidFraction) + ooo);
                    //#else
                    //                        T const B = particleData.solidFraction;
                    //#endif
                                            double const oneMinB = 1. - B;
                    
                                            double const bb0 = fEq[vf::lbm::dir::d000] - fEqSolid[vf::lbm::dir::d000];
                                            f[vf::lbm::dir::d000] = fPre[vf::lbm::dir::d000] + oneMinB * (f[vf::lbm::dir::d000] - fPre[vf::lbm::dir::d000]) + B * bb0;
                    
                                            for (int iPop = D3Q27System::FSTARTDIR; iPop <= D3Q27System::FENDDIR; iPop++) {
                                                int const iOpp = D3Q27System::INVDIR[iPop];
                                                double const bb       = B * ((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop] - fEqSolid[iPop]));
                                                double const bbOpp    = B * ((fPre[iPop] - fEq[iPop]) - (fPre[iOpp] - fEqSolid[iOpp]));
                    
                                                f[iPop] = fPre[iPop] + oneMinB * (f[iPop] - fPre[iPop]) + bb;
                                                f[iOpp] = fPre[iOpp] + oneMinB * (f[iOpp] - fPre[iOpp]) + bbOpp;
                    
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[0] -= D3Q27System::DX1[iPop] * (bb - bbOpp);
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[1] -= D3Q27System::DX2[iPop] * (bb - bbOpp);
                                                (*particleData)(x1, x2, x3)->hydrodynamicForce[2] -= D3Q27System::DX3[iPop] * (bb - bbOpp);
                                            }
                                        } /* if solidFraction > SOLFRAC_MAX */

                    (*this->restDistributionsF)(x1, x2, x3) = f[vf::lbm::dir::d000];

                    (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = f[vf::lbm::dir::dM00];
                    (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = f[vf::lbm::dir::DIR_0M0];
                    (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = f[vf::lbm::dir::DIR_00M];
                    (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = f[vf::lbm::dir::DIR_MM0];
                    (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = f[vf::lbm::dir::DIR_PM0];
                    (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = f[vf::lbm::dir::DIR_M0M];
                    (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = f[vf::lbm::dir::DIR_P0M];
                    (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = f[vf::lbm::dir::DIR_0MM];
                    (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = f[vf::lbm::dir::DIR_0PM];
                    (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = f[vf::lbm::dir::DIR_MMM];
                    (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = f[vf::lbm::dir::DIR_PMM];
                    (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = f[vf::lbm::dir::DIR_MPM];
                    (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = f[vf::lbm::dir::DIR_PPM];

                    (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = f[vf::lbm::dir::dP00];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = f[vf::lbm::dir::DIR_0P0];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = f[vf::lbm::dir::DIR_00P];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = f[vf::lbm::dir::DIR_PP0];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = f[vf::lbm::dir::DIR_MP0];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = f[vf::lbm::dir::DIR_P0P];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = f[vf::lbm::dir::DIR_M0P];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = f[vf::lbm::dir::DIR_0PP];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = f[vf::lbm::dir::DIR_0MP];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = f[vf::lbm::dir::DIR_PPP];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = f[vf::lbm::dir::DIR_MPP];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = f[vf::lbm::dir::DIR_PMP];
                    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = f[vf::lbm::dir::DIR_MMP];
                    }
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////

