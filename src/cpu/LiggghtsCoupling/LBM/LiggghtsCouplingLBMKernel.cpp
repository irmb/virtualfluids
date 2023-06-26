#include "LiggghtsCouplingLBMKernel.h"
#include "D3Q27System.h"

//void LiggghtsCouplingLBMKernel::collisionOperator(int x1, int x2, int x3, real collFactorM, real fPre[])
//{
//    //if ((*particleData)(x1, x2, x3)->solidFraction >= SOLFRAC_MIN) {
//        LBMReal f[D3Q27System::ENDF + 1];
//        LBMReal fEq[D3Q27System::ENDF + 1];
//        LBMReal fEqSolid[D3Q27System::ENDF + 1];
//        LBMReal vx1, vx2, vx3, drho;
//        D3Q27System::calcIncompMacroscopicValues(f, drho, vx1, vx2, vx3);
//        D3Q27System::calcIncompFeq(fEq, drho, vx1, vx2, vx3);
//
//        std::array<double, 3> uPart;
//        uPart[0] = (*particleData)(x1, x2, x3)->uPart[0];
//        uPart[1] = (*particleData)(x1, x2, x3)->uPart[1];
//        uPart[2] = (*particleData)(x1, x2, x3)->uPart[2];
//
//        D3Q27System::calcIncompFeq(fEqSolid, drho, uPart[0], uPart[1], uPart[2]);
//        real rhoPhaseField = (phi[DIR_000] > c1o2) ? c1o1 : c1o1 / densityRatio;
//        if ((*particleData)(x1, x2, x3)->solidFraction > SOLFRAC_MAX) {
//            double const bb0 = fEq[vf::lbm::dir::DIR_000] - fEqSolid[vf::lbm::dir::DIR_000];
//            f[vf::lbm::dir::DIR_000] = fPre[vf::lbm::dir::DIR_000] + bb0;
//            for (int iPop = D3Q27System::FSTARTDIR; iPop <= D3Q27System::FENDDIR; iPop++) {
//                const int iOpp = D3Q27System::INVDIR[iPop];
//                double const bb = ((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop] - fEqSolid[iPop]));
//                double const bbOpp = ((fPre[iPop] - fEq[iPop]) - (fPre[iOpp] - fEqSolid[iOpp]));
//
//                f[iPop] = fPre[iPop] + bb;
//                f[iOpp] = fPre[iOpp] + bbOpp;
//
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[0] -= D3Q27System::DX1[iPop] * (bb - bbOpp) * rhoPhaseField;
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[1] -= D3Q27System::DX2[iPop] * (bb - bbOpp) * rhoPhaseField;
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[2] -= D3Q27System::DX3[iPop] * (bb - bbOpp) * rhoPhaseField;
//            }
//        } else { /* particleData.solidFraction < SOLFRAC_MAX */
//                 // #ifdef LBDEM_USE_WEIGHING
//            double const ooo = 1. / collFactorM - 0.5;
//            double const B = (*particleData)(x1, x2, x3)->solidFraction * ooo / ((1. - (*particleData)(x1, x2, x3)->solidFraction) + ooo);
//            // #else
//            //                         T const B = particleData.solidFraction;
//            // #endif
//            double const oneMinB = 1. - B;
//
//            double const bb0 = fEq[vf::lbm::dir::DIR_000] - fEqSolid[vf::lbm::dir::DIR_000];
//            f[vf::lbm::dir::DIR_000] = fPre[vf::lbm::dir::DIR_000] + oneMinB * (f[vf::lbm::dir::DIR_000] - fPre[vf::lbm::dir::DIR_000]) + B * bb0;
//
//            for (int iPop = D3Q27System::FSTARTDIR; iPop <= D3Q27System::FENDDIR; iPop++) {
//                int const iOpp = D3Q27System::INVDIR[iPop];
//                double const bb = B * ((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop] - fEqSolid[iPop]));
//                double const bbOpp = B * ((fPre[iPop] - fEq[iPop]) - (fPre[iOpp] - fEqSolid[iOpp]));
//
//                f[iPop] = fPre[iPop] + oneMinB * (f[iPop] - fPre[iPop]) + bb;
//                f[iOpp] = fPre[iOpp] + oneMinB * (f[iOpp] - fPre[iOpp]) + bbOpp;
//
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[0] -= D3Q27System::DX1[iPop] * (bb - bbOpp) * rhoPhaseField;
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[1] -= D3Q27System::DX2[iPop] * (bb - bbOpp) * rhoPhaseField;
//                (*particleData)(x1, x2, x3)->hydrodynamicForce[2] -= D3Q27System::DX3[iPop] * (bb - bbOpp) * rhoPhaseField;
//            }
//        } /* if solidFraction > SOLFRAC_MAX */
//
//    //    (*this->restDistributionsF)(x1, x2, x3) = f[vf::lbm::dir::DIR_000];
//
//    //    (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = f[vf::lbm::dir::DIR_M00];
//    //    (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = f[vf::lbm::dir::DIR_0M0];
//    //    (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = f[vf::lbm::dir::DIR_00M];
//    //    (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = f[vf::lbm::dir::DIR_MM0];
//    //    (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = f[vf::lbm::dir::DIR_PM0];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = f[vf::lbm::dir::DIR_M0M];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = f[vf::lbm::dir::DIR_P0M];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = f[vf::lbm::dir::DIR_0MM];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = f[vf::lbm::dir::DIR_0PM];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = f[vf::lbm::dir::DIR_MMM];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = f[vf::lbm::dir::DIR_PMM];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = f[vf::lbm::dir::DIR_MPM];
//    //    (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = f[vf::lbm::dir::DIR_PPM];
//
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = f[vf::lbm::dir::DIR_P00];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = f[vf::lbm::dir::DIR_0P0];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = f[vf::lbm::dir::DIR_00P];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = f[vf::lbm::dir::DIR_PP0];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = f[vf::lbm::dir::DIR_MP0];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = f[vf::lbm::dir::DIR_P0P];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = f[vf::lbm::dir::DIR_M0P];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = f[vf::lbm::dir::DIR_0PP];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = f[vf::lbm::dir::DIR_0MP];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = f[vf::lbm::dir::DIR_PPP];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = f[vf::lbm::dir::DIR_MPP];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = f[vf::lbm::dir::DIR_PMP];
//    //    (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = f[vf::lbm::dir::DIR_MMP];
//    //}
//}