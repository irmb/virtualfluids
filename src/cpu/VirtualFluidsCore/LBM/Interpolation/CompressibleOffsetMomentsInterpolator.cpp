#include "CompressibleOffsetMomentsInterpolator.h"

#include <algorithm>

#include <basics/constants/NumericConstants.h>

#include <lbm/refinement/InterpolationCF.h>
#include <lbm/refinement/InterpolationFC.h>
#include <lbm/interpolation/InterpolationCoefficients.h>


void calculateCoefficients(vf::lbm::InterpolationCoefficients& coefficients, const D3Q27ICell& icell, real omega, real xoff, real yoff, real zoff)
{
    vf::lbm::MomentsOnSourceNodeSet momentsSet;

    momentsSet.calculateMMM(icell.BSW, omega);
    momentsSet.calculateMMP(icell.TSW, omega);
    momentsSet.calculateMPP(icell.TNW, omega);
    momentsSet.calculateMPM(icell.BNW, omega);
    momentsSet.calculatePMM(icell.BSE, omega);
    momentsSet.calculatePPP(icell.TNE, omega);
    momentsSet.calculatePMP(icell.TSE, omega);
    momentsSet.calculatePPM(icell.BNE, omega);

    momentsSet.calculateCoefficients(coefficients, xoff, yoff, zoff);
}

CompressibleOffsetMomentsInterpolator::CompressibleOffsetMomentsInterpolator(real omegaC, real omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{
}

InterpolationProcessorPtr CompressibleOffsetMomentsInterpolator::clone()
{
   return InterpolationProcessorPtr (new CompressibleOffsetMomentsInterpolator(this->omegaC, this->omegaF));
}

void CompressibleOffsetMomentsInterpolator::setOmegas(real omegaC, real omegaF)
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}

void CompressibleOffsetMomentsInterpolator::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff)
{
    vf::lbm::InterpolationCoefficients coefficients;
    calculateCoefficients(coefficients, icellC, omegaC, xoff, yoff, zoff);

     vf::lbm::interpolateCF(icellF.BSW, omegaF, vf::basics::constant::c1o2, coefficients, -0.25, -0.25, -0.25);
     vf::lbm::interpolateCF(icellF.BNE, omegaF, vf::basics::constant::c1o2, coefficients,  0.25,  0.25, -0.25);
     vf::lbm::interpolateCF(icellF.TNW, omegaF, vf::basics::constant::c1o2, coefficients, -0.25,  0.25,  0.25);
     vf::lbm::interpolateCF(icellF.TSE, omegaF, vf::basics::constant::c1o2, coefficients,  0.25, -0.25,  0.25);
     vf::lbm::interpolateCF(icellF.BNW, omegaF, vf::basics::constant::c1o2, coefficients, -0.25,  0.25, -0.25);
     vf::lbm::interpolateCF(icellF.BSE, omegaF, vf::basics::constant::c1o2, coefficients,  0.25, -0.25, -0.25);
     vf::lbm::interpolateCF(icellF.TSW, omegaF, vf::basics::constant::c1o2, coefficients, -0.25, -0.25,  0.25);
     vf::lbm::interpolateCF(icellF.TNE, omegaF, vf::basics::constant::c1o2, coefficients,  0.25,  0.25,  0.25);
}

void CompressibleOffsetMomentsInterpolator::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
    vf::lbm::InterpolationCoefficients coefficients;
    calculateCoefficients(coefficients, icellF, omegaF, xoff, yoff, zoff);

    vf::lbm::interpolateFC(icellC, vf::basics::constant::c2o1, omegaC, coefficients);
}
