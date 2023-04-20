#include "InterpolationHelper.h"

InterpolationHelper::InterpolationHelper(InterpolationProcessorPtr iProcessor) : iProcessor(iProcessor) {}
//////////////////////////////////////////////////////////////////////////
InterpolationHelper::~InterpolationHelper() = default;
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1(D3Q27ICell &icellF, real *icellC, real /*x1*/, real /*x2*/,
                                          real /*x3*/, real omega)
{
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, vf::basics::constant::c1o1);
    iProcessor->calcInterpolatedNodeFC(icellC, omega);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1WithVelocity(D3Q27ICell &icellF, real x1, real x2, real x3,
                                                      real omega, real &vx1, real &vx2, real &vx3)
{
    iProcessor->setOffsets(vf::basics::constant::c0o1, vf::basics::constant::c0o1, vf::basics::constant::c0o1);
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, vf::basics::constant::c0o1);
    iProcessor->calcInterpolatedVelocity(x1, x2, x3, vx1, vx2, vx3);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1WithVelocityWithShearStress(D3Q27ICell &icellF, real x1, real x2,
                                                                     real x3, real omega, real &vx1,
                                                                     real &vx2, real &vx3, real &tauxx,
                                                                     real &tauyy, real &tauzz, real &tauxy,
                                                                     real &tauxz, real &tauyz)
{
    iProcessor->setOffsets(vf::basics::constant::c0o1, vf::basics::constant::c0o1, vf::basics::constant::c0o1);
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, vf::basics::constant::c0o1);
    iProcessor->calcInterpolatedVelocity(x1, x2, x3, vx1, vx2, vx3);
    iProcessor->calcInterpolatedShearStress(x1, x2, x3, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz);
}

//////////////////////////////////////////////////////////////////////////
