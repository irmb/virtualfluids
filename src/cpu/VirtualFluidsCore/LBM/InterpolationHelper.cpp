#include "InterpolationHelper.h"

InterpolationHelper::InterpolationHelper(InterpolationProcessorPtr iProcessor) : iProcessor(iProcessor) {}
//////////////////////////////////////////////////////////////////////////
InterpolationHelper::~InterpolationHelper() = default;
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1(D3Q27ICell &icellF, real *icellC, double /*x1*/, double /*x2*/,
                                          double /*x3*/, real omega)
{
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, 1.0);
    iProcessor->calcInterpolatedNodeFC(icellC, omega);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1WithVelocity(D3Q27ICell &icellF, double x1, double x2, double x3,
                                                      real omega, real &vx1, real &vx2, real &vx3)
{
    iProcessor->setOffsets(0.0, 0.0, 0.0);
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, 0.0);
    iProcessor->calcInterpolatedVelocity(x1, x2, x3, vx1, vx2, vx3);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationHelper::interpolate8to1WithVelocityWithShearStress(D3Q27ICell &icellF, double x1, double x2,
                                                                     double x3, real omega, real &vx1,
                                                                     real &vx2, real &vx3, real &tauxx,
                                                                     real &tauyy, real &tauzz, real &tauxy,
                                                                     real &tauxz, real &tauyz)
{
    iProcessor->setOffsets(0.0, 0.0, 0.0);
    iProcessor->calcInterpolatedCoefficiets(icellF, omega, 0.0);
    iProcessor->calcInterpolatedVelocity(x1, x2, x3, vx1, vx2, vx3);
    iProcessor->calcInterpolatedShearStress(x1, x2, x3, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz);
}

//////////////////////////////////////////////////////////////////////////
