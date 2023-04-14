#ifndef D3Q27InterpolationHelper_H_
#define D3Q27InterpolationHelper_H_

#include "InterpolationProcessor.h"

class InterpolationHelper;
using InterpolationHelperPtr = SPtr<InterpolationHelper>;

class InterpolationHelper
{
public:
    InterpolationHelper(InterpolationProcessorPtr iProcessor);
    ~InterpolationHelper();
    void interpolate8to1(D3Q27ICell &icellF, real *icellC, real x1, real x2, real x3, real omega);
    void interpolate8to1WithVelocity(D3Q27ICell &icellF, real x1, real x2, real x3, real omega, real &vx1,
                                     real &vx2, real &vx3);
    void interpolate8to1WithVelocityWithShearStress(D3Q27ICell &icellF, real x1, real x2, real x3, real omega,
                                                    real &vx1, real &vx2, real &vx3, real &tauxx,
                                                    real &tauyy, real &tauzz, real &tauxy, real &tauxz,
                                                    real &tauyz);

protected:
private:
    InterpolationProcessorPtr iProcessor;
};

#endif
