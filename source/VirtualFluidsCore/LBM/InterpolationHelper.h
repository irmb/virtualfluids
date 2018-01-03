#ifndef D3Q27InterpolationHelper_H_
#define D3Q27InterpolationHelper_H_

#include "InterpolationProcessor.h"

class InterpolationHelper;
typedef std::shared_ptr<InterpolationHelper> InterpolationHelperPtr;

class InterpolationHelper
{
public:
   InterpolationHelper(InterpolationProcessorPtr iProcessor);
   ~InterpolationHelper();
   void interpolate8to1(D3Q27ICell& icellF, LBMReal* icellC, double x1, double x2, double x3, LBMReal omega);
   void interpolate8to1WithVelocity(D3Q27ICell& icellF, double x1, double x2, double x3, LBMReal omega, LBMReal &vx1, LBMReal &vx2, LBMReal &vx3);
   void interpolate8to1WithVelocityWithShearStress(D3Q27ICell& icellF, double x1, double x2, double x3, LBMReal omega, 
                                             LBMReal &vx1, LBMReal &vx2, LBMReal &vx3, 
                                             LBMReal &tauxx, LBMReal &tauyy, LBMReal &tauzz,LBMReal &tauxy, LBMReal &tauxz, LBMReal &tauyz);
protected:
private:
   InterpolationProcessorPtr iProcessor;
};

#endif

