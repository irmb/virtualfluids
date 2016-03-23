#ifndef D3Q27InterpolationHelper_H_
#define D3Q27InterpolationHelper_H_

#include "D3Q27InterpolationProcessor.h"

class D3Q27InterpolationHelper;
typedef boost::shared_ptr<D3Q27InterpolationHelper> D3Q27InterpolationHelperPtr;

class D3Q27InterpolationHelper
{
public:
   D3Q27InterpolationHelper(D3Q27InterpolationProcessorPtr iProcessor);
   ~D3Q27InterpolationHelper();
   void interpolate8to1(D3Q27ICell& icellF, LBMReal* icellC, double x1, double x2, double x3, LBMReal omega);
   void interpolate8to1WithVelocity(D3Q27ICell& icellF, double x1, double x2, double x3, LBMReal omega, LBMReal &vx1, LBMReal &vx2, LBMReal &vx3);
   void interpolate8to1WithVelocityWithShearStress(D3Q27ICell& icellF, double x1, double x2, double x3, LBMReal omega, 
                                             LBMReal &vx1, LBMReal &vx2, LBMReal &vx3, 
                                             LBMReal &tauxx, LBMReal &tauyy, LBMReal &tauzz,LBMReal &tauxy, LBMReal &tauxz, LBMReal &tauyz);
protected:
private:
   D3Q27InterpolationProcessorPtr iProcessor;
};

#endif

