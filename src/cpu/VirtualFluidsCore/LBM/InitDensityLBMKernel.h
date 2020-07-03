#ifndef InitDensityLBMKernel_h__
#define InitDensityLBMKernel_h__

#include "LBMKernel.h"
#include "basics/utilities/UbTiming.h"
#include "CbArray4D.h"
#include "D3Q27System.h"
#include "CbArray3D.h"

class InitDensityLBMKernel :  public LBMKernel
{
public:
   InitDensityLBMKernel();
   ~InitDensityLBMKernel();
   void calculate(int step);
   SPtr<LBMKernel> clone();
   void setVelocity(int x1, int x2, int x3, LBMReal vvx, LBMReal vvy, LBMReal vvz);
   double getCalculationTime();
protected:
   void initDataSet();
private:
   LBMReal f[D3Q27System::ENDF+1];
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;
   LBMReal OxyyMxzz;
   CbArray4D<LBMReal, IndexerX4X3X2X1> v;
};

#endif // InitDensityLBMKernel_h__

