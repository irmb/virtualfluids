#ifndef InitDensityLBMKernel_h__
#define InitDensityLBMKernel_h__

#include "LBMKernelETD3Q27.h"
#include "basics/utilities/UbTiming.h"

class InitDensityLBMKernel :  public LBMKernelETD3Q27
{
public:
   InitDensityLBMKernel();
   ~InitDensityLBMKernel();
   InitDensityLBMKernel(int nx1, int nx2, int nx3);
   void calculate();
   LBMKernel3DPtr clone();
   void setVelocity(int x1, int x2, int x3, LBMReal vvx, LBMReal vvy, LBMReal vvz);
   double getCallculationTime();
protected:
   void init();
   void collideAll();
private:
   LBMReal f[D3Q27System::ENDF+1];
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;
   LBMReal OxyyMxzz;
   CbArray4D<LBMReal, IndexerX4X3X2X1> v;
};

#endif // InitDensityLBMKernel_h__

