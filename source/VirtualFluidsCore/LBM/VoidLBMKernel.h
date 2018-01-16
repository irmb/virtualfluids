#ifndef VoidLBMKernel_h__
#define VoidLBMKernel_h__

#include "LBMKernel.h"

class VoidLBMKernel : public LBMKernel
{
public:
   VoidLBMKernel();
   ~VoidLBMKernel();
   SPtr<LBMKernel> clone();
   void calculate();
   double getCalculationTime();
   void initDataSet();
protected:

};
#endif // VoidLBMKernel_h__
