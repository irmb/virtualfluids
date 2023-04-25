#ifndef VoidLBMKernel_h__
#define VoidLBMKernel_h__

#include "LBMKernel.h"

class VoidLBMKernel : public LBMKernel
{
public:
    VoidLBMKernel();
    ~VoidLBMKernel() override;
    SPtr<LBMKernel> clone() override;
    void calculate(int step) override;
    real getCalculationTime() override;
    void initDataSet();

protected:
};
#endif // VoidLBMKernel_h__
