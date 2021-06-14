#ifndef BGKLBMKernel_H
#define BGKLBMKernel_H

#include "LBMKernel.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

class BGKLBMKernel : public LBMKernel
{
public:
    BGKLBMKernel();
    ~BGKLBMKernel() override;
    void calculate(int step) override;
    SPtr<LBMKernel> clone() override;
    double getCalculationTime() override;

private:
    void initDataSet();
    // void collideAllCompressible();
    // void collideAllIncompressible();

    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;

    mu::value_type muX1, muX2, muX3;
    LBMReal forcingX1;
    LBMReal forcingX2;
    LBMReal forcingX3;
};

#endif
