#ifndef HighViscosityNoSlipBCAlgorithm_h__
#define HighViscosityNoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class HighViscosityNoSlipBCAlgorithm : public BCAlgorithm
{
public:
    HighViscosityNoSlipBCAlgorithm();
    ~HighViscosityNoSlipBCAlgorithm() override;
    SPtr<BCAlgorithm> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void applyBC() override;
};
#endif // HighViscosityNoSlipBCAlgorithm_h__
