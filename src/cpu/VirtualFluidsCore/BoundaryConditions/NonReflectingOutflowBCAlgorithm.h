#ifndef NonReflectingOutflowBCAlgorithm_h__
#define NonReflectingOutflowBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class NonReflectingOutflowBCAlgorithm : public BCAlgorithm
{
public:
    NonReflectingOutflowBCAlgorithm();
    ~NonReflectingOutflowBCAlgorithm() override;
    SPtr<BCAlgorithm> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void applyBC() override;
};
#endif // NonReflectingDensityBCAlgorithm_h__
