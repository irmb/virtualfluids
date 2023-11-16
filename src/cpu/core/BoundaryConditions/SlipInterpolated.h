#ifndef SlipInterpolated_h__
#define SlipInterpolated_h__

#include "BCStrategy.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class SlipInterpolated : public BCStrategy
{
public:
    SlipInterpolated();
    ~SlipInterpolated() override;
    SPtr<BCStrategy> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void applyBC() override;
};
#endif // SlipBCStrategy_h__
