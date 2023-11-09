#ifndef SlipBCStrategy_h__
#define SlipBCStrategy_h__

#include "BCStrategy.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class SlipBCStrategy : public BCStrategy
{
public:
    SlipBCStrategy();
    ~SlipBCStrategy() override;
    SPtr<BCStrategy> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void applyBC() override;
};
#endif // SlipBCStrategy_h__
