#ifndef ThinWallNoSlipBCAlgorithm_h__
#define ThinWallNoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"
#include <PointerDefinitions.h>

class DistributionArray3D;

class ThinWallNoSlipBCAlgorithm : public BCAlgorithm
{
public:
    ThinWallNoSlipBCAlgorithm();
    ~ThinWallNoSlipBCAlgorithm() override;
    SPtr<BCAlgorithm> clone() override;
    void addDistributions(SPtr<DistributionArray3D> distributions) override;
    void setPass(int pass);
    void applyBC() override;

protected:
    SPtr<DistributionArray3D> distributionsTemp;

private:
    int pass;
    LBMReal fTemp[D3Q27System::ENDF + 1];
};
#endif // ThinWallNoSlipBCAlgorithm_h__
