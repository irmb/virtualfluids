#ifndef ActuatorFarmStandalone_H
#define ActuatorFarmStandalone_H

#include "ActuatorFarm.h"
#include "basics/DataTypes.h"

class ActuatorFarmStandalone : public ActuatorFarm
{
public:
    ActuatorFarmStandalone(
        const real diameter,
        const uint numberOfNodesPerBlade,
        const std::vector<real> turbinePositionsX,
        const std::vector<real> turbinePositionsY,
        const std::vector<real> turbinePositionsZ,
        const std::vector<real> rotorSpeeds,
        const real density,
        const real smearingWidth,
        const int level,
        const real deltaT,
        const real deltaX
    ) : rotorSpeeds(rotorSpeeds),
        ActuatorFarm(diameter, computeBladeRadii(diameter, numberOfNodesPerBlade), turbinePositionsX, turbinePositionsY, turbinePositionsZ, density, smearingWidth, level, deltaT, deltaX, true)
    {}

    ~ActuatorFarmStandalone() = default;

    void updateForcesAndCoordinates() override;
    static std::vector<real> computeBladeRadii(const real diameter, const uint numberOfNodesPerBlade);
    
private:
    std::vector<real> rotorSpeeds;

};

#endif
