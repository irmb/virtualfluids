#ifndef ActuatorFarm_H
#define ActuatorFarm_H

#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <stdexcept>

class Parameter;
class GridProvider;
using namespace vf::basics::constant;

class ActuatorFarm : public PreCollisionInteractor
{
public:
    ActuatorFarm(
        const real diameter,
        const std::vector<real> bladeRadii,
        const std::vector<real> turbinePositionsX,
        const std::vector<real> turbinePositionsY,
        const std::vector<real> turbinePositionsZ,
        const real density,
        const real smearingWidth,
        const int level,
        const real deltaT,
        const real deltaX,
        const bool useHostArrays
    ) :
        diameter(diameter),
        bladeRadii(bladeRadii),
        numberOfNodesPerBlade(static_cast<uint>(bladeRadii.size())),
        numberOfNodesPerTurbine(numberOfNodesPerBlade*numberOfBlades),
        numberOfTurbines(static_cast<uint>(turbinePositionsX.size())),
        initialTurbinePositionsX(turbinePositionsX),
        initialTurbinePositionsY(turbinePositionsY),
        initialTurbinePositionsZ(turbinePositionsZ),
        density(density),
        smearingWidth(smearingWidth),
        level(level),
        useHostArrays(useHostArrays),
        deltaT(deltaT*exp2(-level)),
        deltaX(deltaX*exp2(-level)),
        invDeltaX(c1o1/deltaX)
    {
        if(this->smearingWidth < this->deltaX)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: smearing width needs to be larger than dx!");
        if(numberOfTurbines != turbinePositionsY.size() || numberOfTurbines != turbinePositionsZ.size())
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: turbine positions need to have the same length!");
        azimuths = std::vector<real>(numberOfTurbines, 0.0);
    }

    ~ActuatorFarm() override = default;
    void init(Parameter* para, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;
    void getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider) override;

    void enableOutput(const std::string outputName, uint tStart, uint tOut) {
        this->outputName = outputName;
        this->writeOutput = true;
        this->tStartOut = tStart;
        this->tOut = tOut;
    }

    void write(const std::string& filename) const;

    real getDensity() const { return this->density; };
    real getDeltaT() const { return this->deltaT; };
    real getDeltaX() const { return this->deltaX; };

    uint getNumberOfTurbines() const { return this->numberOfTurbines; };
    uint getNumberOfNodesPerTurbine() const { return this->numberOfNodesPerTurbine; };
    uint getNumberOfNodesPerBlade() const { return this->numberOfNodesPerBlade; };
    uint getNumberOfBladesPerTurbine() const { return ActuatorFarm::numberOfBlades; };

    uint getNumberOfIndices() const { return this->numberOfIndices; };
    uint getNumberOfGridNodes() const { return this->numberOfGridNodes; };

    real* getAllTurbinePosX() const { return turbinePosXH; };
    real* getAllTurbinePosY() const { return turbinePosYH; };
    real* getAllTurbinePosZ() const { return turbinePosZH; };

    real getTurbinePosX(uint turbine) const { return turbinePosXH[turbine]; };
    real getTurbinePosY(uint turbine) const { return turbinePosYH[turbine]; };
    real getTurbinePosZ(uint turbine) const { return turbinePosZH[turbine]; };

    real* getAllBladeCoordsX() const { return this->bladeCoordsXH; };
    real* getAllBladeCoordsY() const { return this->bladeCoordsYH; };
    real* getAllBladeCoordsZ() const { return this->bladeCoordsZH; };
    real* getAllBladeVelocitiesX() const { return this->bladeVelocitiesXH; };
    real* getAllBladeVelocitiesY() const { return this->bladeVelocitiesYH; };
    real* getAllBladeVelocitiesZ() const { return this->bladeVelocitiesZH; };
    real* getAllBladeForcesX() const { return this->bladeForcesXH; };
    real* getAllBladeForcesY() const { return this->bladeForcesYH; };
    real* getAllBladeForcesZ() const { return this->bladeForcesZH; };

    real* getTurbineBladeCoordsX(uint turbine) const { return &this->bladeCoordsXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsY(uint turbine) const { return &this->bladeCoordsYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsZ(uint turbine) const { return &this->bladeCoordsZH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesX(uint turbine) const { return &this->bladeVelocitiesXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesY(uint turbine) const { return &this->bladeVelocitiesYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesZ(uint turbine) const { return &this->bladeVelocitiesZH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesX(uint turbine) const { return &this->bladeForcesXH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesY(uint turbine) const { return &this->bladeForcesYH[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesZ(uint turbine) const { return &this->bladeForcesZH[turbine*numberOfNodesPerTurbine]; };

    real* getAllBladeCoordsXDevice() const { return this->bladeCoordsXDCurrentTimestep; };
    real* getAllBladeCoordsYDevice() const { return this->bladeCoordsYDCurrentTimestep; };
    real* getAllBladeCoordsZDevice() const { return this->bladeCoordsZDCurrentTimestep; };
    real* getAllBladeVelocitiesXDevice() const { return this->bladeVelocitiesXDCurrentTimestep; };
    real* getAllBladeVelocitiesYDevice() const { return this->bladeVelocitiesYDCurrentTimestep; };
    real* getAllBladeVelocitiesZDevice() const { return this->bladeVelocitiesZDCurrentTimestep; };
    real* getAllBladeForcesXDevice() const { return this->bladeForcesXDCurrentTimestep; };
    real* getAllBladeForcesYDevice() const { return this->bladeForcesYDCurrentTimestep; };
    real* getAllBladeForcesZDevice() const { return this->bladeForcesZDCurrentTimestep; };

    real* getTurbineBladeCoordsXDevice(uint turbine) const { return &this->bladeCoordsXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsYDevice(uint turbine) const { return &this->bladeCoordsYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeCoordsZDevice(uint turbine) const { return &this->bladeCoordsZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesXDevice(uint turbine) const { return &this->bladeVelocitiesXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesYDevice(uint turbine) const { return &this->bladeVelocitiesYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeVelocitiesZDevice(uint turbine) const { return &this->bladeVelocitiesZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesXDevice(uint turbine) const { return &this->bladeForcesXDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesYDevice(uint turbine) const { return &this->bladeForcesYDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };
    real* getTurbineBladeForcesZDevice(uint turbine) const { return &this->bladeForcesZDCurrentTimestep[turbine*numberOfNodesPerTurbine]; };

    void setAllBladeCoords(const real* _bladeCoordsX, const real* _bladeCoordsY, const real* _bladeCoordsZ);
    void setAllBladeVelocities(const real* _bladeVelocitiesX, const real* _bladeVelocitiesY, const real* _bladeVelocitiesZ);
    void setAllBladeForces(const real* _bladeForcesX, const real* _bladeForcesY, const real* _bladeForcesZ);

    void setTurbineBladeCoords(uint turbine, const real* _bladeCoordsX, const real* _bladeCoordsY, const real* _bladeCoordsZ);
    void setTurbineBladeVelocities(uint turbine, const real* _bladeVelocitiesX, const real* _bladeVelocitiesY, const real* _bladeVelocitiesZ);
    void setTurbineBladeForces(uint turbine, const real* _bladeForcesX, const real* _bladeForcesY, const real* _bladeForcesZ);

    void setTurbineAzimuth(uint turbine, real azimuth){azimuths[turbine] = azimuth;}

    virtual void updateForcesAndCoordinates()=0;

private:
    void initTurbineGeometries(CudaMemoryManager* cudaManager);
    void initBoundingSpheres(Parameter* para, CudaMemoryManager* cudaManager);
    void initBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeForces(CudaMemoryManager* cudaManager);
    void initBladeIndices(CudaMemoryManager* cudaManager);
    std::string getFilename(Parameter* para, uint t) const;
    void swapDeviceArrays();

public:
    real* bladeCoordsXH, * bladeCoordsYH, * bladeCoordsZH;
    real* bladeCoordsXDPreviousTimestep, * bladeCoordsYDPreviousTimestep, * bladeCoordsZDPreviousTimestep;
    real* bladeCoordsXDCurrentTimestep, * bladeCoordsYDCurrentTimestep, * bladeCoordsZDCurrentTimestep;    
    real* bladeVelocitiesXH, * bladeVelocitiesYH, * bladeVelocitiesZH;
    real* bladeVelocitiesXDPreviousTimestep, * bladeVelocitiesYDPreviousTimestep, * bladeVelocitiesZDPreviousTimestep;
    real* bladeVelocitiesXDCurrentTimestep, * bladeVelocitiesYDCurrentTimestep, * bladeVelocitiesZDCurrentTimestep;
    real* bladeForcesXH, * bladeForcesYH, * bladeForcesZH;
    real* bladeForcesXDPreviousTimestep, * bladeForcesYDPreviousTimestep, * bladeForcesZDPreviousTimestep;
    real* bladeForcesXDCurrentTimestep, * bladeForcesYDCurrentTimestep, * bladeForcesZDCurrentTimestep;
    uint* bladeIndicesH;
    uint* bladeIndicesD; 
    uint* boundingSphereIndicesH;
    uint* boundingSphereIndicesD;
    real* turbinePosXH, *turbinePosYH, *turbinePosZH;
    real* turbinePosXD, *turbinePosYD, *turbinePosZD;

protected:
    static constexpr uint numberOfBlades{3};
    std::vector<real> bladeRadii, initialTurbinePositionsX, initialTurbinePositionsY, initialTurbinePositionsZ;
    std::vector<real> azimuths;
    const real diameter;
    const bool useHostArrays;
    const real density;
    const real deltaT, deltaX, invDeltaX;
    const uint numberOfTurbines, numberOfNodesPerBlade, numberOfNodesPerTurbine;
    const real smearingWidth; // in m
    const int level;
    uint numberOfIndices{0};
    uint numberOfGridNodes{0};

    real forceRatio, factorGaussian;
    int streamIndex;

    bool writeOutput{false};
    std::string outputName;
    uint tOut{0};
    uint tStartOut{0};
};

#endif
