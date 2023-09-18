#ifndef ActuatorFarm_H
#define ActuatorFarm_H

#include "PreCollisionInteractor.h"
#include "PointerDefinitions.h"
#include "basics/constants/NumericConstants.h"
#include <stdexcept>

using namespace vf::basics::constant;

class Parameter;
class GridProvider;
using namespace vf::basics::constant;

class ActuatorFarm : public PreCollisionInteractor
{
public:
    ActuatorFarm(
        const uint _nBlades,
        const real _density,
        const uint _nBladeNodes,
        const real _epsilon,
        int _level,
        const real _deltaT,
        const real _deltaX,
        const bool _useHostArrays
    ) :
        numberOfBlades(_nBlades),
        density(_density),
        numberOfBladeNodes(_nBladeNodes), 
        epsilon(_epsilon),
        level(_level),
        useHostArrays(_useHostArrays),
        numberOfTurbines(0),
        numberOfGridNodes(0),
        PreCollisionInteractor()
    {
        this->deltaT = _deltaT*exp2(-this->level);
        this->deltaX = _deltaX*exp2(-this->level);
        this->invEpsilonSqrd = 1/(epsilon*epsilon);
        this->invDeltaX = c1o1/this->deltaX;
     
        if(this->epsilon<this->deltaX)
            throw std::runtime_error("ActuatorFarm::ActuatorFarm: epsilon needs to be larger than dx!");
    }

    ~ActuatorFarm() override = default;
    void addTurbine(real turbinePosX, real turbinePosY, real turbinePosZ, real diameter, real omega, real azimuth, real yaw, std::vector<real> bladeRadii);
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;
    void getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider) override;

    void write(uint t);

    real getDensity(){ return this->density; };
    real getDeltaT(){ return this->deltaT; };
    real getDeltaX(){ return this->deltaX; };

    uint getNumberOfTurbines(){ return this->numberOfTurbines; };
    uint getNumberOfNodesPerBlade(){ return this->numberOfBladeNodes; };
    uint getNumberOfBladesPerTurbine(){ return this->numberOfBlades; };

    uint getNumberOfIndices(){ return this->numberOfIndices; };
    uint getNumberOfGridNodes(){ return this->numberOfGridNodes; };

    real* getAllAzimuths(){ return azimuthsH; };
    real* getAllOmegas(){ return omegasH; };
    real* getAllYaws(){ return yawsH; };

    real* getAllTurbinePosX(){ return turbinePosXH; };
    real* getAllTurbinePosY(){ return turbinePosYH; };
    real* getAllTurbinePosZ(){ return turbinePosZH; };

    real getTurbineAzimuth(uint turbine){ return azimuthsH[turbine]; };
    real getTurbineOmega  (uint turbine){ return omegasH[turbine];   };
    real getTurbineYaw    (uint turbine){ return yawsH[turbine];     };

    real getTurbinePosX(uint turbine){ return turbinePosXH[turbine]; };
    real getTurbinePosY(uint turbine){ return turbinePosYH[turbine]; };
    real getTurbinePosZ(uint turbine){ return turbinePosZH[turbine]; };

    real* getAllBladeRadii(){ return this->bladeRadiiH; };
    real* getAllBladeCoordsX(){ return this->bladeCoordsXH; };
    real* getAllBladeCoordsY(){ return this->bladeCoordsYH; };
    real* getAllBladeCoordsZ(){ return this->bladeCoordsZH; };
    real* getAllBladeVelocitiesX(){ return this->bladeVelocitiesXH; };
    real* getAllBladeVelocitiesY(){ return this->bladeVelocitiesYH; };
    real* getAllBladeVelocitiesZ(){ return this->bladeVelocitiesZH; };
    real* getAllBladeForcesX(){ return this->bladeForcesXH; };
    real* getAllBladeForcesY(){ return this->bladeForcesYH; };
    real* getAllBladeForcesZ(){ return this->bladeForcesZH; };

    real* getTurbineBladeRadii(uint turbine){ return &this->bladeRadiiH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeCoordsX(uint turbine){ return &this->bladeCoordsXH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeCoordsY(uint turbine){ return &this->bladeCoordsYH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeCoordsZ(uint turbine){ return &this->bladeCoordsZH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesX(uint turbine){ return &this->bladeVelocitiesXH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesY(uint turbine){ return &this->bladeVelocitiesYH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesZ(uint turbine){ return &this->bladeVelocitiesZH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesX(uint turbine){ return &this->bladeForcesXH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesY(uint turbine){ return &this->bladeForcesYH[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesZ(uint turbine){ return &this->bladeForcesZH[turbine*numberOfBladeNodes*numberOfBlades]; };

    real* getAllBladeRadiiDevice(){ return this->bladeRadiiD; };
    real* getAllBladeCoordsXDevice(){ return this->bladeCoordsXDCurrentTimestep; };
    real* getAllBladeCoordsYDevice(){ return this->bladeCoordsYDCurrentTimestep; };
    real* getAllBladeCoordsZDevice(){ return this->bladeCoordsZDCurrentTimestep; };
    real* getAllBladeVelocitiesXDevice(){ return this->bladeVelocitiesXDCurrentTimestep; };
    real* getAllBladeVelocitiesYDevice(){ return this->bladeVelocitiesYDCurrentTimestep; };
    real* getAllBladeVelocitiesZDevice(){ return this->bladeVelocitiesZDCurrentTimestep; };
    real* getAllBladeForcesXDevice(){ return this->bladeForcesXDCurrentTimestep; };
    real* getAllBladeForcesYDevice(){ return this->bladeForcesYDCurrentTimestep; };
    real* getAllBladeForcesZDevice(){ return this->bladeForcesZDCurrentTimestep; };

    real* getTurbineBladeRadiiDevice(uint turbine){ return &this->bladeRadiiD[turbine*numberOfBladeNodes]; };
    real* getTurbineBladeCoordsXDevice(uint turbine){ return &this->bladeCoordsXDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeCoordsYDevice(uint turbine){ return &this->bladeCoordsYDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeCoordsZDevice(uint turbine){ return &this->bladeCoordsZDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesXDevice(uint turbine){ return &this->bladeVelocitiesXDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesYDevice(uint turbine){ return &this->bladeVelocitiesYDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeVelocitiesZDevice(uint turbine){ return &this->bladeVelocitiesZDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesXDevice(uint turbine){ return &this->bladeForcesXDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesYDevice(uint turbine){ return &this->bladeForcesYDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };
    real* getTurbineBladeForcesZDevice(uint turbine){ return &this->bladeForcesZDCurrentTimestep[turbine*numberOfBladeNodes*numberOfBlades]; };

    void setAllAzimuths(real* _azimuth);
    void setAllOmegas(real* _omegas);
    void setAllYaws(real* yaws);
    
    void setTurbineAzimuth(uint turbine, real azimuth){ azimuthsH[turbine] = azimuth; };
    void setTurbineYaw(uint turbine, real yaw){ yawsH[turbine] = yaw; };
    void setTurbineOmega(uint turbine, real omega){ omegasH[turbine] = omega; };

    void setAllBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ);
    void setAllBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ);
    void setAllBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ);

    void setTurbineBladeCoords(uint turbine, real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ);
    void setTurbineBladeVelocities(uint turbine, real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ);
    void setTurbineBladeForces(uint turbine, real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ);

    virtual void calcBladeForces();

private:
    void initTurbineGeometries(CudaMemoryManager* cudaManager);
    void initBoundingSpheres(Parameter* para, CudaMemoryManager* cudaManager);
    void initBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeForces(CudaMemoryManager* cudaManager);
    void initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager);

    void calcForcesEllipticWing();
    void rotateBlades(real angle, uint turbineID);

    void writeBladeCoords(uint t);
    void writeBladeCoordsToVtkForDebug(const std::string& filename);
    void writeBladeForces(uint t);
    void writeBladeVelocities(uint t);

    void swapDeviceArrays();

public:
    real* bladeRadiiH;
    real* bladeRadiiD;
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
    real* turbinePosXH, *turbinePosYH, *turbinePosZH, *omegasH, *azimuthsH, *yawsH, *diametersH;
    real* turbinePosXD, *turbinePosYD, *turbinePosZD, *omegasD, *azimuthsD, *yawsD, *diametersD;
    
private:
    std::vector<real> preInitPosX, preInitPosY, preInitPosZ, preInitDiameters, preInitOmegas, preInitAzimuths, preInitYaws;
    std::vector<std::vector<real>> preInitBladeRadii;
    const bool useHostArrays;
    const real density;
    real deltaT, deltaX;
    const uint numberOfBladeNodes, numberOfBlades;
    uint numberOfTurbines;
    const real epsilon; // in m
    const int level;
    uint numberOfIndices;
    uint numberOfGridNodes;
    real forceRatio, factorGaussian, invEpsilonSqrd, invDeltaX;
    int streamIndex;
};

#endif
