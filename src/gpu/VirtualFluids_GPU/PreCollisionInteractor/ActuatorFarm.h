#ifndef ActuatorFarm_H
#define ActuatorFarm_H

#include "PreCollisionInteractor.h"
#include "PointerDefinitions.h"

class Parameter;
class GridProvider;

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
        nBlades(_nBlades),
        density(_density),
        nBladeNodes(_nBladeNodes), 
        epsilon(_epsilon),
        level(_level),
        useHostArrays(_useHostArrays),
        nTurbines(0),
        numberOfNodes(0),
        PreCollisionInteractor()
    {
        this->deltaT = _deltaX/pow(2,this->level);
        this->deltaX = _deltaX/pow(2,this->level);
    }

    virtual  ~ActuatorFarm()
    {
        
    }
    void addTurbine(real turbinePosX, real turbinePosY, real turbinePosZ, real diameter, real omega, real azimuth, real yaw, std::vector<real> bladeRadii);
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);
    void write(uint t);

    uint getNTurbines(){return this->nTurbines; };
    uint getNBladeNodes(){return this->nBladeNodes;};
    uint getNBlades(){return this->nBlades;};

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

    uint getNumberOfIndices(){return this->numberOfIndices;};
    uint getNumberOfNodes(){return this->numberOfNodes;};
    real* getAllBladeCoordsX(){return this->bladeCoordsXH;};
    real* getAllBladeCoordsY(){return this->bladeCoordsYH;};
    real* getAllBladeCoordsZ(){return this->bladeCoordsZH;};
    real* getAllBladeVelocitiesX(){return this->bladeVelocitiesXH;};
    real* getAllBladeVelocitiesY(){return this->bladeVelocitiesYH;};
    real* getAllBladeVelocitiesZ(){return this->bladeVelocitiesZH;};
    real* getAllBladeForcesX(){return this->bladeForcesXH;};
    real* getAllBladeForcesY(){return this->bladeForcesYH;};
    real* getAllBladeForcesZ(){return this->bladeForcesZH;};

    real* getTurbineBladeCoordsX(uint turbine){return &this->bladeCoordsXH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeCoordsY(uint turbine){return &this->bladeCoordsYH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeCoordsZ(uint turbine){return &this->bladeCoordsZH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesX(uint turbine){return &this->bladeVelocitiesXH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesY(uint turbine){return &this->bladeVelocitiesYH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesZ(uint turbine){return &this->bladeVelocitiesZH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesX(uint turbine){return &this->bladeForcesXH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesY(uint turbine){return &this->bladeForcesYH[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesZ(uint turbine){return &this->bladeForcesZH[turbine*nBladeNodes*nBlades];};

    real* getAllBladeCoordsXDevice(){return this->bladeCoordsXD;};
    real* getAllBladeCoordsYDevice(){return this->bladeCoordsYD;};
    real* getAllBladeCoordsZDevice(){return this->bladeCoordsZD;};
    real* getAllBladeVelocitiesXDevice(){return this->bladeVelocitiesXD;};
    real* getAllBladeVelocitiesYDevice(){return this->bladeVelocitiesYD;};
    real* getAllBladeVelocitiesZDevice(){return this->bladeVelocitiesZD;};
    real* getAllBladeForcesXDevice(){return this->bladeForcesXD;};
    real* getAllBladeForcesYDevice(){return this->bladeForcesYD;};
    real* getAllBladeForcesZDevice(){return this->bladeForcesZD;};

    real* getTurbineBladeCoordsXDevice(uint turbine){return &this->bladeCoordsXD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeCoordsYDevice(uint turbine){return &this->bladeCoordsYD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeCoordsZDevice(uint turbine){return &this->bladeCoordsZD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesXDevice(uint turbine){return &this->bladeVelocitiesXD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesYDevice(uint turbine){return &this->bladeVelocitiesYD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeVelocitiesZDevice(uint turbine){return &this->bladeVelocitiesZD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesXDevice(uint turbine){return &this->bladeForcesXD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesYDevice(uint turbine){return &this->bladeForcesYD[turbine*nBladeNodes*nBlades];};
    real* getTurbineBladeForcesZDevice(uint turbine){return &this->bladeForcesZD[turbine*nBladeNodes*nBlades];};

    void setAllBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ);
    void setAllBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ);
    void setAllBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ);

    void setTurbineBladeCoords(uint turbine, real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ);
    void setTurbineBladeVelocities(uint turbine, real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ);
    void setTurbineBladeForces(uint turbine, real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ);


private:
    void initTurbineGeometries(CudaMemoryManager* cudaManager);
    void initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager);
    void initBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeForces(CudaMemoryManager* cudaManager);
    void initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager);

    void calcForcesEllipticWing();
    void calcBladeForces();
    void rotateBlades(real angle, uint turbineID);

    void writeBladeCoords(uint t);
    void writeBladeForces(uint t);
    void writeBladeVelocities(uint t);


    

public:
    real* bladeRadiiH;
    real* bladeRadiiD;
    real* bladeCoordsXH, * bladeCoordsYH, * bladeCoordsZH;
    real* bladeCoordsXD, * bladeCoordsYD, * bladeCoordsZD;
    real* bladeVelocitiesXH, * bladeVelocitiesYH, * bladeVelocitiesZH;
    real* bladeVelocitiesXD, * bladeVelocitiesYD, * bladeVelocitiesZD;
    real* bladeForcesXH, * bladeForcesYH, * bladeForcesZH;
    real* bladeForcesXD, * bladeForcesYD, * bladeForcesZD;
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
    const uint nBladeNodes, nBlades;
    uint nTurbines;
    const real epsilon; // in m
    const int level;
    uint numberOfIndices;
    uint numberOfNodes;
    real forceRatio, factorGaussian, invEpsilonSqrd, invDeltaX;
};

#endif