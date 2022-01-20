#ifndef ActuatorLine_H
#define ActuatorLine_H

#include "PreCollisionInteractor.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

class Parameter;
class GridProvider;

class VIRTUALFLUIDS_GPU_EXPORT ActuatorLine : public PreCollisionInteractor
{
public:
    ActuatorLine(
        const uint _nBlades,
        const real _density,
        const uint _nBladeNodes,
        const real _epsilon,
        real _turbinePosX, real _turbinePosY, real _turbinePosZ,
        const real _diameter,
        int _level,
        const real _delta_t,
        const real _delta_x
    ) : nBlades(_nBlades),
        density(_density),
        nBladeNodes(_nBladeNodes), 
        epsilon(_epsilon),
        turbinePosX(_turbinePosX), turbinePosY(_turbinePosY), turbinePosZ(_turbinePosZ),
        diameter(_diameter),
        level(_level),
        delta_x(_delta_x),
        PreCollisionInteractor()
    {
        this->delta_t = _delta_t/pow(2,this->level);
        this->delta_x = _delta_x/pow(2,this->level);
        this->numberOfNodes = this->nBladeNodes*this->nBlades;
        this->omega = 1.0f;
        this->azimuth = 0.0f;
        this->yaw = 0.0f;
    };

    virtual ~ActuatorLine(){};

    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;
    void write(uint t);

    uint getNBladeNodes(){return this->nBladeNodes;};
    uint getNBlades(){return this->nBlades;};
    uint getNumberOfIndices(){return this->numberOfIndices;};
    uint getNumberOfNodes(){return this->numberOfNodes;};
    real getOmega(){ return this->omega; };
    real getAzimuth(){ return this->azimuth; };
    real getDensity(){ return this->density; };
    real getPositionX(){ return this->turbinePosX; };
    real getPositionY(){ return this->turbinePosY; };
    real getPositionZ(){ return this->turbinePosZ; };
    real* getBladeRadii(){return this->bladeRadiiH;};
    real* getBladeCoordsX(){return this->bladeCoordsXH;};
    real* getBladeCoordsY(){return this->bladeCoordsYH;};
    real* getBladeCoordsZ(){return this->bladeCoordsZH;};
    real* getBladeVelocitiesX(){return this->bladeVelocitiesXH;};
    real* getBladeVelocitiesY(){return this->bladeVelocitiesYH;};
    real* getBladeVelocitiesZ(){return this->bladeVelocitiesZH;};
    real* getBladeForcesX(){return this->bladeForcesXH;};
    real* getBladeForcesY(){return this->bladeForcesYH;};
    real* getBladeForcesZ(){return this->bladeForcesZH;};

    void setOmega(real _omega){ this->omega = _omega; };
    void setAzimuth(real _azimuth){ this->azimuth = _azimuth; };
    void setYaw(real _yaw){ this->yaw = _yaw; };
    void setBladeCoords(real* _bladeCoordsX, real* _bladeCoordsY, real* _bladeCoordsZ);
    void setBladeVelocities(real* _bladeVelocitiesX, real* _bladeVelocitiesY, real* _bladeVelocitiesZ);
    void setBladeForces(real* _bladeForcesX, real* _bladeForcesY, real* _bladeForcesZ);
    virtual void calcBladeForces();

private:
    void initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager);

    void initBladeRadii(CudaMemoryManager* cudaManager);
    void initBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeForces(CudaMemoryManager* cudaManager);
    void initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager);

    void calcForcesEllipticWing();
    void rotateBlades(real angle);

    void writeBladeCoords(uint t){};
    void writeBladeForces(uint t){};
    void writeBladeVelocities(uint t){};
    

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
    
private:
    const real density;
    real turbinePosX, turbinePosY, turbinePosZ;
    real omega, azimuth, yaw, delta_t, delta_x;
    const real diameter;
    const uint nBladeNodes;
    const uint nBlades;
    const real epsilon; // in m
    const int level;
    uint numberOfIndices;
    uint numberOfNodes;
};

#endif