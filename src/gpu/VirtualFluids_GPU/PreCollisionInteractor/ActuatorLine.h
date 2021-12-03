#ifndef ActuatorLine_H
#define ActuatorLine_H

#include "PreCollisionInteractor.h"
#include "PointerDefinitions.h"

class Parameter;
class GridProvider;

class ActuatorLine : public PreCollisionInteractor
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

    }

    virtual  ~ActuatorLine()
    {
        
    }

    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void visit(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);
    void write(uint t);

    uint getNBladeNodes(){return this->nBladeNodes;};
    uint getNBlades(){return this->nBlades;};
    uint getNumberOfIndices(){return this->numberOfIndices;};
    uint getNumberOfNodes(){return this->numberOfNodes;};
    real* getBladeCoordsX(){return this->bladeCoordsXH;};
    real* getBladeCoordsY(){return this->bladeCoordsYH;};
    real* getBladeCoordsZ(){return this->bladeCoordsZH;};
    real* getBladeVelocitiesX(){return this->bladeVelocitiesXH;};
    real* getBladeVelocitiesY(){return this->bladeVelocitiesYH;};
    real* getBladeVelocitiesZ(){return this->bladeVelocitiesZH;};
    real* getBladeForcesX(){return this->bladeForcesXH;};
    real* getBladeForcesY(){return this->bladeForcesYH;};
    real* getBladeForcesZ(){return this->bladeForcesZH;};
    void setBladeCoordsX(real* _bladeCoordsX){ this->bladeCoordsXH = _bladeCoordsX;};
    void setBladeCoordsY(real* _bladeCoordsY){ this->bladeCoordsYH = _bladeCoordsY;};
    void setBladeCoordsZ(real* _bladeCoordsZ){ this->bladeCoordsZH = _bladeCoordsZ;};
    void setBladeVelocitiesX(real* _bladeVelocitiesX){ this->bladeVelocitiesXH = _bladeVelocitiesX;};
    void setBladeVelocitiesY(real* _bladeVelocitiesY){ this->bladeVelocitiesYH = _bladeVelocitiesY;};
    void setBladeVelocitiesZ(real* _bladeVelocitiesZ){ this->bladeVelocitiesZH = _bladeVelocitiesZ;};
    void setBladeForcesX(real* _bladeForcesX){ this->bladeForcesXH = _bladeForcesX;};
    void setBladeForcesY(real* _bladeForcesY){ this->bladeForcesYH = _bladeForcesY;};
    void setBladeForcesZ(real* _bladeForcesZ){ this->bladeForcesZH = _bladeForcesZ;};

private:
    void initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager);

    void initBladeRadii(CudaMemoryManager* cudaManager);
    void initBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeForces(CudaMemoryManager* cudaManager);
    void initBladeIndices(Parameter* para, CudaMemoryManager* cudaManager);

    void calcForcesEllipticWing();
    void rotateBlades(real angle);
    void calcBladeForces();

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
    
private:
    const real density;
    real turbinePosX, turbinePosY, turbinePosZ;
    real omega, azimuth, delta_t, delta_x;
    const real diameter;
    const uint nBladeNodes;
    const uint nBlades;
    const real epsilon; // in m
    const int level;
    uint numberOfIndices;
    uint numberOfNodes;
    };

#endif