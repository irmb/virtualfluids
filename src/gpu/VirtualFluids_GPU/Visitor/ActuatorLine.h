#ifndef ActuatorLine_H
#define ActuatorLine_H

#include "Visitor.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"
#include "GridGenerator/grid/GridBuilder/GridBuilder.h"

class ActuatorLine : public Visitor
{
public:
    ActuatorLine(
        const unsigned int &_nBlades,
        const real &_density,
        const unsigned int &_nBladeNodes,
        const real &_epsilon,
        real &_turbinePosX, real &_turbinePosY, real &_turbinePosZ,
        const real &_diameter,
        int &_level
    ) : nBlades(_nBlades),
        density(_density),
        nBladeNodes(_nBladeNodes), 
        epsilon(_epsilon),
        turbinePosX(_turbinePosX), turbinePosY(_turbinePosY), turbinePosZ(_turbinePosZ),
        diameter(_diameter),
        level(_level),
        Visitor()
    {
        this->numberOfNodes = this->nBladeNodes*this->nBlades;
        this->omega = 0.0f;
        this->azimuth = 0.0f;
    }

    virtual  ~ActuatorLine()
    {
        this->freeBladeCoords();
        this->freeBladeVelocities();
        this->freeBladeForces();
        this->freeBladeIndices();
        this->freeSphereIndices();
    }

    void visit(Parameter* para, int level, unsigned int t);
    void initBoundingSphere(Parameter* para, CudaMemoryManager* cudaManager);
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);

private:

    void allocBladeCoords(CudaMemoryManager* cudaManager);
    void initBladeCoords();
    void copyBladeCoordsHtoD();
    void copyBladeCoordsDtoH();
    void freeBladeCoords();
    void rotateBlades(real angle);

    void allocBladeVelocities(CudaMemoryManager* cudaManager);
    void initBladeVelocities();
    void copyBladeVelocitiesHtoD();
    void copyBladeVelocitiesDtoH();
    void freeBladeVelocities();

    void allocBladeForces(CudaMemoryManager* cudaManager);
    void initBladeForces();
    void copyBladeForcesHtoD();
    void copyBladeForcesDtoH();
    void freeBladeForces();

    void allocBladeIndices(CudaMemoryManager* cudaManager);
    void initBladeIndices(Parameter* para);
    void copyBladeIndicesHtoD();
    void freeBladeIndices();

    void allocSphereIndices(CudaMemoryManager* cudaManager);
    void copySphereIndices();
    void freeSphereIndices();
    
private:
    const real density;
    real turbinePosX, turbinePosY, turbinePosZ;
    real* bladeRadii;
    real* bladeCoordsXH, * bladeCoordsYH, * bladeCoordsZH;
    real* bladeCoordsXD, * bladeCoordsYD, * bladeCoordsZD;
    real* bladeVelocitiesXH, * bladeVelocitiesYH, * bladeVelocitiesZH;
    real* bladeVelocitiesXD, * bladeVelocitiesYD, * bladeVelocitiesZD;
    real* bladeForcesXH, * bladeForcesYH, * bladeForcesZH;
    real* bladeForcesXD, * bladeForcesYD, * bladeForcesZD;
    uint* bladeIndicesH;
    uint* bladeIndicesD; 
    int* boundingSphereIndicesH, * boundingSphereIndicesD;
    real omega, azimuth;
    const real diameter;
    const unsigned int nBladeNodes;
    const unsigned int nBlades;
    const real epsilon; // in m
    const int level;
    int numberOfIndices;
    int numberOfNodes;
    
    unsigned int* indicesBoundingBox;


};

#endif