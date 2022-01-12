#ifndef Probe_H
#define Probe_H

#include <cuda.h>

#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "PointerDefinitions.h"

enum class PostProcessingVariable{ 
    // HowTo add new PostProcessingVariable: Add enum here, LAST has to stay last
    // In interpQuantities add computation of quantity in switch statement
    // In writeGridFiles add lb->rw conversion factor
    // In getPostProcessingVariableNames add names
    // If new quantity depends on other quantities i.e. mean, catch in addPostProcessingVariable
    Means,
    Variances,
    LAST,
};

struct ProbeStruct{
    uint nPoints, nArrays, vals;
    uint *pointIndicesH, *pointIndicesD;
    real *pointCoordsX, *pointCoordsY, *pointCoordsZ;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    real *quantitiesArrayH, *quantitiesArrayD;
    bool *quantitiesH, *quantitiesD;
    uint *arrayOffsetsH, *arrayOffsetsD;
};

__global__ void interpQuantities(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* distX, real* distY, real* distZ,
                                    real* distributions, 
                                    uint size_Mat, bool isEvenTimestep,
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray,
                                    bool interpolate
                                );


class Probe : public PreCollisionInteractor 
{
public:
    Probe(
        const std::string _probeName,
        uint _tStartAvg,
        uint _tStartOut,
        uint _tOut
    ):  probeName(_probeName),
        tStartAvg(_tStartAvg),
        tStartOut(_tStartOut),
        tOut(_tOut),
        PreCollisionInteractor()
    {
        assert("Output starts before averaging!" && tStartOut>=tStartAvg);
    }
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);

    SPtr<ProbeStruct> getProbeStruct(int level){ return this->probeParams[level]; }

    void addPostProcessingVariable(PostProcessingVariable _variable);

private:
    virtual void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                       std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                       std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                       int level) = 0;
    void addProbeStruct(CudaMemoryManager* cudaManager, std::vector<int>& probeIndices,
                        std::vector<real>& distX, std::vector<real>& distY, std::vector<real>& distZ,   
                        std::vector<real>& pointCoordsX, std::vector<real>& pointCoordsY, std::vector<real>& pointCoordsZ,
                        int level);
    virtual void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level) = 0;

    void write(Parameter* para, int level, int t);
    void writeCollectionFile(Parameter* para, int t);
    void writeGridFiles(Parameter* para, int level, std::vector<std::string >& fnames, int t);
    std::vector<std::string> getVarNames();
    
private:
    const std::string probeName;

    std::vector<SPtr<ProbeStruct>> probeParams;
    bool quantities[int(PostProcessingVariable::LAST)];
    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> varNames;

    uint tStartAvg;
    uint tStartOut;
    uint tOut;
};

#endif