#ifndef Probe_H
#define Probe_H

#include <cuda.h>

#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "PointerDefinitions.h"

enum class PostProcessingVariable{ 
    // HowTo add new PostProcessingVariable: Add enum here, LAST has to stay last
    // In calculatePointwiseQuantities add computation of quantity in switch statement (for point-wise statistics) 
    //  or, alternatively, in .... (for spatial statistics)
    // In writeGridFiles add lb->rw conversion factor
    // In getPostProcessingVariableNames add names
    // In isAvailableProcessingVariableNames add name in switch statement or set statement to true if the variable is 
    //  only new for the specific probe type
    // If new quantity depends on other quantities i.e. mean, catch in addPostProcessingVariable
    
    // Variables available in Point and Plane probe (all temporal pointwise statistics)
    Instantaneous,
    Means,
    Variances,

    // Variables available in PlanarAverage probe
    SpatialMeans,
    SpatioTemporalMeans,
    SpatialCovariances,
    SpatioTemporalCovariances,
    SpatialSkewness,
    SpatioTemporalSkewness,
    SpatialFlatness,
    SpatioTemporalFlatness,
    LAST,
};

struct ProbeStruct{
    uint nPoints, nIndices, nArrays, vals;
    uint *pointIndicesH, *pointIndicesD;
    real *pointCoordsX, *pointCoordsY, *pointCoordsZ;
    bool hasDistances=false;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    real *quantitiesArrayH, *quantitiesArrayD;
    bool *quantitiesH, *quantitiesD;
    uint *arrayOffsetsH, *arrayOffsetsD;
};

__global__ void calcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );

__global__ void interpAndCalcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );


class Probe : public PreCollisionInteractor 
{
public:
    Probe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut
    ):  probeName(_probeName),
        outputPath(_outputPath),
        tStartAvg(_tStartAvg),
        tAvg(_tAvg),
        tStartOut(_tStartOut),
        tOut(_tOut),
        PreCollisionInteractor()
    {
        assert("Output starts before averaging!" && tStartOut>=tStartAvg);
    }
    
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    SPtr<ProbeStruct> getProbeStruct(int level){ return this->probeParams[level]; }

    void addPostProcessingVariable(PostProcessingVariable _variable);

private:
    virtual bool isAvailablePostProcessingVariable(PostProcessingVariable _variable) = 0;

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
    const std::string outputPath;

    std::vector<SPtr<ProbeStruct>> probeParams;
    bool quantities[int(PostProcessingVariable::LAST)] = {};
    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> varNames;

    uint tStartAvg;
    uint tAvg;
    uint tStartOut;
    uint tOut;
};

#endif