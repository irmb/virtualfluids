#ifndef Probe_H
#define Probe_H

#include "Visitor.h"
#include "PointerDefinitions.h"


enum class PostProcessingVariable{ 
    // LAST is for counting total number of arrays
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


class Probe : public Visitor 
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
        Visitor()
    {
        assert("Output starts before averaging!" && tStartOut>=tStartAvg);
    }

    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void visit(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);

    ProbeStruct* getProbeStruct(int level){ return this->probeParams[level]; }

    void setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ);
    void setProbePointsFromXNormalPlane(real pos_x, real pos0_y, real pos0_z, real pos1_y, real pos1_z, real delta_y, real delta_z);
    void addPostProcessingVariable(PostProcessingVariable _variable);

    void write(Parameter* para, int level, int t);
    void writeCollectionFile(Parameter* para, int t);
    void writeGridFiles(Parameter* para, int level, std::vector<std::string >& fnames, int t);
    std::vector<std::string> getVarNames();

    
private:
    const std::string probeName;
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ; 
    uint nProbePoints;

    std::vector<ProbeStruct*> probeParams;
    bool quantities[int(PostProcessingVariable::LAST)];
    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> varNames;

    uint tStartAvg;
    uint tStartOut;
    uint tOut;
};

#endif