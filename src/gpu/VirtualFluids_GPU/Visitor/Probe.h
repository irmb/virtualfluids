#ifndef Probe_H
#define Probe_H

#include "Visitor.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"
#include "GridGenerator/grid/GridBuilder/GridBuilder.h"

struct ProbeStruct{
    int nPoints;
    int *pointIndicesH, *pointIndicesD;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    std::vector<void*> quantitiesH, quantitiesD;
};

enum class PostProcessingVariable{ 
    //Enum val is index in pointer array -> increment between enum1 and enum2 is number of quantities allocated for enum1
    Means = 0,
    Variances = 4
};

class Probe : public Visitor 
{
public:
    Probe(
        const std::string _probeName

    ):  probeName(_probeName)

    {
        
    }
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void visit(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);

    ProbeStruct* getProbeStruct(int level){ return this->probeParams[level]; }

    void setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ);
    void addPostProcessingVariable(PostProcessingVariable _variable);

    
private:
    const std::string probeName;
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ; 
    uint nProbePoints;

    std::vector<ProbeStruct*> probeParams;
    std::vector<PostProcessingVariable> postProcessingVariables;
    // int* pointIndicesH, *pointIndicesD;
    // std::vector< std::vector<int> > probeIndices;
    // std::vector< std::vector<real> > distX, distY, distZ;
};





#endif