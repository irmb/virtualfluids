#ifndef Probe_H
#define Probe_H

#include "Visitor.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"
#include "GridGenerator/grid/GridBuilder/GridBuilder.h"

enum class PostProcessingVariable{ 
    // Enum val is index in pointer array -> increment between enum1 and enum2 is number of quantities allocated for enum1
    // LAST is for counting total number of arrays
    // HowTo add new PostProcessingVariable: Add enum here, assign it value of LAST, assign LAST previous number+ number of quantities needed for new postProc Variable
    // In interpQuantities add computation of quantity in switch statement
    // In init add number of arrays + offset in switch statement
    // If new quantity depends on other quantities i.e. mean, catch in addPostProcessingVariable
    Means = 0,
    Variances = 4,
    LAST = 8,
};
struct ProbeStruct{
    int nPoints, nArrays;
    int *pointIndicesH, *pointIndicesD;
    real *pointCoordsX, *pointCoordsY, *pointCoordsZ;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    real* quantitiesArrayH, *quantitiesArrayD;
    PostProcessingVariable* quantitiesH,* quantitiesD; 
    int* arrayOffsetsH, *arrayOffsetsD;
};



class Probe : public Visitor 
{
public:
    Probe(
        const std::string _probeName,
        uint _tStart,
        uint _tOut

    ):  probeName(_probeName),
        tStart(_tStart),
        tOut(_tOut),
        Visitor()

    {
        
    }
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);
    void visit(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t);
    void free(Parameter* para, CudaMemoryManager* cudaManager);

    ProbeStruct* getProbeStruct(int level){ return this->probeParams[level]; }
    std::vector<PostProcessingVariable> getPostProcessingVariables(){return this->postProcessingVariables; }

    void setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ);
    void addPostProcessingVariable(PostProcessingVariable _variable);

    void write(Parameter* para, int level, int t);
    void writeCollectionFile(Parameter* para, int t);
    void writeGridFile(Parameter* para, int level, std::vector<std::string >& fnames, int t);
    std::vector<std::string> getVarNames();

    
private:
    const std::string probeName;
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ; 
    uint nProbePoints;

    std::vector<ProbeStruct*> probeParams;
    std::vector<PostProcessingVariable> postProcessingVariables;
    std::vector<std::string> fileNamesForCollectionFile;

    uint tStart;
    uint tOut;
    // int* pointIndicesH, *pointIndicesD;
    // std::vector< std::vector<int> > probeIndices;
    // std::vector< std::vector<real> > distX, distY, distZ;
};





#endif