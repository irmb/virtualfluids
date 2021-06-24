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
};

class Probe : public Visitor 
{
public:
    Probe(
        const std::string _probeName

    ):  probeName(_probeName)

    {
        
    }

    void visit(Parameter* para, int level, unsigned int t);
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager);

    void setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ);

    
private:
    const std::string probeName;
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ; 
    uint nProbePoints;

    std::vector<ProbeStruct*> probeParams;

    // int* pointIndicesH, *pointIndicesD;
    // std::vector< std::vector<int> > probeIndices;
    // std::vector< std::vector<real> > distX, distY, distZ;
};





#endif