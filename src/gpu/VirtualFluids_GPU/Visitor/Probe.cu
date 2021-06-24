#include "Probe.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>


void Probe::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{

    probeParams.resize(para->getMaxLevel()+1);

    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        std::vector<int> probeIndices_level;
        std::vector<real> distX_level;
        std::vector<real> distY_level;
        std::vector<real> distZ_level;
        real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
        for(uint j=0; j<para->getParH(level)->size_Mat_SP; j++ )
        {    
            for(uint point=0; point<this->nProbePoints; point++)
            {
                real distX = this->pointCoordsX[point]-para->getParH(level)->coordX_SP[j];
                real distY = this->pointCoordsY[point]-para->getParH(level)->coordY_SP[j];
                real distZ = this->pointCoordsZ[point]-para->getParH(level)->coordZ_SP[j];
                if( distX <=dx && distY <=dx && distZ <=dx &&
                    distX >0.f && distY >0.f && distZ >0.f)
                {
                    probeIndices_level.push_back(j);
                    distX_level.push_back( distX/dx );
                    distY_level.push_back( distY/dx );
                    distZ_level.push_back( distZ/dx );
                    // printf("Found Point %i, x: %f, y: %f, z: %f, \n For %f %f %f, \n distx: %f, disty: %f, distz: %f \n", j, para->getParH(level)->coordX_SP[j],para->getParH(level)->coordY_SP[j],para->getParH(level)->coordZ_SP[j],
                    // this->pointCoordsX[point], this->pointCoordsY[point], this->pointCoordsZ[point], 
                    // distX, distY, distZ);
                }
            }
        }

        probeParams[level] = new ProbeStruct;
        probeParams[level]->nPoints = probeIndices_level.size();
        // Might have to catch nPoints=0 ?!?!
        checkCudaErrors( cudaMallocHost((void**) &probeParams[level]->distXH,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMallocHost((void**) &probeParams[level]->distYH,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMallocHost((void**) &probeParams[level]->distZH,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMallocHost((void**) &probeParams[level]->pointIndicesH, sizeof(int)*probeParams[level]->nPoints) );

        checkCudaErrors( cudaMalloc    ((void**) &probeParams[level]->distXD,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMalloc    ((void**) &probeParams[level]->distYD,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMalloc    ((void**) &probeParams[level]->distZD,        sizeof(real)*probeParams[level]->nPoints) );
        checkCudaErrors( cudaMalloc    ((void**) &probeParams[level]->pointIndicesD, sizeof(int)*probeParams[level]->nPoints) );

        std::copy(distX_level.begin(), distX_level.end(), probeParams[level]->distXH);
        std::copy(distY_level.begin(), distY_level.end(), probeParams[level]->distYH);
        std::copy(distZ_level.begin(), distZ_level.end(), probeParams[level]->distZH);
        std::copy(probeIndices_level.begin(), probeIndices_level.end(), probeParams[level]->pointIndicesH);
    }
}


void Probe::visit(Parameter* para, int level, unsigned int t)
{

}

void Probe::setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ)
{
    bool isSameLength = ( (_pointCoordsX.size()==_pointCoordsY.size()) && (_pointCoordsY.size()==_pointCoordsZ.size()));
    assert("Probe: point lists have different lengths" && isSameLength);
    this->pointCoordsX = _pointCoordsX;
    this->pointCoordsY = _pointCoordsY;
    this->pointCoordsZ = _pointCoordsZ;
    this->nProbePoints = _pointCoordsX.size();
    printf("Adde list of %u  points", this->nProbePoints );
}