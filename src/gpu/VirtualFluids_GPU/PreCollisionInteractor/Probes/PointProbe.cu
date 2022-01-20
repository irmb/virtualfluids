#include "PointProbe.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda/CudaGrid.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

void PointProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                       std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                       std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                       int level)
{

    real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {    
        for(uint point=0; point<this->pointCoordsX.size(); point++)
        {
            real pointCoordX = this->pointCoordsX[point];
            real pointCoordY = this->pointCoordsY[point];
            real pointCoordZ = this->pointCoordsZ[point];
            real distX = pointCoordX-para->getParH(level)->coordX_SP[j];
            real distY = pointCoordY-para->getParH(level)->coordY_SP[j];
            real distZ = pointCoordZ-para->getParH(level)->coordZ_SP[j];
            if( distX <=dx && distY <=dx && distZ <=dx &&
                distX >0.f && distY >0.f && distZ >0.f)
            {
                probeIndices_level.push_back(j);
                distX_level.push_back( distX/dx );
                distY_level.push_back( distY/dx );
                distZ_level.push_back( distZ/dx );
                pointCoordsX_level.push_back( pointCoordX );
                pointCoordsY_level.push_back( pointCoordY );
                pointCoordsZ_level.push_back( pointCoordZ );
            }
        }
    }
}

void PointProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, probeStruct->nPoints);

    interpQuantities<<<grid.grid, grid.threads>>>(  probeStruct->pointIndicesD, probeStruct->nPoints, probeStruct->vals,
                                                    probeStruct->distXD, probeStruct->distYD, probeStruct->distZD,
                                                    para->getParD(level)->d0SP.f[0], para->getParD(level)->size_Mat_SP, para->getParD(level)->evenOrOdd,
                                                    para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP, 
                                                    probeStruct->quantitiesD, probeStruct->arrayOffsetsD, probeStruct->quantitiesArrayD, true);
}

void PointProbe::addProbePointsFromList(std::vector<real>& _pointCoordsX, std::vector<real>& _pointCoordsY, std::vector<real>& _pointCoordsZ)
{
    bool isSameLength = ( (_pointCoordsX.size()==_pointCoordsY.size()) && (_pointCoordsY.size()==_pointCoordsZ.size()));
    assert("Probe: point lists have different lengths" && isSameLength);
    this->pointCoordsX.insert(this->pointCoordsX.end(), _pointCoordsX.begin(),  _pointCoordsX.end());
    this->pointCoordsY.insert(this->pointCoordsY.end(), _pointCoordsY.begin(),  _pointCoordsY.end());
    this->pointCoordsZ.insert(this->pointCoordsZ.end(), _pointCoordsZ.begin(),  _pointCoordsZ.end());
    printf("Added list of %u  points \n", uint(_pointCoordsX.size()) );
}

void PointProbe::addProbePointsFromXNormalPlane(real pos_x, real pos0_y, real pos0_z, real pos1_y, real pos1_z, uint n_y, uint n_z)
{
    int delta_y = (pos1_y-pos0_y)/(n_y-1);
    int delta_z = (pos1_z-pos0_z)/(n_z-1);

    this->pointCoordsX.reserve(this->pointCoordsX.size()+n_y*n_z);
    this->pointCoordsY.reserve(this->pointCoordsY.size()+n_y*n_z);
    this->pointCoordsZ.reserve(this->pointCoordsZ.size()+n_y*n_z);

    for(int n_y=0; n_y<n_y; n_y++)
    {
        for(int n_z=0; n_z<n_z; n_z++)
        {
            this->pointCoordsX.push_back(pos_x);
            this->pointCoordsY.push_back(pos0_y+delta_y*n_y);
            this->pointCoordsZ.push_back(pos0_z+delta_z*n_z);
        }
    }
    printf("Added %u  points \n",  n_y*n_z);

}