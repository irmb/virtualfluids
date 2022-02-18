#include "PlanarAverageProbe.h"

#include <cuda/CudaGrid.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

void PlanarAverageProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                            std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                            std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                            int level)
{
    real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {
        real pointCoordX = para->getParH(level)->coordX_SP[j];
        real pointCoordY = para->getParH(level)->coordY_SP[j];
        real pointCoordZ = para->getParH(level)->coordZ_SP[j];
        real distX = pointCoordX - this->posX;
        real distY = pointCoordY - this->posY;
        real distZ = pointCoordZ - this->posZ;

        if( distX <= this->deltaX && distY <= this->deltaY && distZ <= this->deltaZ &&
            distX >=0.f && distY >=0.f && distZ >=0.f)
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

void PlanarAverageProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, probeStruct->nPoints);
    interpQuantities<<<grid.grid, grid.threads>>>(  probeStruct->pointIndicesD, probeStruct->nPoints, probeStruct->vals,
                                                    probeStruct->distXD, probeStruct->distYD, probeStruct->distZD,
                                                    para->getParD(level)->vx_SP, para->getParD(level)->vy_SP, para->getParD(level)->vz_SP, para->getParD(level)->rho_SP, 
                                                    para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP, 
                                                    probeStruct->quantitiesD, probeStruct->arrayOffsetsD, probeStruct->quantitiesArrayD, false);

}