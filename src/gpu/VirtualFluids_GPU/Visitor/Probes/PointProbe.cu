#include "PointProbe.h"

#include "Kernel/Utilities/CudaGrid.h"

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
        for(uint point=0; point<this->nProbePoints; point++)
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
                // printf("x %f y %f z %f", pointCoordX, pointCoordY, pointCoordZ);
            }
        }
    }
}

void PointProbe::calculateQuantities(ProbeStruct* probeStruct, Parameter* para, int level)
{
    vf::gpu::CudaGrid grid = vf::gpu::CudaGrid(128, probeStruct->nPoints);

    interpQuantities<<<grid.grid, grid.threads>>>(  probeStruct->pointIndicesD, probeStruct->nPoints, probeStruct->vals,
                                                    probeStruct->distXD, probeStruct->distYD, probeStruct->distZD,
                                                    para->getParD(level)->vx_SP, para->getParD(level)->vy_SP, para->getParD(level)->vz_SP, para->getParD(level)->rho_SP, 
                                                    para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP, 
                                                    probeStruct->quantitiesD, probeStruct->arrayOffsetsD, probeStruct->quantitiesArrayD, true);
}

void PointProbe::setProbePointsFromList(std::vector<real>& _pointCoordsX, std::vector<real>& _pointCoordsY, std::vector<real>& _pointCoordsZ)
{
    bool isSameLength = ( (_pointCoordsX.size()==_pointCoordsY.size()) && (_pointCoordsY.size()==_pointCoordsZ.size()));
    assert("Probe: point lists have different lengths" && isSameLength);
    this->pointCoordsX = _pointCoordsX;
    this->pointCoordsY = _pointCoordsY;
    this->pointCoordsZ = _pointCoordsZ;
    this->nProbePoints = uint(_pointCoordsX.size());
    printf("Added list of %u  points \n", this->nProbePoints );
}

void PointProbe::setProbePointsFromXNormalPlane(real pos_x, real pos0_y, real pos0_z, real pos1_y, real pos1_z, real delta_y, real delta_z)
{
    int n_points_y = int((pos1_y-pos0_y)/delta_y);
    int n_points_z = int((pos1_z-pos0_z)/delta_z);

    std::vector<real> pointCoordsXtmp, pointCoordsYtmp, pointCoordsZtmp;

    for(int n_y=0; n_y<n_points_y; n_y++)
    {
        for(int n_z=0; n_z<n_points_z; n_z++)
        {
            pointCoordsXtmp.push_back(pos_x);
            pointCoordsYtmp.push_back(pos0_y+delta_y*n_y);
            pointCoordsZtmp.push_back(pos0_z+delta_z*n_z);
        }
    }
    this->setProbePointsFromList(pointCoordsXtmp, pointCoordsYtmp, pointCoordsZtmp);
}