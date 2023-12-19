#include "Probe.h"
#include "PointProbe.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Cuda/CudaMemoryManager.h"

bool PointProbe::isAvailableStatistic(Statistic _variable)
{
    bool isAvailable;
    switch (_variable)
    {
        case Statistic::Instantaneous:
        case Statistic::Means:
        case Statistic::Variances:
            isAvailable = true;
            break;
        case Statistic::SpatialMeans:
        case Statistic::SpatioTemporalMeans:
        case Statistic::SpatialCovariances:
        case Statistic::SpatioTemporalCovariances:
        case Statistic::SpatialSkewness:
        case Statistic::SpatioTemporalSkewness:
        case Statistic::SpatialFlatness:
        case Statistic::SpatioTemporalFlatness:
            isAvailable = false;
            break;
        default:
            isAvailable = false;
    }
    return isAvailable;
}

std::vector<PostProcessingVariable> PointProbe::getPostProcessingVariables(Statistic statistic)
{
    std::vector<PostProcessingVariable> postProcessingVariables;
    switch (statistic)
    {
    case Statistic::Instantaneous:
        postProcessingVariables.push_back( PostProcessingVariable("vx",  velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("rho", this->densityRatio ) );
        break;
    case Statistic::Means:
        postProcessingVariables.push_back( PostProcessingVariable("vx_mean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_mean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_mean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("rho_mean", this->densityRatio ) );
        break;
    case Statistic::Variances:
        postProcessingVariables.push_back( PostProcessingVariable("vx_var",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_var",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_var",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("rho_var", this->densityRatio) );
        break;

    default:
        throw std::runtime_error("PointProbe::getPostProcessingVariables: Statistic unavailable!");
        break;
    }
    return postProcessingVariables;
}

void PointProbe::findPoints(Parameter* para, std::vector<int>& probeIndices_level,
                       std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                       std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                       int level)
{

    real dx = abs(para->getParH(level)->coordinateX[1]-para->getParH(level)->coordinateX[para->getParH(level)->neighborX[1]]);
    for(size_t pos = 1; pos < para->getParH(level)->numberOfNodes; pos++ )
    {    
        for(uint point=0; point<this->pointCoordsX.size(); point++)
        {
            real pointCoordX = this->pointCoordsX[point];
            real pointCoordY = this->pointCoordsY[point];
            real pointCoordZ = this->pointCoordsZ[point];
            real distX = pointCoordX-para->getParH(level)->coordinateX[pos];
            real distY = pointCoordY-para->getParH(level)->coordinateY[pos];
            real distZ = pointCoordZ-para->getParH(level)->coordinateZ[pos];
            if( distX <=dx && distY <=dx && distZ <=dx &&
                distX >0.f && distY >0.f && distZ >0.f)
            {
                probeIndices_level.push_back((int)pos);
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

void PointProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, probeStruct->nPoints);
    int oldTimestepInTimeseries = this->outputTimeSeries ? calcOldTimestep(probeStruct->timestepInTimeseries, probeStruct->lastTimestepInOldTimeseries) : 0;
    int currentTimestep = this->outputTimeSeries ? probeStruct->timestepInTimeseries : 0;
    interpAndCalcQuantitiesKernel<<<grid.grid, grid.threads>>>(  probeStruct->pointIndicesD, probeStruct->nPoints, oldTimestepInTimeseries, currentTimestep, probeStruct->timestepInTimeAverage, probeStruct->nTimesteps,
                                                probeStruct->distXD, probeStruct->distYD, probeStruct->distZD,
                                                para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ, para->getParD(level)->rho, 
                                                para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ, 
                                                probeStruct->quantitiesD, probeStruct->arrayOffsetsD, probeStruct->quantitiesArrayD);
}

void PointProbe::addProbePoint(real pointCoordX, real pointCoordY, real pointCoordZ)
{
    this->pointCoordsX.push_back(pointCoordX);
    this->pointCoordsY.push_back(pointCoordY);
    this->pointCoordsZ.push_back(pointCoordZ);
}

void PointProbe::addProbePointsFromList(std::vector<real>& _pointCoordsX, std::vector<real>& _pointCoordsY, std::vector<real>& _pointCoordsZ)
{
    bool isSameLength = ( (_pointCoordsX.size()==_pointCoordsY.size()) && (_pointCoordsY.size()==_pointCoordsZ.size()));
    if (!isSameLength) throw std::runtime_error("Probe::addProbePointsFromList(): point lists have different lengths!");
    this->pointCoordsX.insert(this->pointCoordsX.end(), _pointCoordsX.begin(),  _pointCoordsX.end());
    this->pointCoordsY.insert(this->pointCoordsY.end(), _pointCoordsY.begin(),  _pointCoordsY.end());
    this->pointCoordsZ.insert(this->pointCoordsZ.end(), _pointCoordsZ.begin(),  _pointCoordsZ.end());
    printf("Added list of %u  points \n", uint(_pointCoordsX.size()) );
}

void PointProbe::getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);
        std::vector<uint> probeIndices( probeStruct->pointIndicesH, probeStruct->pointIndicesH+probeStruct->nIndices);
        gridProvider->tagFluidNodeIndices( probeIndices, CollisionTemplate::WriteMacroVars, level);
    }
}
