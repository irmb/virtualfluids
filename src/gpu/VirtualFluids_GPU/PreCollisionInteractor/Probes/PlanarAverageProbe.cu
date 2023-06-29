#include "Probe.h"
#include "PlanarAverageProbe.h"

#include <cuda/CudaGrid.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/inner_product.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "basics/constants/NumericConstants.h"

#include <algorithm>

using namespace vf::basics::constant;
///////////////////////////////////////////////////////////////////////////////////
/// Functors for thrust reductions
///////////////////////////////////////////////////////////////////////////////////

template<typename T>
struct pow2 : public thrust::unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x * x;
  }
};

template<typename T>
struct pow3 : public thrust::unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x * x * x;
  }
};

template<typename T>
struct pow4 : public thrust::unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x * x * x * x;
  }
};

struct nth_moment
{
    const float mean;
    const int n;

    nth_moment(float _mean, int _n) : mean(_mean), n(_n) {}

    __host__ __device__
        float operator()(const float& x) const { 
            
            real fluctuation = x-mean;
            real moment = fluctuation;
            for(int i = 1; i<n; i++) moment *= fluctuation;
            
            return moment;
        }
};


///////////////////////////////////////////////////////////////////////////////////

__global__ void moveIndicesInPosNormalDir( uint* pointIndices, uint nPoints, uint* neighborNormal, real* coordsX, real* coordsY, real* coordsZ )
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=nPoints) return;

    uint k = pointIndices[node];

    pointIndices[node] = neighborNormal[k];
}

__global__ void moveIndicesInNegNormalDir( uint* pointIndices, uint nPoints, uint* neighborWSB, uint* neighborInplane1, uint* neighborInplane2, real* coordsX, real* coordsY, real* coordsZ )
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=nPoints) return;

    uint k = pointIndices[node];

    pointIndices[node] = neighborWSB[neighborInplane1[neighborInplane2[k]]];
}

///////////////////////////////////////////////////////////////////////////////////

bool PlanarAverageProbe::isAvailableStatistic(Statistic _variable)
{
    bool isAvailable;

    switch (_variable)
    {
        case Statistic::Instantaneous:
        case Statistic::Means:
        case Statistic::Variances:
            isAvailable = false;
            break;
        case Statistic::SpatialMeans:
        case Statistic::SpatioTemporalMeans:
        case Statistic::SpatialCovariances:
        case Statistic::SpatioTemporalCovariances:
        case Statistic::SpatialSkewness:
        case Statistic::SpatioTemporalSkewness:
        case Statistic::SpatialFlatness:
        case Statistic::SpatioTemporalFlatness:
            isAvailable =  true;
            break;
        default:
            isAvailable =  false;
    }
    return isAvailable;
}

///////////////////////////////////////////////////////////////////////////////////
std::vector<PostProcessingVariable> PlanarAverageProbe::getPostProcessingVariables(Statistic statistic)
{
    std::vector<PostProcessingVariable> postProcessingVariables;
    switch (statistic)
    {
    case Statistic::SpatialMeans:
        postProcessingVariables.push_back( PostProcessingVariable("vx_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("nut_spatMean", this->viscosityRatio) );
        break;
    case Statistic::SpatioTemporalMeans:
        postProcessingVariables.push_back( PostProcessingVariable("vx_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("nut_spatTmpMean", this->viscosityRatio) );
        break;
    case Statistic::SpatialCovariances:
        postProcessingVariables.push_back( PostProcessingVariable("vxvx_spatMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vyvy_spatMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vzvz_spatMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vxvy_spatMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vxvz_spatMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vyvz_spatMean",  this->stressRatio) );
        break;
    case Statistic::SpatioTemporalCovariances:
        postProcessingVariables.push_back( PostProcessingVariable("vxvx_spatTmpMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vyvy_spatTmpMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vzvz_spatTmpMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vxvy_spatTmpMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vxvz_spatTmpMean",  this->stressRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vyvz_spatTmpMean",  this->stressRatio) );
        break;
    case Statistic::SpatialSkewness:
        postProcessingVariables.push_back( PostProcessingVariable("Sx_spatMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Sy_spatMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Sz_spatMean",  this->nondimensional) );
        break;
    case Statistic::SpatioTemporalSkewness:
        postProcessingVariables.push_back( PostProcessingVariable("Sx_spatTmpMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Sy_spatTmpMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Sz_spatTmpMean",  this->nondimensional) );
        break;
    case Statistic::SpatialFlatness:
        postProcessingVariables.push_back( PostProcessingVariable("Fx_spatMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Fy_spatMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Fz_spatMean",  this->nondimensional) );
        break;
    case Statistic::SpatioTemporalFlatness:
        postProcessingVariables.push_back( PostProcessingVariable("Fx_spatTmpMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Fy_spatTmpMean",  this->nondimensional) );
        postProcessingVariables.push_back( PostProcessingVariable("Fz_spatTmpMean",  this->nondimensional) );
        break;

    default:
        throw std::runtime_error("PlanarAverageProbe::getPostProcessingVariables: Statistic unavailable!");
        break;
    }
    return postProcessingVariables;
}

///////////////////////////////////////////////////////////////////////////////////

void PlanarAverageProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                            std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                            std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                            int level)
{
    real dx = abs(para->getParH(level)->coordinateX[1]-para->getParH(level)->coordinateX[para->getParH(level)->neighborX[1]]);
    
    real /* *pointCoordsInplane1_par, *pointCoordsInplane2_par,*/ *pointCoordsNormal_par;
    std::vector<real> *pointCoordsInplane1, *pointCoordsInplane2, *pointCoordsNormal;
    
    if(this->planeNormal == 'x'){  
                                    pointCoordsNormal       = &pointCoordsX_level; 
                                    pointCoordsInplane1     = &pointCoordsY_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordinateX; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordY_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'y'){  
                                    pointCoordsNormal       = &pointCoordsY_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordinateY; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'z'){  
                                    pointCoordsNormal       = &pointCoordsZ_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsY_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordinateZ; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordY_SP;
                                }

    // Find all points along the normal direction
    for(size_t j = 1; j < para->getParH(level)->numberOfNodes; j++ )
    {
        if(para->getParH(level)->typeOfGridNode[j] == GEO_FLUID)
        {   
            if( std::find(pointCoordsNormal->begin(), pointCoordsNormal->end(), pointCoordsNormal_par[j]) == pointCoordsNormal->end())  
            {
                pointCoordsNormal->push_back( pointCoordsNormal_par[j] );
                pointCoordsInplane1->push_back(999999.);
                pointCoordsInplane2->push_back(999999.);
            }
        }
    }
    std::sort(pointCoordsNormal->begin(), pointCoordsNormal->end());
    
    // Find all pointCoords in the first plane 
    for(size_t pos = 1; pos < para->getParH(level)->numberOfNodes; pos++ )
    {
        if( para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID && pointCoordsNormal_par[pos] == pointCoordsNormal->at(0)) 
        {
            //not needed in current state, might become relevant for two-point correlations
            // pointCoordsNormal->push_back( pointCoordsNormal_par[j] ); 
            // pointCoordsInplane1->push_back( pointCoordsInplane1_par[j] );
            // pointCoordsInplane2->push_back( pointCoordsInplane2_par[j] );

            probeIndices_level.push_back((int)pos);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

void PlanarAverageProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t_level, int level)
{   
    // Compute macroscopic variables in entire domain
    CalcMacCompSP27(
        para->getParD(level)->velocityX, 
        para->getParD(level)->velocityY, 
        para->getParD(level)->velocityZ,
        para->getParD(level)->rho, 
        para->getParD(level)->pressure, 
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, 
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->numberofthreads, 
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("In PlanarAverageProbe Kernel CalcMacSP27 execution failed");

    // Definition of normal and inplane directions for moveIndices kernels
    uint *neighborNormal, *neighborInplane1, *neighborInplane2;
    if( this->planeNormal == 'x' )
    {
        neighborNormal   = para->getParD(level)->neighborX;
        neighborInplane1 = para->getParD(level)->neighborY;
        neighborInplane2 = para->getParD(level)->neighborZ;
    }
    if( this->planeNormal == 'y' )
    {
        neighborNormal   = para->getParD(level)->neighborY;
        neighborInplane1 = para->getParD(level)->neighborX;
        neighborInplane2 = para->getParD(level)->neighborZ;
    }
    if( this->planeNormal == 'z' )
    {
        neighborNormal   = para->getParD(level)->neighborZ;
        neighborInplane1 = para->getParD(level)->neighborX;
        neighborInplane2 = para->getParD(level)->neighborY;
    }

    bool doTmpAveraging = t_level>=(this->getTStartTmpAveraging()*exp2(level));

    // Pointer casts to use device arrays in thrust reductions
    thrust::device_ptr<uint> indices_thrust = thrust::device_pointer_cast(probeStruct->pointIndicesD);
    thrust::device_ptr<real> vx_thrust = thrust::device_pointer_cast(para->getParD(level)->velocityX);
    thrust::device_ptr<real> vy_thrust = thrust::device_pointer_cast(para->getParD(level)->velocityY);
    thrust::device_ptr<real> vz_thrust = thrust::device_pointer_cast(para->getParD(level)->velocityZ);
    thrust::device_ptr<real> nut_thrust = thrust::device_pointer_cast(para->getParD(level)->turbViscosity);

    real N = (real)probeStruct->nIndices;
    real invNumberOfTimestepsInTmpAvg = c1o1/real(probeStruct->timestepInTimeAverage+1);
    uint nPoints = probeStruct->nPoints;
    // Permutation iterators for direct iteration over the velocities of the planes
    typedef thrust::device_vector<real>::iterator valIterator;
    typedef thrust::device_vector<uint>::iterator indIterator;
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_begin(vx_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_end  (vx_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_begin(vy_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_end  (vy_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_begin(vz_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_end  (vz_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> nut_iter_begin(nut_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> nut_iter_end  (nut_thrust, indices_thrust+probeStruct->nIndices);

    for( uint i=0; i<nPoints; i++ )
    {
        uint node = probeStruct->isEvenTAvg? i : nPoints-1-i; // Note, loop moves in positive normal dir at even calls and in negative normal dir in odd calls

        if(probeStruct->quantitiesH[int(Statistic::SpatialMeans)])
        {
            // Compute the instantaneous spatial means of the velocity moments 
            real spatMean_vx = thrust::reduce(vx_iter_begin, vx_iter_end)/N;
            real spatMean_vy = thrust::reduce(vy_iter_begin, vy_iter_end)/N;
            real spatMean_vz = thrust::reduce(vz_iter_begin, vz_iter_end)/N;

            real spatMean_nut;
            if(para->getUseTurbulentViscosity()) spatMean_nut = thrust::reduce(nut_iter_begin, nut_iter_end)/N;

            uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialMeans)];
            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_vx;
            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_vy;
            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_vz;
            if(para->getUseTurbulentViscosity()) probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] = spatMean_nut;

            if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalMeans)] && doTmpAveraging)
            {
                uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalMeans)];
                real spatTmpMean_vx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                real spatTmpMean_vy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                real spatTmpMean_vz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];
                real spatTmpMean_nut_old;
                if(para->getUseTurbulentViscosity()) spatTmpMean_nut_old = probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node];;

                probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_vx-spatTmpMean_vx_old)*invNumberOfTimestepsInTmpAvg;
                probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_vy-spatTmpMean_vy_old)*invNumberOfTimestepsInTmpAvg;
                probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_vz-spatTmpMean_vz_old)*invNumberOfTimestepsInTmpAvg;
                if(para->getUseTurbulentViscosity()) probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] += (spatMean_nut-spatTmpMean_nut_old)*invNumberOfTimestepsInTmpAvg;

            }
        
            if(probeStruct->quantitiesH[int(Statistic::SpatialCovariances)])
            {   // <u_i' u_j'> = <u_i u_j> - <u_i>*<u_i> 
                real vx2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vy2 = thrust::transform_reduce(vy_iter_begin, vy_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vz2 = thrust::transform_reduce(vz_iter_begin, vz_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vxvy = thrust::inner_product(vx_iter_begin, vx_iter_end, vy_iter_begin, 0.f)/N;
                real vxvz = thrust::inner_product(vx_iter_begin, vx_iter_end, vz_iter_begin, 0.f)/N;
                real vyvz = thrust::inner_product(vy_iter_begin, vy_iter_end, vz_iter_begin, 0.f)/N;
                real spatMean_vxvx = vx2-spatMean_vx*spatMean_vx;
                real spatMean_vyvy = vy2-spatMean_vy*spatMean_vy;
                real spatMean_vzvz = vz2-spatMean_vz*spatMean_vz;
                real spatMean_vxvy = vxvy-spatMean_vx*spatMean_vy;
                real spatMean_vxvz = vxvz-spatMean_vx*spatMean_vz;
                real spatMean_vyvz = vyvz-spatMean_vy*spatMean_vz;

                uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialCovariances)];
                probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_vxvx;
                probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_vyvy;
                probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_vzvz;
                probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] = spatMean_vxvy;
                probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node] = spatMean_vxvz;
                probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node] = spatMean_vyvz;

                if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalCovariances)] && doTmpAveraging)
                {
                    uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalCovariances)];
                    real spatTmpMean_vxvx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                    real spatTmpMean_vyvy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                    real spatTmpMean_vzvz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];
                    real spatTmpMean_vxvy_old = probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node];
                    real spatTmpMean_vxvz_old = probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node];
                    real spatTmpMean_vyvz_old = probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node];

                    probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_vxvx-spatTmpMean_vxvx_old)*invNumberOfTimestepsInTmpAvg;
                    probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_vyvy-spatTmpMean_vyvy_old)*invNumberOfTimestepsInTmpAvg;
                    probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_vzvz-spatTmpMean_vzvz_old)*invNumberOfTimestepsInTmpAvg;
                    probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] += (spatMean_vxvy-spatTmpMean_vxvy_old)*invNumberOfTimestepsInTmpAvg;
                    probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node] += (spatMean_vxvz-spatTmpMean_vxvz_old)*invNumberOfTimestepsInTmpAvg;
                    probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node] += (spatMean_vyvz-spatTmpMean_vyvz_old)*invNumberOfTimestepsInTmpAvg;
                }

                if(probeStruct->quantitiesH[int(Statistic::SpatialSkewness)])
                {   // <u_i'^3> = <u_i^3> - <u_i>^3 - 3 <u_i> <u_i'^2>
                    // real vx3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    // real vy3 = thrust::transform_reduce(vy_iter_begin, vy_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    // real vz3 = thrust::transform_reduce(vz_iter_begin, vz_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    real spatMean_vxvxvx = thrust::transform_reduce(vx_iter_begin, vx_iter_end, nth_moment(spatMean_vx, 3), 0.f, thrust::plus<real>())/N; 
                    //vx3 - spatMean_vx*spatMean_vx*spatMean_vx - 3*spatMean_vx*spatMean_vxvx; -> alternative only using vx3, etc. but containing some bug. Potentially better in terms of round-off errors.
                    real spatMean_vyvyvy = thrust::transform_reduce(vy_iter_begin, vy_iter_end, nth_moment(spatMean_vy, 3), 0.f, thrust::plus<real>())/N; 
                    //vy3 - spatMean_vy*spatMean_vy*spatMean_vy - 3*spatMean_vy*spatMean_vzvz;
                    real spatMean_vzvzvz = thrust::transform_reduce(vz_iter_begin, vz_iter_end, nth_moment(spatMean_vz, 3), 0.f, thrust::plus<real>())/N; 
                    //vz3 - spatMean_vz*spatMean_vz*spatMean_vz - 3*spatMean_vz*spatMean_vzvz;
                    real spatMean_Sx = spatMean_vxvxvx/pow(spatMean_vxvx, 1.5f);
                    real spatMean_Sy = spatMean_vyvyvy/pow(spatMean_vyvy, 1.5f);
                    real spatMean_Sz = spatMean_vzvzvz/pow(spatMean_vzvz, 1.5f);

                    uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialSkewness)];
                    probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_Sx;
                    probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_Sy;
                    probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_Sz;

                    if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalSkewness)] && doTmpAveraging)
                    {
                        uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalSkewness)];
                        real spatTmpMean_Sx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                        real spatTmpMean_Sy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                        real spatTmpMean_Sz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];

                        probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_Sx-spatTmpMean_Sx_old)*invNumberOfTimestepsInTmpAvg;
                        probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_Sy-spatTmpMean_Sy_old)*invNumberOfTimestepsInTmpAvg;
                        probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_Sz-spatTmpMean_Sz_old)*invNumberOfTimestepsInTmpAvg;
                    }

                    if(probeStruct->quantitiesH[int(Statistic::SpatialFlatness)])
                    {   // <u_i'^4> = <u_i^4> - <u_i>^4 - 6 <u_i>^2 <u_i'^2> - 4 <u> <u'^3>
                        // real vx4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        // real vy4 = thrust::transform_reduce(vy_iter_begin, vy_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        // real vz4 = thrust::transform_reduce(vz_iter_begin, vz_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        real spatMean_vxvxvxvx = thrust::transform_reduce(vx_iter_begin, vx_iter_end, nth_moment(spatMean_vx, 4), 0.f, thrust::plus<real>())/N; //vx4 - spatMean_vx*spatMean_vx*spatMean_vx*spatMean_vx - 6*spatMean_vx*spatMean_vx*vx2 - 4*spatMean_vx*vx3;
                        real spatMean_vyvyvyvy = thrust::transform_reduce(vy_iter_begin, vy_iter_end, nth_moment(spatMean_vy, 4), 0.f, thrust::plus<real>())/N; //vy4 - spatMean_vy*spatMean_vy*spatMean_vy*spatMean_vy - 6*spatMean_vy*spatMean_vx*vy2 - 4*spatMean_vy*vy3;
                        real spatMean_vzvzvzvz = thrust::transform_reduce(vz_iter_begin, vz_iter_end, nth_moment(spatMean_vz, 4), 0.f, thrust::plus<real>())/N; //vz4 - spatMean_vz*spatMean_vz*spatMean_vz*spatMean_vz - 6*spatMean_vz*spatMean_vx*vz2 - 4*spatMean_vz*vz3;
                        real spatMean_Fx = spatMean_vxvxvxvx/(spatMean_vxvx*spatMean_vxvx);
                        real spatMean_Fy = spatMean_vyvyvyvy/(spatMean_vyvy*spatMean_vyvy);
                        real spatMean_Fz = spatMean_vzvzvzvz/(spatMean_vzvz*spatMean_vzvz);

                        uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialFlatness)];
                        probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_Fx;
                        probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_Fy;
                        probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_Fz;

                        if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalFlatness)] && doTmpAveraging)
                        {
                            uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalFlatness)];
                            real spatTmpMean_Fx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                            real spatTmpMean_Fy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                            real spatTmpMean_Fz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];

                            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_Fx-spatTmpMean_Fx_old)*invNumberOfTimestepsInTmpAvg;
                            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_Fy-spatTmpMean_Fy_old)*invNumberOfTimestepsInTmpAvg;
                            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_Fz-spatTmpMean_Fz_old)*invNumberOfTimestepsInTmpAvg;
                        }
                    }
                }
        }
    }
        if(i<probeStruct->nPoints-1)
        {
            vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, probeStruct->nIndices);
            if(probeStruct->isEvenTAvg) 
                moveIndicesInPosNormalDir<<<grid.grid, grid.threads>>>( probeStruct->pointIndicesD, probeStruct->nIndices, neighborNormal, para->getParD(level)->coordinateX, para->getParD(level)->coordinateY, para->getParD(level)->coordinateZ );
            else 
                moveIndicesInNegNormalDir<<<grid.grid, grid.threads>>>( probeStruct->pointIndicesD, probeStruct->nIndices, para->getParD(level)->neighborInverse, neighborInplane1, neighborInplane2, para->getParD(level)->coordinateX, para->getParD(level)->coordinateY, para->getParD(level)->coordinateZ ); 
        } 
    }
    probeStruct->isEvenTAvg=!(probeStruct->isEvenTAvg);

    getLastCudaError("PlanarAverageProbe::calculateQuantities execution failed");
}
