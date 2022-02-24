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

bool PlanarAverageProbe::isAvailablePostProcessingVariable(PostProcessingVariable _variable)
{
    bool isAvailable = false;

    switch (_variable)
    {
    case PostProcessingVariable::Instantaneous:
        isAvailable = false;
        break;
    case PostProcessingVariable::Means:
        isAvailable = false;
        break;
    case PostProcessingVariable::Variances:
        isAvailable = false;
        break;
    case PostProcessingVariable::SpatialMeans:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatioTemporalMeans:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatialCovariances:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatioTemporalCovariances:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatialSkewness:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatioTemporalSkewness:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatialFlatness:
        isAvailable = true;
        break;
    case PostProcessingVariable::SpatioTemporalFlatness:
        isAvailable = true;
        break;
    }
    return isAvailable;
}


void PlanarAverageProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                            std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                            std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                            int level)
{
    real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
    
    real /* *pointCoordsInplane1_par, *pointCoordsInplane2_par,*/ *pointCoordsNormal_par;
    std::vector<real> *pointCoordsInplane1, *pointCoordsInplane2, *pointCoordsNormal;
    
    if(this->planeNormal == 'x'){  
                                    pointCoordsNormal       = &pointCoordsX_level; 
                                    pointCoordsInplane1     = &pointCoordsY_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordX_SP; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordY_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'y'){  
                                    pointCoordsNormal       = &pointCoordsY_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordY_SP; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'z'){  
                                    pointCoordsNormal       = &pointCoordsZ_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsY_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordZ_SP; 
                                    // pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    // pointCoordsInplane2_par = para->getParH(level)->coordY_SP;
                                }

    // Find all points along the normal direction
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {
        if(para->getParH(level)->geoSP[j] == GEO_FLUID)
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
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {
        if( para->getParH(level)->geoSP[j] == GEO_FLUID && pointCoordsNormal_par[j] == pointCoordsNormal->at(0)) 
        {
            //not needed in current state, might become relevant for two-point correlations
            // pointCoordsNormal->push_back( pointCoordsNormal_par[j] ); 
            // pointCoordsInplane1->push_back( pointCoordsInplane1_par[j] );
            // pointCoordsInplane2->push_back( pointCoordsInplane2_par[j] );

            probeIndices_level.push_back(j);
        }
    }
}

void PlanarAverageProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level)
{   
    // Definition of normal and inplane directions for moveIndices kernels
    uint *neighborNormal, *neighborInplane1, *neighborInplane2;
    if( this->planeNormal == 'x' )
    {
        neighborNormal   = para->getParD(level)->neighborX_SP;
        neighborInplane1 = para->getParD(level)->neighborY_SP;
        neighborInplane2 = para->getParD(level)->neighborZ_SP;
    }
    if( this->planeNormal == 'y' )
    {
        neighborNormal   = para->getParD(level)->neighborY_SP;
        neighborInplane1 = para->getParD(level)->neighborX_SP;
        neighborInplane2 = para->getParD(level)->neighborZ_SP;
    }
    if( this->planeNormal == 'z' )
    {
        neighborNormal   = para->getParD(level)->neighborZ_SP;
        neighborInplane1 = para->getParD(level)->neighborX_SP;
        neighborInplane2 = para->getParD(level)->neighborY_SP;
    }

    // Pointer casts to use device arrays in thrust reductions
    thrust::device_ptr<uint> indices_thrust = thrust::device_pointer_cast(probeStruct->pointIndicesD);
    thrust::device_ptr<real> vx_thrust = thrust::device_pointer_cast(para->getParD(level)->vx_SP);
    thrust::device_ptr<real> vy_thrust = thrust::device_pointer_cast(para->getParD(level)->vy_SP);
    thrust::device_ptr<real> vz_thrust = thrust::device_pointer_cast(para->getParD(level)->vz_SP);

    real N = (real)probeStruct->nIndices;

    // Permutation iterators for direct iteration over the velocities of the planes
    typedef thrust::device_vector<real>::iterator valIterator;
    typedef thrust::device_vector<uint>::iterator indIterator;
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_begin(vx_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_end  (vx_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_begin(vy_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_end  (vy_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_begin(vz_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_end  (vz_thrust, indices_thrust+probeStruct->nIndices);


    printf("isEven: %u \n ", this->isEvenTAvg);
    for( int i=0; i<probeStruct->nPoints; i++ )
    {
        uint idx = this->isEvenTAvg? i : probeStruct->nPoints-1-i;

        // Compute the instantaneous spatial means of the velocity moments 
        real mean_vx = thrust::reduce(vx_iter_begin, vx_iter_end)/N;
        real mean_vy = thrust::reduce(vy_iter_begin, vy_iter_end)/N;
        real mean_vz = thrust::reduce(vz_iter_begin, vz_iter_end)/N;

        real vx2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
        real vy2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
        real vz2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;

        real vxvy = thrust::inner_product(vx_iter_begin, vx_iter_end, vy_iter_begin, 0.f)/N;
        real vxvz = thrust::inner_product(vx_iter_begin, vx_iter_end, vz_iter_begin, 0.f)/N;
        real vyvz = thrust::inner_product(vy_iter_begin, vy_iter_end, vz_iter_begin, 0.f)/N;

        real vx3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
        real vy3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
        real vz3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;

        real vx4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
        real vy4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
        real vz4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;

        
        // checkCudaErrors( cudaMemcpy(probeStruct->pointIndicesH, probeStruct->pointIndicesD, sizeof(int)*probeStruct->nIndices, cudaMemcpyDeviceToHost) );
        // if(i<2)
        {
            for( int j=0; j<9; j++ )
            {   //if((real)vx_iter_begin[j] != (real)vx_thrust[indices_thrust[j]])
                {
                printf("z: %u \t %f  \t N %f \n", idx, probeStruct->pointCoordsZ[idx], N);
                printf("IDX: \t %u  \t %u \n", probeStruct->pointIndicesH[j], (uint)indices_thrust[j]);  
                printf("vx[%u] permutation: \t %f \t %f \t %f \t %f \n\n", j, (real)vx_iter_begin[j], (real)vx_iter_begin[8], mean_vx, vx2);
                }
            }
        }

        if(i<probeStruct->nPoints-1)
        {
            vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, probeStruct->nIndices);
            if(this->isEvenTAvg) 
                moveIndicesInPosNormalDir<<<grid.grid, grid.threads>>>( probeStruct->pointIndicesD, probeStruct->nIndices, neighborNormal, para->getParD(level)->coordX_SP, para->getParD(level)->coordY_SP, para->getParD(level)->coordZ_SP );
            else 
                moveIndicesInNegNormalDir<<<grid.grid, grid.threads>>>( probeStruct->pointIndicesD, probeStruct->nIndices, para->getParD(level)->neighborWSB_SP, neighborInplane1, neighborInplane2, para->getParD(level)->coordX_SP, para->getParD(level)->coordY_SP, para->getParD(level)->coordZ_SP ); 
        } 
    }
    this->isEvenTAvg=!this->isEvenTAvg;
}
