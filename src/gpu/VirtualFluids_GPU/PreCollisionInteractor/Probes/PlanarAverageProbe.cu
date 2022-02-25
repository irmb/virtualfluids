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
    bool isAvailable;

    switch (_variable)
    {
        case PostProcessingVariable::Instantaneous:
        case PostProcessingVariable::Means:
        case PostProcessingVariable::Variances:
            isAvailable = false;
            break;
        case PostProcessingVariable::SpatialMeans:
        case PostProcessingVariable::SpatioTemporalMeans:
        case PostProcessingVariable::SpatialCovariances:
        case PostProcessingVariable::SpatioTemporalCovariances:
        case PostProcessingVariable::SpatialSkewness:
        case PostProcessingVariable::SpatioTemporalSkewness:
        case PostProcessingVariable::SpatialFlatness:
        case PostProcessingVariable::SpatioTemporalFlatness:
            isAvailable =  true;
            break;
        default:
            isAvailable =  false;
    }
    return isAvailable;
}

///////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////

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
    real n = (real)probeStruct->vals;
    int nPoints = probeStruct->nPoints;
    // Permutation iterators for direct iteration over the velocities of the planes
    typedef thrust::device_vector<real>::iterator valIterator;
    typedef thrust::device_vector<uint>::iterator indIterator;
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_begin(vx_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vx_iter_end  (vx_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_begin(vy_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vy_iter_end  (vy_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_begin(vz_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> vz_iter_end  (vz_thrust, indices_thrust+probeStruct->nIndices);

    for( int i=0; i<nPoints; i++ )
    {
        uint node = this->isEvenTAvg? i : nPoints-1-i; // Note, loop moves in positive normal dir at even calls and in negative normal dir in odd calls

        if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatialMeans)])
        {
            // Compute the instantaneous spatial means of the velocity moments 
            real spatMean_vx = thrust::reduce(vx_iter_begin, vx_iter_end)/N;
            real spatMean_vy = thrust::reduce(vy_iter_begin, vy_iter_end)/N;
            real spatMean_vz = thrust::reduce(vz_iter_begin, vz_iter_end)/N;

            uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatialMeans)];
            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_vx;
            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_vy;
            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_vz;

            if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatioTemporalMeans)])
            {
            uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatioTemporalMeans)];
            real spatTmpMean_vx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
            real spatTmpMean_vy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
            real spatTmpMean_vz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];

            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_vx-spatTmpMean_vx_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_vy-spatTmpMean_vy_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_vz-spatTmpMean_vz_old)/n;
            }
        
            if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatialCovariances)])
            {   // <u_i' u_j'> = <u_i u_j> - <u_i>*<u_i> 
                real vx2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vy2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vz2 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow2<real>(), 0.f, thrust::plus<real>())/N;
                real vxvy = thrust::inner_product(vx_iter_begin, vx_iter_end, vy_iter_begin, 0.f)/N;
                real vxvz = thrust::inner_product(vx_iter_begin, vx_iter_end, vz_iter_begin, 0.f)/N;
                real vyvz = thrust::inner_product(vy_iter_begin, vy_iter_end, vz_iter_begin, 0.f)/N;
                real spatMean_vxvx = vx2-spatMean_vx*spatMean_vx;
                real spatMean_vyvy = vy2-spatMean_vy*spatMean_vy;
                real spatMean_vzvz = vz2-spatMean_vz*spatMean_vz;
                real spatMean_vxvy = vxvy-spatMean_vx*spatMean_vy;
                real spatMean_vxvz = vxvy-spatMean_vx*spatMean_vz;
                real spatMean_vyvz = vxvy-spatMean_vy*spatMean_vz;

                uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatialCovariances)];
                probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_vxvx;
                probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_vyvy;
                probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_vzvz;
                probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] = spatMean_vxvy;
                probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node] = spatMean_vxvz;
                probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node] = spatMean_vyvz;

                if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatioTemporalCovariances)])
                {
                    uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatioTemporalCovariances)];
                    real spatTmpMean_vxvx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                    real spatTmpMean_vyvy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                    real spatTmpMean_vzvz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];
                    real spatTmpMean_vxvy_old = probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node];
                    real spatTmpMean_vxvz_old = probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node];
                    real spatTmpMean_vyvz_old = probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node];

                    probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_vxvx-spatTmpMean_vxvx_old)/n;
                    probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_vyvy-spatTmpMean_vyvy_old)/n;
                    probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_vzvz-spatTmpMean_vzvz_old)/n;
                    probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+node] += (spatMean_vxvy-spatTmpMean_vxvy_old)/n;
                    probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+node] += (spatMean_vxvz-spatTmpMean_vxvz_old)/n;
                    probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+node] += (spatMean_vyvz-spatTmpMean_vyvz_old)/n;
                }

                if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatialSkewness)])
                {   // <u_i'^3> = <u_i^3> - <u_i>^3 - 3 <u_i> <u_i'^2>
                    real vx3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    real vy3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    real vz3 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow3<real>(), 0.f, thrust::plus<real>())/N;
                    real spatMean_vxvxvx = vx3 - spatMean_vx*spatMean_vx*spatMean_vx - 3*spatMean_vx*spatMean_vxvx;
                    real spatMean_vyvyvy = vy3 - spatMean_vy*spatMean_vy*spatMean_vy - 3*spatMean_vy*spatMean_vyvy;
                    real spatMean_vzvzvz = vz3 - spatMean_vz*spatMean_vz*spatMean_vz - 3*spatMean_vz*spatMean_vzvz;
                    real spatMean_Sx = ( vx3 - 3*spatMean_vx*spatMean_vxvx - spatMean_vx*spatMean_vx*spatMean_vx )/pow(spatMean_vxvx, 1.5f);
                    real spatMean_Sy = ( vy3 - 3*spatMean_vy*spatMean_vyvy - spatMean_vy*spatMean_vy*spatMean_vy )/pow(spatMean_vyvy, 1.5f);
                    real spatMean_Sz = ( vz3 - 3*spatMean_vz*spatMean_vzvz - spatMean_vz*spatMean_vz*spatMean_vz )/pow(spatMean_vzvz, 1.5f);

                    uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatialSkewness)];
                    probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_Sx;
                    probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_Sy;
                    probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_Sz;

                    if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatioTemporalSkewness)])
                    {
                        uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatioTemporalSkewness)];
                        real spatTmpMean_Sx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                        real spatTmpMean_Sy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                        real spatTmpMean_Sz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];

                        probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_Sx-spatTmpMean_Sx_old)/n;
                        probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_Sy-spatTmpMean_Sy_old)/n;
                        probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_Sz-spatTmpMean_Sz_old)/n;
                    }

                    if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatialFlatness)])
                    {   // <u_i'^4> = <u_i^4> - <u_i>^4 - 6 <u_i>^2 <u_i'^2> - 4 <u> <u'^3>
                        real vx4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        real vy4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        real vz4 = thrust::transform_reduce(vx_iter_begin, vx_iter_end, pow4<real>(), 0.f, thrust::plus<real>())/N;
                        real spatMean_vxvxvxvx = vx4 - spatMean_vx*spatMean_vx*spatMean_vx*spatMean_vx - 6*spatMean_vx*spatMean_vx*spatMean_vxvx - 4*spatMean_vx*spatMean_vxvxvx;
                        real spatMean_vyvyvyvy = vy4 - spatMean_vy*spatMean_vy*spatMean_vy*spatMean_vy - 6*spatMean_vy*spatMean_vx*spatMean_vyvy - 4*spatMean_vy*spatMean_vyvyvy;
                        real spatMean_vzvzvzvz = vz4 - spatMean_vz*spatMean_vz*spatMean_vz*spatMean_vz - 6*spatMean_vz*spatMean_vx*spatMean_vzvz - 4*spatMean_vz*spatMean_vzvzvz;
                        real spatMean_Fx = ( vx4 - 3*spatMean_vx*spatMean_vxvx - spatMean_vx*spatMean_vx*spatMean_vx*spatMean_vx )/(spatMean_vxvx*spatMean_vxvx);
                        real spatMean_Fy = ( vy4 - 3*spatMean_vy*spatMean_vyvy - spatMean_vx*spatMean_vy*spatMean_vy*spatMean_vy )/(spatMean_vyvy*spatMean_vyvy);
                        real spatMean_Fz = ( vz4 - 3*spatMean_vz*spatMean_vzvz - spatMean_vx*spatMean_vz*spatMean_vz*spatMean_vz )/(spatMean_vzvz*spatMean_vzvz);

                        uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatialFlatness)];
                        probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] = spatMean_Fx;
                        probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] = spatMean_Fy;
                        probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] = spatMean_Fz;

                        if(probeStruct->quantitiesH[int(PostProcessingVariable::SpatioTemporalFlatness)])
                        {
                            uint arrOff = probeStruct->arrayOffsetsH[int(PostProcessingVariable::SpatioTemporalFlatness)];
                            real spatTmpMean_Fx_old = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node];
                            real spatTmpMean_Fy_old = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node];
                            real spatTmpMean_Fz_old = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node];

                            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+node] += (spatMean_Fx-spatTmpMean_Fx_old)/n;
                            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+node] += (spatMean_Fy-spatTmpMean_Fy_old)/n;
                            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+node] += (spatMean_Fz-spatTmpMean_Fz_old)/n;
                        }
                    }
                }
        }
    }
        // checkCudaErrors( cudaMemcpy(probeStruct->pointIndicesH, probeStruct->pointIndicesD, sizeof(int)*probeStruct->nIndices, cudaMemcpyDeviceToHost) );
        // if(i<2)
        // {
        //     for( int j=0; j<9; j++ )
        //     {   //if((real)vx_iter_begin[j] != (real)vx_thrust[indices_thrust[j]])
        //         {
        //         printf("z: %u \t %f  \t N %f \n", idx, probeStruct->pointCoordsZ[idx], N);
        //         printf("IDX: \t %u  \t %u \n", probeStruct->pointIndicesH[j], (uint)indices_thrust[j]);  
        //         printf("vx[%u] permutation: \t %f \t %f \t %f \t %f \n\n", j, (real)vx_iter_begin[j], (real)vx_iter_begin[8], mean_vx, vx2);
        //         }
        //     }
        // }

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
