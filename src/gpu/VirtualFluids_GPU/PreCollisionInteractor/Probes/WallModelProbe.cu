#include "Probe.h"
#include "WallModelProbe.h"

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
bool WallModelProbe::isAvailableStatistic(Statistic _variable)
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
            isAvailable =  true;
            break;
        case Statistic::SpatialCovariances:
        case Statistic::SpatioTemporalCovariances:
        case Statistic::SpatialSkewness:
        case Statistic::SpatioTemporalSkewness:
        case Statistic::SpatialFlatness:
        case Statistic::SpatioTemporalFlatness:
            isAvailable =  false;
            break;
        default:
            isAvailable =  false;
    }
    return isAvailable;
}

///////////////////////////////////////////////////////////////////////////////////

std::vector<PostProcessingVariable> WallModelProbe::getPostProcessingVariables(Statistic statistic)
{
    std::vector<PostProcessingVariable> postProcessingVariables;
    switch (statistic)
    {
    case Statistic::SpatialMeans:
        postProcessingVariables.push_back( PostProcessingVariable("vx_el_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_el_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_el_spatMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vx1_spatMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy1_spatMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz1_spatMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("u_star_spatMean", this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("Fx_spatMean",     this->outputStress? this->stressRatio: this->forceRatio) ); 
        postProcessingVariables.push_back( PostProcessingVariable("Fy_spatMean",     this->outputStress? this->stressRatio: this->forceRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("Fz_spatMean",     this->outputStress? this->stressRatio: this->forceRatio) );
        if(this->evaluatePressureGradient)
        {
            postProcessingVariables.push_back( PostProcessingVariable("dpdx_spatMean",     this->accelerationRatio) ); 
            postProcessingVariables.push_back( PostProcessingVariable("dpdy_spatMean",     this->accelerationRatio) );
            postProcessingVariables.push_back( PostProcessingVariable("dpdz_spatMean",     this->accelerationRatio) );
        }
        break;
    case Statistic::SpatioTemporalMeans:
        postProcessingVariables.push_back( PostProcessingVariable("vx_el_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy_el_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz_el_spatTmpMean",  this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vx1_spatTmpMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vy1_spatTmpMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("vz1_spatTmpMean",    this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("u_star_spatTmpMean", this->velocityRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("Fx_spatTmpMean",     this->outputStress? this->stressRatio: this->forceRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("Fy_spatTmpMean",     this->outputStress? this->stressRatio: this->forceRatio) );
        postProcessingVariables.push_back( PostProcessingVariable("Fz_spatTmpMean",     this->outputStress? this->stressRatio: this->forceRatio) );
        if(this->evaluatePressureGradient)
        {
            postProcessingVariables.push_back( PostProcessingVariable("dpdx_spatTmpMean",     this->accelerationRatio) ); 
            postProcessingVariables.push_back( PostProcessingVariable("dpdy_spatTmpMean",     this->accelerationRatio) );
            postProcessingVariables.push_back( PostProcessingVariable("dpdz_spatTmpMean",     this->accelerationRatio) );
        }
        break;

    default:
        printf("Statistic unavailable in WallModelProbe\n");
        assert(false);
        break;
    }
    return postProcessingVariables;
}

///////////////////////////////////////////////////////////////////////////////////

void WallModelProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                            std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                            std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                            int level)
{
    assert( para->getParD(level)->stressBC.numberOfBCnodes > 0 && para->getHasWallModelMonitor() );

    real dt = para->getTimeRatio();
    uint nt = uint((para->getTEnd()-this->tStartAvg)/this->tAvg);
    
    for(uint t=0; t<nt; t++)
    {
        pointCoordsX_level.push_back(dt*(t*this->tAvg)+this->tStartAvg); // x coord will serve as time in this probe
        pointCoordsY_level.push_back(0);
        pointCoordsZ_level.push_back(0);
    }

    if(this->evaluatePressureGradient)
    {
        assert(para->getIsBodyForce());
        // Find all fluid nodes
        for(uint j=1; j<para->getParH(level)->numberOfNodes; j++ )
        {
            if( para->getParH(level)->typeOfGridNode[j] == GEO_FLUID) 
            {
                probeIndices_level.push_back(j);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

void WallModelProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level)
{   
    bool doTmpAveraging = (t>this->getTStartTmpAveraging());

    // Pointer casts to use device arrays in thrust reductions
    thrust::device_ptr<real> u_el_thrust    = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vx);
    thrust::device_ptr<real> v_el_thrust    = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vy);
    thrust::device_ptr<real> w_el_thrust    = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vz);
    thrust::device_ptr<real> u1_thrust      = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vx1);
    thrust::device_ptr<real> v1_thrust      = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vy1);
    thrust::device_ptr<real> w1_thrust      = thrust::device_pointer_cast(para->getParD(level)->stressBC.Vz1);
    thrust::device_ptr<real> u_star_thrust  = thrust::device_pointer_cast(para->getParD(level)->wallModel.u_star);
    thrust::device_ptr<real> Fx_thrust      = thrust::device_pointer_cast(para->getParD(level)->wallModel.Fx);
    thrust::device_ptr<real> Fy_thrust      = thrust::device_pointer_cast(para->getParD(level)->wallModel.Fy);
    thrust::device_ptr<real> Fz_thrust      = thrust::device_pointer_cast(para->getParD(level)->wallModel.Fz);
    thrust::device_ptr<real> dpdx_thrust    = thrust::device_pointer_cast(para->getParD(level)->forceX_SP);
    thrust::device_ptr<real> dpdy_thrust    = thrust::device_pointer_cast(para->getParD(level)->forceY_SP);
    thrust::device_ptr<real> dpdz_thrust    = thrust::device_pointer_cast(para->getParD(level)->forceZ_SP);

    thrust::device_ptr<uint> indices_thrust = thrust::device_pointer_cast(probeStruct->pointIndicesD);
    typedef thrust::device_vector<real>::iterator valIterator;
    typedef thrust::device_vector<uint>::iterator indIterator;
    thrust::permutation_iterator<valIterator, indIterator> dpdx_iter_begin(dpdx_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> dpdx_iter_end  (dpdx_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> dpdy_iter_begin(dpdy_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> dpdy_iter_end  (dpdy_thrust, indices_thrust+probeStruct->nIndices);
    thrust::permutation_iterator<valIterator, indIterator> dpdz_iter_begin(dpdz_thrust, indices_thrust);
    thrust::permutation_iterator<valIterator, indIterator> dpdz_iter_end  (dpdz_thrust, indices_thrust+probeStruct->nIndices);

    real N = para->getParD(level)->stressBC.numberOfBCnodes;
    real n = (real)probeStruct->vals;
    int nPoints = probeStruct->nPoints;

    if(probeStruct->quantitiesH[int(Statistic::SpatialMeans)])
    {
        // Compute the instantaneous spatial means of the velocity moments 
        real spatMean_u_el      = thrust::reduce(u_el_thrust, u_el_thrust+N)/N;
        real spatMean_v_el      = thrust::reduce(v_el_thrust, v_el_thrust+N)/N;
        real spatMean_w_el      = thrust::reduce(w_el_thrust, w_el_thrust+N)/N;
        real spatMean_u1        = thrust::reduce(u1_thrust, u1_thrust+N)/N;
        real spatMean_v1        = thrust::reduce(v1_thrust, v1_thrust+N)/N;
        real spatMean_w1        = thrust::reduce(w1_thrust, w1_thrust+N)/N;
        real spatMean_u_star    = thrust::reduce(u_star_thrust, u_star_thrust+N)/N;
        real spatMean_Fx        = thrust::reduce(Fx_thrust, Fx_thrust+N)/N;
        real spatMean_Fy        = thrust::reduce(Fy_thrust, Fy_thrust+N)/N;
        real spatMean_Fz        = thrust::reduce(Fz_thrust, Fz_thrust+N)/N;

        uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialMeans)];
        probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+tProbe] = spatMean_u_el;
        probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+tProbe] = spatMean_v_el;
        probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+tProbe] = spatMean_w_el;
        probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+tProbe] = spatMean_u1;
        probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+tProbe] = spatMean_v1;
        probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+tProbe] = spatMean_w1;
        probeStruct->quantitiesArrayH[(arrOff+6)*nPoints+tProbe] = spatMean_u_star;
        probeStruct->quantitiesArrayH[(arrOff+7)*nPoints+tProbe] = spatMean_Fx;
        probeStruct->quantitiesArrayH[(arrOff+8)*nPoints+tProbe] = spatMean_Fy;
        probeStruct->quantitiesArrayH[(arrOff+9)*nPoints+tProbe] = spatMean_Fz;

        real spatMean_dpdx;
        real spatMean_dpdy;
        real spatMean_dpdz;
        if(this->evaluatePressureGradient)
        {
            real N_fluid = (real)probeStruct->nIndices;
            spatMean_dpdx      = thrust::reduce(dpdx_iter_begin, dpdx_iter_end)/N_fluid;
            spatMean_dpdy      = thrust::reduce(dpdy_iter_begin, dpdy_iter_end)/N_fluid;
            spatMean_dpdz      = thrust::reduce(dpdz_iter_begin, dpdz_iter_end)/N_fluid;
            probeStruct->quantitiesArrayH[(arrOff+10)*nPoints+tProbe] = spatMean_dpdx;
            probeStruct->quantitiesArrayH[(arrOff+11)*nPoints+tProbe] = spatMean_dpdy;
            probeStruct->quantitiesArrayH[(arrOff+12)*nPoints+tProbe] = spatMean_dpdz;
        }

        if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalMeans)] && doTmpAveraging)
        {
            uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalMeans)];
            real spatMean_u_el_old      = probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+tProbe-1];
            real spatMean_v_el_old      = probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+tProbe-1];
            real spatMean_w_el_old      = probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+tProbe-1];
            real spatMean_u1_old        = probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+tProbe-1];
            real spatMean_v1_old        = probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+tProbe-1];
            real spatMean_w1_old        = probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+tProbe-1];
            real spatMean_u_star_old    = probeStruct->quantitiesArrayH[(arrOff+6)*nPoints+tProbe-1];
            real spatMean_Fx_old        = probeStruct->quantitiesArrayH[(arrOff+7)*nPoints+tProbe-1];
            real spatMean_Fy_old        = probeStruct->quantitiesArrayH[(arrOff+8)*nPoints+tProbe-1];
            real spatMean_Fz_old        = probeStruct->quantitiesArrayH[(arrOff+9)*nPoints+tProbe-1];

            probeStruct->quantitiesArrayH[(arrOff+0)*nPoints+tProbe] = spatMean_u_el_old + (spatMean_u_el-spatMean_u_el_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+1)*nPoints+tProbe] = spatMean_v_el_old + (spatMean_v_el-spatMean_v_el_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+2)*nPoints+tProbe] = spatMean_w_el_old + (spatMean_w_el-spatMean_w_el_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+3)*nPoints+tProbe] = spatMean_u1_old + (spatMean_u1-spatMean_u1_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+4)*nPoints+tProbe] = spatMean_v1_old + (spatMean_v1-spatMean_v1_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+5)*nPoints+tProbe] = spatMean_w1_old + (spatMean_w1-spatMean_w1_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+6)*nPoints+tProbe] = spatMean_u_star_old +(spatMean_u_star-spatMean_u_star_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+7)*nPoints+tProbe] = spatMean_Fx_old + (spatMean_Fx-spatMean_Fx_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+8)*nPoints+tProbe] = spatMean_Fy_old + (spatMean_Fy-spatMean_Fy_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+9)*nPoints+tProbe] = spatMean_Fz_old + (spatMean_Fz-spatMean_Fz_old)/n;

            if(this->evaluatePressureGradient)
            {
            real spatMean_dpdx_old     = probeStruct->quantitiesArrayH[(arrOff+10)*nPoints+tProbe-1];
            real spatMean_dpdy_old     = probeStruct->quantitiesArrayH[(arrOff+11)*nPoints+tProbe-1];
            real spatMean_dpdz_old     = probeStruct->quantitiesArrayH[(arrOff+12)*nPoints+tProbe-1];
            probeStruct->quantitiesArrayH[(arrOff+10)*nPoints+tProbe] = spatMean_dpdx_old + (spatMean_dpdx-spatMean_dpdx_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+11)*nPoints+tProbe] = spatMean_dpdy_old + (spatMean_dpdy-spatMean_dpdy_old)/n;
            probeStruct->quantitiesArrayH[(arrOff+12)*nPoints+tProbe] = spatMean_dpdz_old + (spatMean_dpdz-spatMean_dpdz_old)/n;
            }
        }    
    }

    this->tProbe += 1;
    getLastCudaError("WallModelProbe::calculateQuantities execution failed");
}

