//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \ingroup gpu_core core
//! \{
#include "Probe.h"
#include "WallModelProbe.h"

#include <cuda_helper/CudaGrid.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Cuda/CudaMemoryManager.h"
#include "basics/constants/NumericConstants.h"

using namespace vf::basics::constant;
typedef thrust::device_vector<real>::iterator valIterator;
typedef thrust::device_vector<uint>::iterator indIterator;
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
            postProcessingVariables.push_back( PostProcessingVariable("dpdx_spatMean",     this->forceRatio) ); 
            postProcessingVariables.push_back( PostProcessingVariable("dpdy_spatMean",     this->forceRatio) );
            postProcessingVariables.push_back( PostProcessingVariable("dpdz_spatMean",     this->forceRatio) );
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
            postProcessingVariables.push_back( PostProcessingVariable("dpdx_spatTmpMean",     this->forceRatio) ); 
            postProcessingVariables.push_back( PostProcessingVariable("dpdy_spatTmpMean",     this->forceRatio) );
            postProcessingVariables.push_back( PostProcessingVariable("dpdz_spatTmpMean",     this->forceRatio) );
        }
        break;

    default:
        throw std::runtime_error("WallModelProbe::getPostProcessingVariables: Statistic unavailable!");
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
    if ( !para->getHasWallModelMonitor())                    throw std::runtime_error("WallModelProbe::findPoints(): !para->getHasWallModelMonitor() !");
    
    pointCoordsX_level.push_back(0); 
    pointCoordsY_level.push_back(0);
    pointCoordsZ_level.push_back(0);

    if(this->evaluatePressureGradient)
    {
        if (!para->getIsBodyForce()) throw std::runtime_error("WallModelProbe::findPoints(): bodyforce not allocated!");
        // Find all fluid nodes
        for(size_t pos = 1; pos < para->getParH(level)->numberOfNodes; pos++ )
        {
            if( para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID) 
            {
                probeIndices_level.push_back((int)pos);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

template<typename T>
T spatial_mean(T* device_pointer, uint numberOfPoints)
{
    thrust::device_ptr<T> thrust_pointer = thrust::device_pointer_cast(device_pointer);
    return thrust::reduce(thrust_pointer, thrust_pointer+numberOfPoints)/real(numberOfPoints);
}

template<typename T>
T index_based_spatial_mean(T* device_pointer, thrust::device_ptr<uint> indeces_ptr, uint numberOfIndeces)
{
    thrust::device_ptr<T> thrust_pointer = thrust::device_pointer_cast(device_pointer);

    thrust::permutation_iterator<valIterator, indIterator> iter_begin(thrust_pointer, indeces_ptr);
    thrust::permutation_iterator<valIterator, indIterator> iter_end  (thrust_pointer, indeces_ptr+numberOfIndeces);
    return thrust::reduce(iter_begin, iter_end)/real(numberOfIndeces);
}

template<typename T>
T compute_and_save_mean(T* device_pointer, uint numberOfPoints, T* quantitiesArray, uint timestep, uint numberOfTimesteps, uint indexOfArray)
{
    T mean = spatial_mean(device_pointer, numberOfPoints);
    quantitiesArray[calcArrayIndex(0, 1, timestep, numberOfTimesteps, indexOfArray)] = mean;
    return mean;
}

template<typename T>
T compute_and_save_index_based_mean(T* device_pointer, thrust::device_ptr<uint> indeces_ptr, uint numberOfIndices, T* quantitiesArray, uint timestep, uint numberOfTimesteps, uint indexOfArray)
{
    T mean = index_based_spatial_mean(device_pointer, indeces_ptr, numberOfIndices);
    quantitiesArray[calcArrayIndex(0, 1, timestep, numberOfTimesteps, indexOfArray)] = mean;
    return mean;
}

template<typename T>
void temporal_average(T* quantitiesArray, T currentValue, uint currentTimestep, uint numberOfTimesteps, uint oldTimestep, uint indexOfArray, real invNumberOfAverages)
{
    T oldMean = quantitiesArray[calcArrayIndex(0, 1, oldTimestep, numberOfTimesteps, indexOfArray)];
    quantitiesArray[calcArrayIndex(0, 1, currentTimestep, numberOfTimesteps, indexOfArray)] = oldMean + (currentValue-oldMean)*invNumberOfAverages;
}

void WallModelProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level)
{   
    bool doTmpAveraging = (t>this->getTStartTmpAveraging());
    uint numberOfStressBCPoints = para->getParD(level)->stressBC.numberOfBCnodes;
    if(numberOfStressBCPoints<1) return; //Skipping levels without StressBC
    uint timestep = probeStruct->timestepInTimeseries;
    real inv_n = c1o1/real(probeStruct->timestepInTimeAverage+1);
    uint oldTimestep = calcOldTimestep(timestep, probeStruct->lastTimestepInOldTimeseries);

    thrust::device_ptr<uint> indices_thrust = thrust::device_pointer_cast(probeStruct->pointIndicesD);

    if(probeStruct->quantitiesH[int(Statistic::SpatialMeans)])
    {
        uint arrOff = probeStruct->arrayOffsetsH[int(Statistic::SpatialMeans)];
        // Compute the instantaneous spatial means of the velocity moments 
        real spatMean_u_el      = compute_and_save_mean(para->getParD(level)->stressBC.Vx     , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+0);
        real spatMean_v_el      = compute_and_save_mean(para->getParD(level)->stressBC.Vy     , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+1);
        real spatMean_w_el      = compute_and_save_mean(para->getParD(level)->stressBC.Vz     , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+2);
        real spatMean_u1        = compute_and_save_mean(para->getParD(level)->stressBC.Vx1    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+3);
        real spatMean_v1        = compute_and_save_mean(para->getParD(level)->stressBC.Vy1    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+4);
        real spatMean_w1        = compute_and_save_mean(para->getParD(level)->stressBC.Vz1    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+5);
        real spatMean_u_star    = compute_and_save_mean(para->getParD(level)->wallModel.u_star, numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+6);
        real spatMean_Fx        = compute_and_save_mean(para->getParD(level)->wallModel.Fx    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+7);
        real spatMean_Fy        = compute_and_save_mean(para->getParD(level)->wallModel.Fy    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+8);
        real spatMean_Fz        = compute_and_save_mean(para->getParD(level)->wallModel.Fz    , numberOfStressBCPoints, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+9);

        real spatMean_dpdx;
        real spatMean_dpdy;
        real spatMean_dpdz;
        if(this->evaluatePressureGradient)
        {
            spatMean_dpdx = compute_and_save_index_based_mean(para->getParD(level)->forceX_SP, indices_thrust, probeStruct->nIndices, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+10);
            spatMean_dpdy = compute_and_save_index_based_mean(para->getParD(level)->forceY_SP, indices_thrust, probeStruct->nIndices, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+11);
            spatMean_dpdz = compute_and_save_index_based_mean(para->getParD(level)->forceZ_SP, indices_thrust, probeStruct->nIndices, probeStruct->quantitiesArrayH, timestep, probeStruct->nTimesteps, arrOff+12);
        }

        if(probeStruct->quantitiesH[int(Statistic::SpatioTemporalMeans)] && doTmpAveraging)
        {
            uint arrOff2 = probeStruct->arrayOffsetsH[int(Statistic::SpatioTemporalMeans)];
            temporal_average(probeStruct->quantitiesArrayH, spatMean_u_el  , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+0, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_v_el  , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+1, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_w_el  , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+2, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_u1    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+3, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_v1    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+4, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_w1    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+5, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_u_star, timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+6, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_Fx    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+7, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_Fy    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+8, inv_n);
            temporal_average(probeStruct->quantitiesArrayH, spatMean_Fz    , timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+9, inv_n);

            if(this->evaluatePressureGradient)
            {
                temporal_average(probeStruct->quantitiesArrayH, spatMean_dpdx, timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+10, inv_n);
                temporal_average(probeStruct->quantitiesArrayH, spatMean_dpdy, timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+11, inv_n);
                temporal_average(probeStruct->quantitiesArrayH, spatMean_dpdz, timestep, probeStruct->nTimesteps, oldTimestep, arrOff2+12, inv_n);
            }
        }    
    }
    getLastCudaError("WallModelProbe::calculateQuantities execution failed");
}

uint WallModelProbe::getNumberOfTimestepsInTimeseries(Parameter* para, int level)
{
    return this->tOut*exp2(level)/this->tAvg+1; 
}


//! \}
