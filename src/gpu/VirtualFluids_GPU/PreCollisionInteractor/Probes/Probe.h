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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Probe.h
//! \author Henry Korb, Henrik Asmuth
//! \date 13/05/2022
//! \brief Base class for probes called in UpdateGrid27
//!
//! Any probe should be initiated in the app and added via para->addProbe( someProbe )
//! Note, that all probes generally require that macroscopic variables have been updated in the 
//! time step they are called. Most collision kernels (atm, all except TurbulentViscosityCumulantK17CompChim )
//! don't do this and would require an explicit call of calcMacroscopicQuantities. It does seem quite 
//! inexpensive though to simply save vx, vy, etc., directly in the collider.
//!
//=======================================================================================

#ifndef Probe_H
#define Probe_H

#include <cuda.h>

#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "PointerDefinitions.h"

//=======================================================================================
//! \note How to add new Statistics 
//! Generally, the Statistic enum refers to the type of statistic to be calculated. 
//! It then depends on the derived probe class, which of these statistics are available. 
//! Some type of statistics are simply only suitable for a certain class, others might 
//! simply not have been implemented, yet.
//! For the same reasons it is also probe-specific, for which quantities (e.g. velocities, rho, etc.) these statistics are computed. 
//! The specific quantity (e.g., mean of vx, or variance of rho) is defined as PostProcessingVariable in getPostProcessingVariables of each respective probe.
//! PostProcessingVariable also holds the name and conversionFactor of the quantity that is required when writing the data to file
//! 
//! To add new Statistics:
//!     1. Add enum here, LAST has to stay last
//!     2. For PointProbe and PlaneProbe: add the computation of the statistic in switch statement in calculatePointwiseQuantities. 
//!     3. For PlanarAverageProbe and WallModelProbe: add the computation directly in calculateQuantities.
//!     4. In getPostProcessingVariables add the static in the switch statement and add the corresponding PostProcessingVariable
//!     5. Add Statistic to isAvailableStatistic of the respective probe
//!     6. When adding new quantities to existing statistics (e.g., add rho to PlanarAverageProbe which currently only computes velocity stats) only do steps 2 to 4
//!

enum class Statistic{ 
    // Variables currently available in Point and Plane probe (all temporal pointwise statistics)
    Instantaneous,
    Means,
    Variances,

    // Variables available in PlanarAverage probe and (partially) in WallModelProbe
    // Spatial statistics are typically computed across fixed spatial subdomains, e.g. a plane of constant height
    // Spatio-temporal statistics additionally average the spatial stats in time
    SpatialMeans,
    SpatioTemporalMeans,
    SpatialCovariances,
    SpatioTemporalCovariances,
    SpatialSkewness,
    SpatioTemporalSkewness,
    SpatialFlatness,
    SpatioTemporalFlatness,
    LAST,
};

typedef struct PostProcessingVariable{
    std::string name;
    real conversionFactor;
    PostProcessingVariable( std::string _name, 
                            real        _conversionFactor): 
    name(_name), conversionFactor(_conversionFactor){};
} PostProcessingVariable;

struct ProbeStruct{
    uint nPoints, nIndices, nArrays, vals;
    uint *pointIndicesH, *pointIndicesD;
    real *pointCoordsX, *pointCoordsY, *pointCoordsZ;
    bool hasDistances=false;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    real *quantitiesArrayH, *quantitiesArrayD;
    bool *quantitiesH, *quantitiesD;
    uint *arrayOffsetsH, *arrayOffsetsD;
};

__global__ void calcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );

__global__ void interpAndCalcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );


class Probe : public PreCollisionInteractor 
{
public:
    Probe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tStartTmpAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut,
        bool _hasDeviceQuantityArray
    ):  probeName(_probeName),
        outputPath(_outputPath),
        tStartAvg(_tStartAvg),
        tAvg(_tAvg),
        tStartOut(_tStartOut),
        tOut(_tOut),
        hasDeviceQuantityArray(_hasDeviceQuantityArray),
        tStartTmpAveraging(_tStartTmpAvg),
        PreCollisionInteractor()
    {
        assert("Output starts before averaging!" && tStartOut>=tStartAvg);
    }
    
    void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager) override;
    void interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t) override;
    void free(Parameter* para, CudaMemoryManager* cudaManager) override;

    SPtr<ProbeStruct> getProbeStruct(int level){ return this->probeParams[level]; }

    void addStatistic(Statistic _variable);
    void addAllAvailableStatistics();
    
    bool getHasDeviceQuantityArray();
    uint getTStartTmpAveraging(){return this->tStartTmpAveraging;}

    void setFileNameToNOut(){this->fileNameLU = false;}
    void setTStartTmpAveraging(uint _tStartTmpAveraging){this->tStartTmpAveraging = _tStartTmpAveraging;}

private:
    virtual bool isAvailableStatistic(Statistic _variable) = 0;

    virtual std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) = 0;

    virtual void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                       std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                       std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                       int level) = 0;
    void addProbeStruct(CudaMemoryManager* cudaManager, std::vector<int>& probeIndices,
                        std::vector<real>& distX, std::vector<real>& distY, std::vector<real>& distZ,   
                        std::vector<real>& pointCoordsX, std::vector<real>& pointCoordsY, std::vector<real>& pointCoordsZ,
                        int level);
    virtual void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, uint t, int level) = 0;

    virtual void write(Parameter* para, int level, int t);
    void writeCollectionFile(Parameter* para, int t);
    void writeGridFiles(Parameter* para, int level, std::vector<std::string >& fnames, int t);
    std::vector<std::string> getVarNames();
    
private:
    const std::string probeName;
    const std::string outputPath;

    std::vector<SPtr<ProbeStruct>> probeParams;
    bool quantities[int(Statistic::LAST)] = {};
    bool hasDeviceQuantityArray; 
    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> varNames;

    bool fileNameLU = true; //!> if true, written file name contains time step in LU, else is the number of the written probe files

protected:
    uint tStartAvg;
    uint tStartTmpAveraging; //!> only non-zero in PlanarAverageProbe to switch on Spatio-temporal averaging (while only doing spatial averaging for t<tStartTmpAveraging) 
    uint tAvg;
    uint tStartOut;
    uint tOut;

    real velocityRatio;
    real densityRatio;
    real forceRatio;
};

#endif