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
//! time step they are called in. Most collision kernels (atm, all except K17CompressibleNavierStokes)
//! don't do this and would require an explicit call of calcMacroscopicQuantities. It does seem quite 
//! inexpensive though to simply save vx, vy, etc., directly in the collider.
//!
//! \todo might have to adapt conversionFactors when using grid refinement
//=======================================================================================

#ifndef Probe_H
#define Probe_H

#include <cuda.h>

#include <basics/PointerDefinitions.h>
#include "basics/writer/WbWriterVtkXmlBinary.h"

#include "gpu/core/PreCollisionInteractor/PreCollisionInteractor.h"

//=======================================================================================
//! \note How to add new Statistics 
//! Generally, the Statistic enum refers to the type of statistic to be calculated. 
//! It then depends on the derived probe class, which of these statistics are available. 
//! Some type of statistics are only suitable for a certain probe class, others might 
//! simply not have been implemented, yet.
//! For the same reasons it is also probe-specific, for which quantities (e.g. velocities, rho, etc.) these statistics are computed. 
//! The specific quantity (e.g., mean of vx, or variance of rho) is defined as PostProcessingVariable in getPostProcessingVariables of each respective probe.
//! PostProcessingVariable also holds the name and conversionFactor of the quantity that is required when writing the data to file
//! 
//! To add new Statistics:
//!     1. Add enum here, LAST has to stay last
//!     2. For PointProbe and PlaneProbe: add the computation of the statistic in switch statement in calculatePointwiseQuantities. 
//!     3. For PlanarAverageProbe and WallModelProbe: add the computation directly in calculateQuantities.
//!     4. In getPostProcessingVariables add the static in the switch statement and add the corresponding PostProcessingVariables
//!     5. Add Statistic to isAvailableStatistic of the respective probe
//!
//!  When adding new quantities to existing statistics (e.g., add rho to PlanarAverageProbe which currently only computes stats of velocity) only do steps 2 to 4
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
    std::function<real(int)> conversionFactor;
    PostProcessingVariable( std::string _name, 
                            std::function<real(int)>  _conversionFactor): 
    name(_name), conversionFactor(_conversionFactor){};
} PostProcessingVariable;

struct ProbeStruct{
    uint nPoints, nIndices, nArrays;
    uint nTimesteps=1;
    uint timestepInTimeseries=0;
    uint timestepInTimeAverage=0;
    uint lastTimestepInOldTimeseries=0;
    uint *pointIndicesH, *pointIndicesD;
    real *pointCoordsX, *pointCoordsY, *pointCoordsZ;
    bool hasDistances=false;
    real *distXH, *distYH, *distZH, *distXD, *distYD, *distZD;
    real *quantitiesArrayH, *quantitiesArrayD;
    bool *quantitiesH, *quantitiesD;
    uint *arrayOffsetsH, *arrayOffsetsD;
    bool isEvenTAvg = true;
};

__host__ __device__ int calcArrayIndex(int node, int nNodes, int timestep, int nTimesteps, int array);

__global__ void calcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint oldTimestepInTimeseries, uint timestepInTimeseries, uint timestepInAverage, uint nTimesteps,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );

__global__ void interpAndCalcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint oldTimestepInTimeseries, uint timestepInTimeseries, uint timestepInAverage, uint nTimesteps,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                );

uint calcOldTimestep(uint currentTimestep, uint lastTimestepInOldSeries);

class Probe : public PreCollisionInteractor 
{
public:
    Probe(
        const std::string probeName,
        const std::string outputPath,
        const uint tStartAvg,
        const uint tStartTmpAvg,
        const uint tAvg,
        const uint tStartOut,
        const uint tOut,
        const bool hasDeviceQuantityArray,
        const bool outputTimeSeries
    ):  probeName(probeName),
        outputPath(outputPath + (outputPath.back() == '/' ? "" : "/")),
        tStartAvg(tStartAvg),
        tStartTmpAveraging(tStartTmpAvg),
        tAvg(tAvg),
        tStartOut(tStartOut),
        tOut(tOut),
        hasDeviceQuantityArray(hasDeviceQuantityArray),
        outputTimeSeries(outputTimeSeries)
    {
        if (tStartOut < tStartAvg) throw std::runtime_error("Probe: tStartOut must be larger than tStartAvg!");
    }

    ~Probe();
    void init() override;
    void interact(int level, uint t) override;

    SPtr<ProbeStruct> getProbeStruct(int level){ return this->probeParams[level]; }

    void addStatistic(Statistic _variable);
    void addAllAvailableStatistics();
    
    bool getHasDeviceQuantityArray();
    uint getTStartTmpAveraging(){return this->tStartTmpAveraging;}

    void setFileNameToNOut(){this->fileNameLU = false;}

protected:
    virtual WbWriterVtkXmlBinary* getWriter(){ return WbWriterVtkXmlBinary::getInstance(); };
    real getNondimensionalConversionFactor(int level);

private:
    virtual bool isAvailableStatistic(Statistic _variable) = 0;

    virtual std::vector<PostProcessingVariable> getPostProcessingVariables(Statistic variable) = 0;

    virtual void findPoints(std::vector<int>& probeIndices_level,
                       std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                       std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                       int level) = 0;
    void addProbeStruct(std::vector<int>& probeIndices,
                        std::vector<real>& distX, std::vector<real>& distY, std::vector<real>& distZ,   
                        std::vector<real>& pointCoordsX, std::vector<real>& pointCoordsY, std::vector<real>& pointCoordsZ,
                        int level);
    virtual void calculateQuantities(SPtr<ProbeStruct> probeStruct, uint t, int level) = 0;

    virtual void write(int level, int t);
    virtual void writeParallelFile(int t);
    virtual void writeGridFile(int level, int t, uint part);
    std::string writeTimeseriesHeader(int level);
    void appendTimeseriesFile(int level, int t);

    std::vector<std::string> getVarNames();
    std::string makeGridFileName(int level, int id, int t, uint part);
    std::string makeParallelFileName(int id, int t);
    std::string makeTimeseriesFileName(int leve, int id);

    virtual uint getNumberOfTimestepsInTimeseries(int level){ (void)para; (void)level; return 1; }

protected:
    const std::string probeName;
    const std::string outputPath;

    std::vector<SPtr<ProbeStruct>> probeParams;
    bool quantities[int(Statistic::LAST)] = {};
    const bool hasDeviceQuantityArray;    //!> flag initiating memCopy in Point and PlaneProbe. Other probes are only based on thrust reduce functions and therefore dont need explict memCopy in interact()
    const bool outputTimeSeries;        //!> flag initiating time series output in Point and WallModelProbe.
    std::vector<std::string> fileNamesForCollectionFile;
    std::vector<std::string> timeseriesFileNames;

    bool fileNameLU = true; //!> if true, written file name contains time step in LU, else is the number of the written probe files

protected:
    const uint tStartAvg;
    const uint tStartTmpAveraging; //!> only non-zero in PlanarAverageProbe and WallModelProbe to switch on Spatio-temporal averaging (while only doing spatial averaging for t<tStartTmpAveraging) 
    const uint tAvg;  //! for tAvg==1 the probe will be evaluated in every sub-timestep of each respective level, else, the probe will only be evaluated in each synchronous time step 
    const uint tStartOut;
    const uint tOut;

    std::function<real(int)> velocityRatio;
    std::function<real(int)> densityRatio;
    std::function<real(int)> forceRatio;
    std::function<real(int)> stressRatio;
    std::function<real(int)> viscosityRatio;
    std::function<real(int)> nondimensional;
};

#endif