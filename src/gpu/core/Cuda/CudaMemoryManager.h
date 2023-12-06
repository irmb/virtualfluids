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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef CudamemoryManager_H
#define CudamemoryManager_H

#include <memory>
#include <string>
#include <vector>

#include "Calculation/Calculation.h"

class Parameter;
class ActuatorFarm;
class Probe;
class VelocitySetter;
class PrecursorWriter;

class CudaMemoryManager
{
public:
    CudaMemoryManager(std::shared_ptr<Parameter> parameter);
    virtual ~CudaMemoryManager() = default;

    void setMemsizeGPU(double admem, bool reset);
    double getMemsizeGPU();

    void cudaCopyPrint(int lev);
    void cudaCopyMeanPrint(int lev);

    void cudaAllocCoord(int lev);
    void cudaCopyCoord(int lev);
    void cudaFreeCoord(int lev);

    void cudaAllocBodyForce(int lev);
    void cudaCopyBodyForce(int lev);
    void cudaFreeBodyForce(int lev);

    void cudaCopyDataToHost(int lev);

    void cudaAllocSP(int lev);
    void cudaCopySP(int lev);
    void cudaFreeSP(int lev);

    void cudaAllocNeighborWSB(int lev);
    void cudaCopyNeighborWSB(int lev);
    void cudaFreeNeighborWSB(int lev);

    void cudaAllocVeloBC(int lev);
    void cudaCopyVeloBC(int lev);
    void cudaFreeVeloBC(int lev);

    void cudaAllocOutflowBC(int lev);
    void cudaCopyOutflowBC(int lev);
    void cudaFreeOutflowBC(int lev);

    void cudaAllocNoSlipBC(int lev);
    void cudaCopyNoSlipBC(int lev);
    void cudaFreeNoSlipBC(int lev);

    void cudaAllocGeomBC(int lev);
    void cudaCopyGeomBC(int lev);
    void cudaFreeGeomBC(int lev);

    void cudaAllocPress(int lev);
    void cudaCopyPress(int lev);
    void cudaFreePress(int lev);

    void cudaAllocForcing();
    void cudaCopyForcingToDevice();
    void cudaCopyForcingToHost();
    void cudaFreeForcing();

    void cudaAllocLevelForcing(int level);
    void cudaCopyLevelForcingToDevice(int level);
    void cudaFreeLevelForcing(int level);

    void cudaAllocQuadricLimiters();
    void cudaCopyQuadricLimitersToDevice();
    void cudaFreeQuadricLimiters();

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    virtual void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv);
    void cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend);
    virtual void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor);
    //
    virtual void cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv);
    void cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend);
    virtual void cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor);
    //
    virtual void cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv);
    void cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend);
    virtual void cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor);

    //////////////////////////////////////////////////////////////////////////

    void cudaAllocTurbulentViscosity(int lev);
    void cudaCopyTurbulentViscosityHD(int lev);
    void cudaCopyTurbulentViscosityDH(int lev);
    void cudaFreeTurbulentViscosity(int lev);

    void cudaAllocTurbulenceIntensity(int lev, uint size);
    void cudaCopyTurbulenceIntensityHD(int lev, uint size);
    void cudaCopyTurbulenceIntensityDH(int lev, uint size);
    void cudaFreeTurbulenceIntensity(int lev);

    void cudaAllocMeanSP(int lev);
    void cudaCopyMeanSP(int lev);
    void cudaFreeMeanSP(int lev);

    void cudaAllocMeanOut(int lev);
    void cudaFreeMeanOut(int lev);

    void cudaAllocMeanOutAD(int lev);
    void cudaFreeMeanOutAD(int lev);

    void cudaAllocInterfaceCF(int lev);
    void cudaCopyInterfaceCF(int lev);
    void cudaFreeInterfaceCF(int lev);

    void cudaAllocInterfaceFC(int lev);
    void cudaCopyInterfaceFC(int lev);
    void cudaCheckInterfaceFCBulk(int lev);
    void cudaFreeInterfaceFC(int lev);

    void cudaAllocInterfaceOffCF(int lev);
    void cudaCopyInterfaceOffCF(int lev);
    void cudaFreeInterfaceOffCF(int lev);

    void cudaAllocInterfaceOffFC(int lev);
    void cudaCopyInterfaceOffFC(int lev);
    void cudaFreeInterfaceOffFC(int lev);

    void cudaAllocSlipBC(int lev);
    void cudaCopySlipBC(int lev);
    void cudaFreeSlipBC(int lev);

    void cudaAllocStressBC(int lev);
    void cudaCopyStressBC(int lev);
    void cudaFreeStressBC(int lev);

    void cudaAllocPrecursorBC(int lev);
    void cudaAllocPrecursorData(int lev);
    void cudaCopyPrecursorBC(int lev);
    void cudaCopyPrecursorData(int lev);
    void cudaFreePrecursorBC(int lev);
    void cudaFreePrecursorData(int lev);

    void cudaAllocWallModel(int lev, bool hasWallModelMonitor);
    void cudaCopyWallModel(int lev, bool hasWallModelMonitor);
    void cudaFreeWallModel(int lev, bool hasWallModelMonitor);

    void cudaAllocGeomValuesBC(int lev);
    void cudaCopyGeomValuesBC(int lev);
    void cudaFreeGeomValuesBC(int lev);

    void cudaAllocTestRE(int lev, unsigned int size);
    void cudaCopyTestREtoDevice(int lev, unsigned int size);
    void cudaCopyTestREtoHost(int lev, unsigned int size);
    void cudaFreeTestRE(int lev);

    void cudaAllocCpTop(int lev);
    void cudaCopyCpTopInit(int lev);
    void cudaCopyCpTop(int lev);
    void cudaFreeCpTop(int lev);

    void cudaAllocCpBottom(int lev);
    void cudaCopyCpBottomInit(int lev);
    void cudaCopyCpBottom(int lev);
    void cudaFreeCpBottom(int lev);

    void cudaAllocCpBottom2(int lev);
    void cudaCopyCpBottom2Init(int lev);
    void cudaCopyCpBottom2(int lev);
    void cudaFreeCpBottom2(int lev);

    void cudaAllocConcFile(int lev);
    void cudaCopyConcFile(int lev);
    void cudaFreeConcFile(int lev);

    void cudaAllocMeasurePointsIndex(int lev);
    void cudaCopyMeasurePointsIndex(int lev);
    void cudaCopyMeasurePointsToHost(int lev);
    void cudaFreeMeasurePointsIndex(int lev);

    void cudaAllocFsForCheckPointAndRestart(int lev) const;
    void cudaAllocFsForAllLevelsOnHost() const;
    //! \brief copy distributions from host to device
    void cudaCopyFsForRestart(int lev) const;
    //! \brief copy distributions from device to host
    void cudaCopyFsForCheckPoint(int lev) const;
    void cudaCopyFsForAllLevelsToHost() const;
    void cudaFreeFsForCheckPointAndRestart(int lev) const;

    void cudaAllocDragLift(int lev, int numofelem);
    void cudaCopyDragLift(int lev, int numofelem);
    void cudaFreeDragLift(int lev);

    void cudaAlloc2ndMoments(int lev, int numofelem);
    void cudaCopy2ndMoments(int lev, int numofelem);
    void cudaFree2ndMoments(int lev);

    void cudaAlloc3rdMoments(int lev, int numofelem);
    void cudaCopy3rdMoments(int lev, int numofelem);
    void cudaFree3rdMoments(int lev);

    void cudaAllocHigherMoments(int lev, int numofelem);
    void cudaCopyHigherMoments(int lev, int numofelem);
    void cudaFreeHigherMoments(int lev);

    void cudaAllocForceVelo(int lev, int numofelem);
    void cudaCopyForceVelo(int lev, int numofelem);
    void cudaFreeForceVelo(int lev);

    void cudaAlloc2ndOrderDerivitivesIsoTest(int lev);
    void cudaCopy2ndOrderDerivitivesIsoTestDH(int lev);
    void cudaCopy2ndOrderDerivitivesIsoTestHD(int lev);
    void cudaFree2ndOrderDerivitivesIsoTest(int lev);

    void cudaAllocConcentration(int lev);
    void cudaCopyConcentrationDeviceToHost(int lev);
    void cudaCopyConcentrationHostToDevice(int lev);
    void cudaFreeConcentration(int lev);

    void cudaAllocConcentrationFs(int lev);

    void cudaAllocConcentrationDirichletBC(int lev);
    void cudaCopyConcentrationDirichletBCHostToDevice(int lev);
    void cudaFreeConcentrationDirichletBC(int lev);

    void cudaAllocConcentrationNoSlipBC(int lev);
    void cudaCopyConcentrationNoSlipBCHD(int lev);
    void cudaFreeConcentrationNoSlipBC(int lev);

    void cudaAllocProcessNeighborADX(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADXFsHD(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADXFsDH(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADXIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborADX(int lev, unsigned int processNeighbor);

    void cudaAllocProcessNeighborADY(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADYFsHD(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADYFsDH(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADYIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborADY(int lev, unsigned int processNeighbor);

    void cudaAllocProcessNeighborADZ(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADZFsHD(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADZFsDH(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborADZIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborADZ(int lev, unsigned int processNeighbor);

    void cudaAllocTaggedFluidNodeIndices(CollisionTemplate tag, int lev);
    void cudaCopyTaggedFluidNodeIndices(CollisionTemplate tag, int lev);
    void cudaFreeTaggedFluidNodeIndices(CollisionTemplate tag, int lev);

    // ActuatorFarm
    void cudaAllocBladeGeometries(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeGeometriesHtoD(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeGeometriesDtoH(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeGeometries(ActuatorFarm* actuatorFarm);

    void cudaAllocBladeOrientations(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeOrientationsHtoD(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeOrientationsDtoH(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeOrientations(ActuatorFarm* actuatorFarm);

    void cudaAllocBladeCoords(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeCoordsHtoD(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeCoordsDtoH(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeCoords(ActuatorFarm* actuatorFarm);

    void cudaAllocBladeIndices(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeIndicesHtoD(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeIndices(ActuatorFarm* actuatorFarm);

    void cudaAllocBladeVelocities(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeVelocitiesHtoD(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeVelocitiesDtoH(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeVelocities(ActuatorFarm* actuatorFarm);

    void cudaAllocBladeForces(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeForcesHtoD(ActuatorFarm* actuatorFarm);
    void cudaCopyBladeForcesDtoH(ActuatorFarm* actuatorFarm);
    void cudaFreeBladeForces(ActuatorFarm* actuatorFarm);

    void cudaAllocSphereIndices(ActuatorFarm* actuatorFarm);
    void cudaCopySphereIndicesHtoD(ActuatorFarm* actuatorFarm);
    void cudaFreeSphereIndices(ActuatorFarm* actuatorFarm);
    // Probes
    void cudaAllocProbeDistances(Probe* probe, int level);
    void cudaCopyProbeDistancesHtoD(Probe* probe, int level);
    void cudaCopyProbeDistancesDtoH(Probe* probe, int level);
    void cudaFreeProbeDistances(Probe* probe, int level);

    void cudaAllocProbeIndices(Probe* probe, int level);
    void cudaCopyProbeIndicesHtoD(Probe* probe, int level);
    void cudaCopyProbeIndicesDtoH(Probe* probe, int level);
    void cudaFreeProbeIndices(Probe* probe, int level);

    void cudaAllocProbeQuantityArray(Probe* probe, int level);
    void cudaCopyProbeQuantityArrayHtoD(Probe* probe, int level);
    void cudaCopyProbeQuantityArrayDtoH(Probe* probe, int level);
    void cudaFreeProbeQuantityArray(Probe* probe, int level);

    void cudaAllocProbeQuantitiesAndOffsets(Probe* probe, int level);
    void cudaCopyProbeQuantitiesAndOffsetsHtoD(Probe* probe, int level);
    void cudaCopyProbeQuantitiesAndOffsetsDtoH(Probe* probe, int level);
    void cudaFreeProbeQuantitiesAndOffsets(Probe* probe, int level);

    // Precursor Writer
    void cudaAllocPrecursorWriter(PrecursorWriter* writer, int level);
    void cudaCopyPrecursorWriterIndicesHtoD(PrecursorWriter* writer, int level);
    void cudaCopyPrecursorWriterOutputVariablesDtoH(PrecursorWriter* writer, int level);
    void cudaFreePrecursorWriter(PrecursorWriter* writer, int level);

private:
    std::shared_ptr<Parameter> parameter;
    double memsizeGPU = 0.0;
};

#endif
