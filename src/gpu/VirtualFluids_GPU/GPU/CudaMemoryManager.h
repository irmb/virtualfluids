#ifndef CudamemoryManager_H
#define CudamemoryManager_H

#include <vector>
#include <string>
#include <memory>
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <curand.h>
#include <curand_kernel.h>

class Parameter;
class PorousMedia;
class ActuatorLine;
class Probe;

class VIRTUALFLUIDS_GPU_EXPORT CudaMemoryManager
{
public:
	CudaMemoryManager(std::shared_ptr<Parameter> parameter);

	void setMemsizeGPU(double admem, bool reset);
	double getMemsizeGPU();

    void cudaAllocFull(int lev);
    void cudaFreeFull(int lev);

    void cudaCopyPrint(int lev);
    void cudaCopyMedianPrint(int lev);

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

    void cudaAllocF3SP(int lev);

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
	//3D domain decomposition
	void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv,
                                      int streamIndex);
    void cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend,
                                      int streamIndex);
	void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv,
                                      int streamIndex);
    void cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend,
                                      int streamIndex);
	void cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsRecv,
                                      int streamIndex);
    void cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor, const unsigned int &memsizeFsSend,
                                      int streamIndex);
	void cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor);

	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//3D domain decomposition F3
	void cudaAllocProcessNeighborF3X(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3XFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3XFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3XIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborF3X(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborF3Y(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3YFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3YFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3YIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborF3Y(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborF3Z(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3ZFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3ZFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborF3ZIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborF3Z(int lev, unsigned int processNeighbor);
	//////////////////////////////////////////////////////////////////////////

	void cudaAllocTurbulentViscosity(int lev);
    void cudaCopyTurbulentViscosityHD(int lev);
    void cudaCopyTurbulentViscosityDH(int lev);
    void cudaFreeTurbulentViscosity(int lev);

    void cudaAllocTurbulenceIntensity(int lev, uint size);
    void cudaCopyTurbulenceIntensityHD(int lev, uint size);
    void cudaCopyTurbulenceIntensityDH(int lev, uint size);
    void cudaFreeTurbulenceIntensity(int lev);

    void cudaAllocMedianSP(int lev);
    void cudaCopyMedianSP(int lev);
    void cudaFreeMedianSP(int lev);

    void cudaAllocMedianOut(int lev);
    void cudaFreeMedianOut(int lev);

    void cudaAllocMedianOutAD(int lev);
    void cudaFreeMedianOutAD(int lev);

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

    void cudaAllocWallModel(int lev, bool hasWallModelMonitor);
    void cudaCopyWallModel(int lev,  bool hasWallModelMonitor);
    void cudaFreeWallModel(int lev,  bool hasWallModelMonitor);

    void cudaAllocGeomValuesBC(int lev);
    void cudaCopyGeomValuesBC(int lev);
    void cudaFreeGeomValuesBC(int lev);

    void cudaAllocGeomNormals(int lev);
    void cudaCopyGeomNormals(int lev);
    void cudaFreeGeomNormals(int lev);

    void cudaAllocInflowNormals(int lev);
    void cudaCopyInflowNormals(int lev);
    void cudaFreeInflowNormals(int lev);

    void cudaAllocOutflowNormals(int lev);
    void cudaCopyOutflowNormals(int lev);
    void cudaFreeOutflowNormals(int lev);

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

    void cudaAllocInlet(int lev);
    void cudaCopyInlet(int lev);
    void cudaFreeInlet(int lev);
    void cudaAllocOutlet(int lev);
    void cudaCopyOutlet(int lev);
    void cudaFreeOutlet(int lev);


    void cudaAllocPressX0(int lev);
    void cudaCopyPressX0(int lev);
    void cudaFreePressX0(int lev);
    void cudaAllocPressX1(int lev);
    void cudaCopyPressX1(int lev);
    void cudaFreePressX1(int lev);

    void cudaAllocVeloPropeller(int lev);
    void cudaCopyVeloPropeller(int lev);
    void cudaFreeVeloPropeller(int lev);

    void cudaAllocMeasurePoints(int lev, int i);
    void cudaCopyMeasurePoints(int lev, int i);
    void cudaFreeMeasurePoints(int lev, int i);

    void cudaAllocMeasurePointsIndex(int lev);
    void cudaCopyMeasurePointsIndex(int lev);
    void cudaCopyMeasurePointsToHost(int lev);
    void cudaFreeMeasurePointsIndex(int lev);

    void cudaAllocFsForCheckPointAndRestart(int lev);
    void cudaCopyFsForRestart(int lev);
    void cudaCopyFsForCheckPoint(int lev);
    void cudaFreeFsForCheckPointAndRestart(int lev);

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


    void cudaAllocParticles(int lev);
    void cudaCopyParticles(int lev);
    void cudaFreeParticles(int lev);

    void cudaAllocRandomValues();

    void cudaAllocPorousMedia(PorousMedia* pm, int lev);
    void cudaCopyPorousMedia(PorousMedia* pm, int lev);
    void cudaFreePorousMedia(PorousMedia* pm, int lev);

    void cudaAllocConcentration(int lev);
    void cudaCopyConcentrationDeviceToHost(int lev);
    void cudaCopyConcentrationHostToDevice(int lev);
    void cudaFreeConcentration(int lev);

    void cudaAllocTempFs(int lev);

    void cudaAllocTempPressBC(int lev);
    void cudaCopyTempPressBCHD(int lev);
    void cudaFreeTempPressBC(int lev);

    void cudaAllocTempVeloBC(int lev);
    void cudaCopyTempVeloBCHD(int lev);
    void cudaFreeTempVeloBC(int lev);

    void cudaAllocTempNoSlipBC(int lev);
    void cudaCopyTempNoSlipBCHD(int lev);
    void cudaFreeTempNoSlipBC(int lev);

    void cudaAllocPlaneConcIn(int lev, int numofelem);
    void cudaCopyPlaneConcIn(int lev, int numofelem);
    void cudaAllocPlaneConcOut1(int lev, int numofelem);
    void cudaCopyPlaneConcOut1(int lev, int numofelem);
    void cudaAllocPlaneConcOut2(int lev, int numofelem);
    void cudaCopyPlaneConcOut2(int lev, int numofelem);
    void cudaFreePlaneConc(int lev);

    void cudaAllocProcessNeighbor(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborFsHD(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborFsDH(int lev, unsigned int processNeighbor);
    void cudaCopyProcessNeighborIndex(int lev, unsigned int processNeighbor);
    void cudaFreeProcessNeighbor(int lev, unsigned int processNeighbor);

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

    void cudaAllocFluidNodeIndices(int lev);
    void cudaCopyFluidNodeIndices(int lev);
    void cudaFreeFluidNodeIndices(int lev);
    void cudaAllocFluidNodeIndicesBorder(int lev);
    void cudaCopyFluidNodeIndicesBorder(int lev);
    void cudaFreeFluidNodeIndicesBorder(int lev);

    // Actuator Line
    void cudaAllocBladeRadii(ActuatorLine* actuatorLine);
    void cudaCopyBladeRadiiHtoD(ActuatorLine* actuatorLine);
    void cudaCopyBladeRadiiDtoH(ActuatorLine* actuatorLine);
    void cudaFreeBladeRadii(ActuatorLine* actuatorLine);

    void cudaAllocBladeCoords(ActuatorLine* actuatorLine);
    void cudaCopyBladeCoordsHtoD(ActuatorLine* actuatorLine);
    void cudaCopyBladeCoordsDtoH(ActuatorLine* actuatorLine);
    void cudaFreeBladeCoords(ActuatorLine* actuatorLine);

    void cudaAllocBladeIndices(ActuatorLine* actuatorLine);
    void cudaCopyBladeIndicesHtoD(ActuatorLine* actuatorLine);
    void cudaFreeBladeIndices(ActuatorLine* actuatorLine);

    void cudaAllocBladeVelocities(ActuatorLine* actuatorLine);
    void cudaCopyBladeVelocitiesHtoD(ActuatorLine* actuatorLine);
    void cudaCopyBladeVelocitiesDtoH(ActuatorLine* actuatorLine);
    void cudaFreeBladeVelocities(ActuatorLine* actuatorLine);

    void cudaAllocBladeForces(ActuatorLine* actuatorLine);
    void cudaCopyBladeForcesHtoD(ActuatorLine* actuatorLine);
    void cudaCopyBladeForcesDtoH(ActuatorLine* actuatorLine);
    void cudaFreeBladeForces(ActuatorLine* actuatorLine);

    void cudaAllocSphereIndices(ActuatorLine* actuatorLine);
    void cudaCopySphereIndicesHtoD(ActuatorLine* actuatorLine);
    void cudaFreeSphereIndices(ActuatorLine* actuatorLine);

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

private:
    std::shared_ptr<Parameter> parameter;
    double memsizeGPU = 0.;

};
#endif
