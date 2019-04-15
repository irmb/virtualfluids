#include "TrafficTimestep.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <sstream>

#include <thrust/reduce.h>

#include <iomanip>

#include "RoadNetwork/RoadNetworkData.h"
#include "Junction/Junction.h"
#include "Sink/Sink.h"
#include "Source/Source.h"
#include "Utilities/safe_casting.h"
#include "Utilities/invalidInput_error.h"

//kernel
__global__   void trafficTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, real * pConcArray, uint* juncInCellIndices, bool* juncCarCanEnter,
	int* juncCarsOnJunction, uint* juncAlreadyMoved, uint* juncOldSpeeds, real* sinkCarBlockedPossibilities,
	uint size_road, uint size_juncInCells, uint maxVelocity, uint maxAcceleration, uint safetyDistance, bool useSlowToStart, real slowStartPossibility, real dawdlePossibility, curandState *state);

__global__	void sourceTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, uint* sourceIndices, real* sinkCarBlockedPossibilities,
	float* sourcePossibilities, real * pConcArray, uint maxVelocity, uint safetyDistance, uint size_sources, curandState *state);

__global__ void junctionTimestepKernel(int* juncCarsOnJunction, uint* juncInCellIndices, int* juncOutCellIndices, uint* juncStartInIncells, uint* juncStartInOutcells,
	uint* juncAlreadyMoved, int* juncCarCanNotEnterThisOutCell, bool* juncOutCellIsOpen, uint* juncOldSpeeds, bool* juncCarCanEnter, uint* juncTrafficLightSwitchTime,
	int* roadCurrent, int* roadNext, int* neighbors, real * pConcArray, real* sinkCarBlockedPossibilities, uint safetyDistance,
	uint size_juncInCells, uint size_juncOutCells, uint size_junctions, uint numTimestep, curandState *state);

__global__ void randomSetupKernel(curandState *state, uint size);

__global__ void resetNext(int* roadNext, uint size_road);

//device functions for movement
__device__ inline uint getJunctionInCellsVectorIndex(uint * juncInCellIndices, uint size_juncInCells, uint cell);

__device__ inline uint getGapAfterOutCell(int* roadCurrent, int* neighbors, real* sinkCarBlockedPossibilities, int sourceIndex, uint speed, uint safetyDistance, curandState* state);

__device__ inline void carStaysOnJunction(int* juncCarsOnJunction, uint* juncInCellIndices, uint* juncAlreadyMoved, uint*juncOldSpeeds, real* pConcArray, uint inCellVectorIndex);



//device functions for concentrations
__device__ inline real calcConcentration(uint oldSpeed, uint newSpeed);

__device__ inline void putConcIntoArray(real * pConcArray, uint oldSpeed, uint newSpeed, uint newIndex);

__device__ inline void addConcToArray(real * pConcArray, uint oldSpeed, uint newSpeed, uint newIndex);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TrafficTimestep::TrafficTimestep(std::shared_ptr<RoadNetworkData> road, real * pConcArray)
{
	//calculate sizes
	this->size_roads = road->roadLength;
	this->size_junctions = castSizeT_Uint(road->junctions.size());
	this->size_sources = castSizeT_Uint(road->sources.size());
	this->size_sinks = castSizeT_Uint(road->sinks.size());

	//set attributes
	this->maxVelocity = road->maxVelocity;
	this->safetyDistance = road->safetyDistance;
	this->dawdlePossibility = road->dawdlePossibility;
	this->useSlowToStart = road->useSlowToStart;
	this->slowStartPossibility = road->slowStartPossibility;
	this->maxAcceleration = road->maxAcceleration;

	//prepare road
	this->neighbors = road->neighbors;
	this->roadCurrent = *road->pcurrent;

	this->roadNext.resize(size_roads);
	thrust::fill(roadNext.begin(), roadNext.end(), -1);

	switchCurrentNext();

	//prepare junctions
	combineJuncInCellIndices(road->junctions);
	combineJuncOutCellIndices(road->junctions);
	combineJuncCarCanNotEnterThisOutCell(road->junctions);
	combineUseTrafficLights(road->junctions);

	initJuncCarCanEnter();
	initJuncCarsOnJunction();
	initJuncAlreadyMoved();
	initjuncOutCellIsOpen();
	initJuncOldSpeeds();

	//prepare sinks
	combineSinkBlockedPossibilities(road->sinks);

	//prepare sources
	combineSourcePossibilities(road->sources);
	combineSourceIndices(road->sources);

	//preapre grid dimensions
	calculateTrafficTimestepKernelDimensions();
	calculateJunctionKernelDimensions();
	calculateSourceKernelDimensions();

	////setUp random numbers
	checkCudaErrors(cudaMalloc((void **)&statesRoad, gridRoad.x * gridRoad.y * threadsRoads.x * sizeof(curandState)));
	randomSetupKernel << <gridRoad, threadsRoads >> > (statesRoad, size_roads);
	cudaDeviceSynchronize();
	getLastCudaError("random_setup_kernel for roads execution failed");

	checkCudaErrors(cudaMalloc((void **)&statesJunctions, gridJunctions.x * gridJunctions.y * threadsJunctions.x * sizeof(curandState)));
	randomSetupKernel << <gridJunctions, threadsJunctions >> > (statesJunctions, size_junctions);
	cudaDeviceSynchronize();
	getLastCudaError("random_setup_kernel for junctions execution failed");

	checkCudaErrors(cudaMalloc((void **)&statesSources, gridSources.x * gridSources.y * threadsSources.x * sizeof(curandState)));
	randomSetupKernel << <gridSources, threadsSources >> > (statesSources, size_sources);
	cudaDeviceSynchronize();
	getLastCudaError("random_setup_kernel for sources execution failed");

	//prepare ConcWriter
	if (pConcArray == nullptr) checkCudaErrors(cudaMalloc((void **)&this->pConcArray, size_roads * sizeof(real)));
	else this->pConcArray = pConcArray;
	checkCudaErrors(cudaMemset(this->pConcArray, 0.0, size_t(size_roads) * sizeof(real)));
}


void TrafficTimestep::run(std::shared_ptr<RoadNetworkData> road)
{
	switchCurrentNext();

	//reset Junction open outcells
	resetOutCellIsOpen();

	callTrafficTimestepKernel();
	getLastCudaError("trafficTimestepKernel execution failed");

	callJunctionTimestepKernel();
	getLastCudaError("junctionTimestepKernel execution failed");

	callSourceTimestepKernel();
	getLastCudaError("sourceTimestepKernel execution failed");

	numTimestep++;
}


void TrafficTimestep::switchCurrentNext()
{
	if (timestepIsEven) timestepIsEven = false;
	else timestepIsEven = true;
}


void TrafficTimestep::copyCurrentDeviceToHost(std::shared_ptr<RoadNetworkData> road)
{
	if (timestepIsEven) thrust::copy(roadCurrent.begin(), roadCurrent.end(), road->pcurrent->begin());
	else thrust::copy(roadNext.begin(), roadNext.end(), road->pcurrent->begin());

	//checkCudaErrors(cudaMemcpy(&(road->conc[0]), pConcArray, size_roads * sizeof(real), cudaMemcpyDeviceToHost)); //dispConcFromGPU
}


void TrafficTimestep::cleanUp()
{
	cudaFree(statesRoad);
	cudaFree(statesJunctions);
	cudaFree(statesSources);
}


void TrafficTimestep::callTrafficTimestepKernel()
{
	int* pRoadCurrent = roadCurrent.data().get();
	int* pRoadNext = roadNext.data().get();
	if (!timestepIsEven) {
		pRoadCurrent = roadNext.data().get();
		pRoadNext = roadCurrent.data().get();
	}

	resetNext << < gridRoad, threadsRoads >> > (
		pRoadNext,
		size_roads
		);

	trafficTimestepKernel << < gridRoad, threadsRoads >> > (
		pRoadCurrent,
		pRoadNext,
		neighbors.data().get(),
		pConcArray,
		juncInCellIndices.data().get(),
		juncCarCanEnter.data().get(),
		juncCarsOnJunction.data().get(),
		juncAlreadyMoved.data().get(),
		juncOldSpeeds.data().get(),
		sinkCarBlockedPossibilities.data().get(),
		size_roads,
		size_juncInCells,
		maxVelocity,
		maxAcceleration,
		safetyDistance,
		useSlowToStart,
		slowStartPossibility,
		dawdlePossibility,
		statesRoad);
}


void TrafficTimestep::callSourceTimestepKernel()
{
	int*pRoadCurrent = roadCurrent.data().get();
	int*pRoadNext = roadNext.data().get();
	if (!timestepIsEven) {
		pRoadCurrent = roadNext.data().get();
		pRoadNext = roadCurrent.data().get();
	}

	sourceTimestepKernel << < gridSources, threadsSources >> > (
		pRoadCurrent,
		pRoadNext,
		neighbors.data().get(),
		sourceIndices.data().get(),
		sinkCarBlockedPossibilities.data().get(),
		sourcePossibilities.data().get(),
		pConcArray,
		maxVelocity,
		safetyDistance,
		size_sources,
		statesSources);
}


void TrafficTimestep::callJunctionTimestepKernel()
{
	int*pRoadCurrent = roadCurrent.data().get();
	int*pRoadNext = roadNext.data().get();
	if (!timestepIsEven) {
		pRoadCurrent = roadNext.data().get();
		pRoadNext = roadCurrent.data().get();
	}

	junctionTimestepKernel << < gridJunctions, threadsJunctions >> > (
		juncCarsOnJunction.data().get(),
		juncInCellIndices.data().get(),
		juncOutCellIndices.data().get(),
		juncStartInIncells.data().get(),
		juncStartInOutcells.data().get(),
		juncAlreadyMoved.data().get(),
		juncCarCanNotEnterThisOutCell.data().get(),
		juncOutCellIsOpen.data().get(),
		juncOldSpeeds.data().get(),
		juncCarCanEnter.data().get(),
		juncTrafficLightSwitchTime.data().get(),
		pRoadCurrent,
		pRoadNext,
		neighbors.data().get(),
		pConcArray,
		sinkCarBlockedPossibilities.data().get(),
		safetyDistance,
		size_juncInCells,
		size_juncOutCells,
		size_junctions,
		numTimestep,
		statesJunctions);
}



__global__ void resetNext(int * roadNext, uint size_road)
{
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;
	const uint index = nx*(ny*z + y) + x;

	if (index >= size_road) return;

	roadNext[index] = -1;
}



__global__ void trafficTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, real * pConcArray, uint* juncInCellIndices, bool* juncCarCanEnter,
	int* juncCarsOnJunction, uint* juncAlreadyMoved, uint* juncOldSpeeds, real* sinkCarBlockedPossibilities,
	uint size_road, uint size_juncInCells, uint maxVelocity, uint maxAcceleration, uint safetyDistance, bool useSlowToStart, real slowStartPossibility, real dawdlePossibility, curandState *state)
{
	//////////////////////////////////////////////////////////////////////////
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;
	const uint index = nx*(ny*z + y) + x;

	if (index >= size_road) return;
	////////////////////////////////////////////////////////////////////////////////
	if (roadCurrent[index] < 0) return;

	//reset concentrations
	pConcArray[index] = 0.0;

	//printf("index %d ", index);

	uint speed = roadCurrent[index];


	//// accelerate car ////////////////////////////////////////////////////////////////////
	if (speed < maxVelocity) {
		if (speed <= maxVelocity - maxAcceleration)
			speed += maxAcceleration;
		else
			speed = maxVelocity;
	}


	//////// brake car /////////////////////////////////////////////////////////////////////////
	//getGapAfterCar
	uint gap = speed;
	uint idx = 0;
	int neighbor = neighbors[index];
	uint currentCell = index;	

	for (uint i = 0; i < (speed + safetyDistance); i++) {

		//sink
		if (neighbor <= -2000) {
			curandState localState = state[index];
			float random = curand_uniform(&localState);

			if (i <= speed && !(random < sinkCarBlockedPossibilities[(neighbor + 2000)*-1]))  gap = speed;
			else gap = i;

			state[index] = localState;
			break;
		}

		//junction
		if (neighbor <= -1000 && neighbor > -2000) {
			idx = getJunctionInCellsVectorIndex(juncInCellIndices, size_juncInCells, currentCell);
			if (juncCarCanEnter[idx] && i <= speed) gap = speed;
			else gap = i;
			break;
		}

		//car in Cell
		if (roadCurrent[neighbor] > -1) {
			if (i <= safetyDistance) gap = 0;
			else gap = i - safetyDistance;
			break;
		}

		//empty cell -> get next neighbor, update currentCell
		currentCell = neighbor;
		neighbor = neighbors[neighbor];
	}

	//brake
	if (speed > gap)
		speed = gap;


	////// dawdleCar ///////////////////////////////////////////////////////////////////////////
	curandState localState = state[index];
	float random = curand_uniform(&localState);
	state[index] = localState;

	//Barlovic / SlowToStart
	if (useSlowToStart == true && roadCurrent[index] == 0) {
		if (random < slowStartPossibility) speed = 0;
	}
	else if (random < dawdlePossibility) //Standard NaSch
		if (speed >= maxAcceleration)
			speed -= maxAcceleration;
		else
			speed = 0;


	////// moveCar /////////////////////////////////////////////////////////////////////////////
	if (speed == 0) {
		(roadNext)[index] = 0;
		putConcIntoArray(pConcArray, roadCurrent[index], speed, index);
		return;
	}

	neighbor = neighbors[index];
	currentCell = index;

	// iterateNeighborsInMove
	uint numberOfCellsMoved = 0;
	for (uint i = 2; i <= speed; i++) {
		if (neighbor >= 0) {
			currentCell = neighbor;
			neighbor = neighbors[neighbor];
			++numberOfCellsMoved;
		}
		else
			break;
	}

	if (neighbor <= -2000) return;

	if (neighbor <= -1000 && neighbor > -2000) {
		//registerCar

		juncCarsOnJunction[idx] = speed;
		//getCarsOnJunction[index] = speed -2; //all cars, which enter the junction have to slow down
		juncOldSpeeds[idx] = roadCurrent[index];
		juncCarCanEnter[idx] = false;
		juncAlreadyMoved[idx] = numberOfCellsMoved;
		return;
	}

	if (neighbor >= 0) {
		roadNext[neighbor] = speed;
		putConcIntoArray(pConcArray, roadCurrent[index], speed, neighbor);
	}
}



__global__ void sourceTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, uint* sourceIndices, real* sinkCarBlockedPossibilities, float* sourcePossibilities, real* pConcArray,
	uint maxVelocity, uint safetyDistance, uint size_sources, curandState *state)
{
	//////////////////////////////////////////////////////////////////////////
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;

	const uint index = nx*(ny*z + y) + x;
	////////////////////////////////////////////////////////////////////////////////

	if (index >= size_sources) return;
	////////////////////////////////////////////////////////////////////////////////

	int sourceIndex = sourceIndices[index];

	uint gap = getGapAfterOutCell(roadCurrent, neighbors, sinkCarBlockedPossibilities, sourceIndex, maxVelocity, safetyDistance, state);
	if (gap > 0) {
		//get car with random speed
		curandState localState = state[index];
		if (curand_uniform(&localState) < sourcePossibilities[index]) {
			unsigned int speed = ceilf(curand_uniform(&localState) * (maxVelocity + 1)) - 1;
			roadNext[sourceIndex] = speed;
			putConcIntoArray(pConcArray, speed, speed, sourceIndex);
		}
		state[index] = localState;
	}
}



__global__ void junctionTimestepKernel(int* juncCarsOnJunction, uint* juncInCellIndices, int* juncOutCellIndices, uint* juncStartInIncells, uint* juncStartInOutcells,
	uint* juncAlreadyMoved, int* juncCarCanNotEnterThisOutCell, bool* juncOutCellIsOpen, uint* juncOldSpeeds, bool* juncCarCanEnter, uint* juncTrafficLightSwitchTime,
	int* roadCurrent, int* roadNext, int* neighbors, real* pConcArray, real* sinkCarBlockedPossibilities, uint safetyDistance,
	uint size_juncInCells, uint size_juncOutCells, uint size_junctions, uint numTimestep, curandState* state) {
	//////////////////////////////////////////////////////////////////////////
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;

	const uint index = nx*(ny*z + y) + x;

	////////////////////////////////////////////////////////////////////////////////

	if (index >= size_junctions) return;
	//printf("junctionKernelIndex %d ", index);

	////////////////////////////////////////////////////////////////////////////////

	//calc indices

	uint inCellsSize = 0;
	uint firstInCellIndex = juncStartInIncells[index];
	if (index < size_junctions - 1) inCellsSize = juncStartInIncells[index + 1] - firstInCellIndex;
	else inCellsSize = size_juncInCells - firstInCellIndex;

	uint outCellSize = 0;
	uint firstOutCellIndex = juncStartInOutcells[index];
	if (index < size_junctions - 1) outCellSize = juncStartInOutcells[index + 1] - firstOutCellIndex;
	else outCellSize = size_juncOutCells - firstOutCellIndex;


	//// loop through all cars /////////////////////////////////////////////////////////////////////////////

	for (uint inCellVectorIndex = firstInCellIndex; inCellVectorIndex < firstInCellIndex + inCellsSize; inCellVectorIndex++) {

		if (juncCarsOnJunction[inCellVectorIndex] >= 0) {

			//// applyRules /////////////////////////////////////////////////////////////////////////////
			uint speed = juncCarsOnJunction[inCellVectorIndex];
			if (speed == 0 && juncAlreadyMoved[inCellVectorIndex] == 0) speed += 1;
			int remainingDist = speed - static_cast<int>(juncAlreadyMoved[inCellVectorIndex]);

			//printf(" incell %d, speed %d, remainingDist %d, alreadyMoved %d ", inCellVectorIndex, speed, remainingDist, juncAlreadyMoved[inCellVectorIndex]);

			if (remainingDist == 0) { //car can't leave the junction
				carStaysOnJunction(juncCarsOnJunction, juncInCellIndices, juncAlreadyMoved, juncOldSpeeds, pConcArray, inCellVectorIndex);
				continue;
			}

			else {

				//// calc numberOfPossibleOutCells ////////////////////////////////////////////////////////////////
				uint numberOfPossibleOutCells = 0;

				for (uint outCellIndex = firstOutCellIndex; outCellIndex < firstOutCellIndex + outCellSize; outCellIndex++) 
					if (juncCarCanNotEnterThisOutCell[inCellVectorIndex] != juncOutCellIndices[outCellIndex] && juncOutCellIsOpen[outCellIndex] == true)
						numberOfPossibleOutCells++;


				if (numberOfPossibleOutCells == 0)  //car can't leave the junction
				{
					carStaysOnJunction(juncCarsOnJunction, juncInCellIndices, juncAlreadyMoved, juncOldSpeeds, pConcArray, inCellVectorIndex);
					continue;
				}

				//// chooseOutCell ///////////////////////////////////////////////////////////////////////////////
				int chosenCell = -1;
				curandState localState = state[index];
				uint random = ceilf(curand_uniform(&localState) * numberOfPossibleOutCells);
				for (uint outCellVectorIndex = firstOutCellIndex; outCellVectorIndex < firstOutCellIndex + outCellSize; outCellVectorIndex++) {
					if (juncCarCanNotEnterThisOutCell[inCellVectorIndex] != juncOutCellIndices[outCellVectorIndex] && juncOutCellIsOpen[outCellVectorIndex] == true) {
						if (random == 1) {
							chosenCell = juncOutCellIndices[outCellVectorIndex];
							juncOutCellIsOpen[outCellVectorIndex] = false;
							break;
						}
						random--;
					}
				}
				state[index] = localState;

				//// brakeCar ////////////////////////////////////////////////////////////////////////////////////
				if (chosenCell < 0);
				uint gap = getGapAfterOutCell(roadCurrent, neighbors, sinkCarBlockedPossibilities, chosenCell, speed, safetyDistance, state);
				if (gap < remainingDist) {
					if (gap > speed) gap = speed;
					speed = speed - remainingDist + gap;
					remainingDist = gap;
				}

				//// moveCar /////////////////////////////////////////////////////////////////////////////////////
				if (remainingDist <= 0) { //car can't leave the junction
					carStaysOnJunction(juncCarsOnJunction, juncInCellIndices, juncAlreadyMoved, juncOldSpeeds, pConcArray, inCellVectorIndex);
					continue;
				}

				if (remainingDist > 0) {

					if (remainingDist == 1) {
						roadNext[chosenCell] = speed;
						putConcIntoArray(pConcArray, juncOldSpeeds[inCellVectorIndex], speed, chosenCell);
						juncCarsOnJunction[inCellVectorIndex] = -1;
						juncCarCanEnter[inCellVectorIndex] = true;
						break;
					}

					//iterate through neighbors
					int neighbor = chosenCell;
					for (uint i = 2; i <= remainingDist; i++) {
						if (neighbor >= 0) {
							chosenCell = neighbor;
							neighbor = neighbors[neighbor];
						}
						else
							break;
					}

					if (neighbor >= 0) {
						roadNext[neighbor] = speed;
						putConcIntoArray(pConcArray, juncOldSpeeds[inCellVectorIndex], speed, neighbor);
					}

					juncCarsOnJunction[inCellVectorIndex] = -1;
					juncCarCanEnter[inCellVectorIndex] = true;
				}
			}
		}
	}


	//generate red TrafficLights

	if (juncTrafficLightSwitchTime[index] > 0) {

		uint halfNumStreets = (uint)(std::floor((float)inCellsSize * 0.5f));

		if ((uint)(std::floor((float)numTimestep / (float)juncTrafficLightSwitchTime[index])) % 2 == 0) {
			for (uint i = firstInCellIndex; i < firstInCellIndex + halfNumStreets; i++)
				juncCarCanEnter[i] = false;

			if (numTimestep % juncTrafficLightSwitchTime[index] == 0) //first timestep with green light --> open the streets that were closed before
				for (uint i = firstInCellIndex + halfNumStreets; i < firstInCellIndex + inCellsSize; i++)
					if (juncCarsOnJunction[i] == -1)
						juncCarCanEnter[i] = true;
		}
		else {
			for (uint i = firstInCellIndex + halfNumStreets; i < firstInCellIndex + inCellsSize; i++)
				juncCarCanEnter[i] = false;

			if (numTimestep % juncTrafficLightSwitchTime[index] == 0) //first timestep with green light --> open the streets that were closed before
				for (uint i = firstInCellIndex; i < firstInCellIndex + halfNumStreets; i++)
					if (juncCarsOnJunction[i] == -1)
						juncCarCanEnter[i] = true;
		}
	}

}


__device__ inline void carStaysOnJunction(int* juncCarsOnJunction, uint* juncInCellIndices, uint* juncAlreadyMoved, uint*juncOldSpeeds, real* pConcArray, uint inCellVectorIndex) {
	addConcToArray(pConcArray, juncOldSpeeds[inCellVectorIndex], 0, juncInCellIndices[inCellVectorIndex]);
	juncCarsOnJunction[inCellVectorIndex] = 0;
	juncAlreadyMoved[inCellVectorIndex] = 0;
	juncOldSpeeds[inCellVectorIndex] = 0;
}


__device__ inline uint getGapAfterOutCell(int* roadCurrent, int* neighbors, real* sinkCarBlockedPossibilities, int outCellIndex, uint speed, uint safetyDistance, curandState* state)
{
	if (roadCurrent[outCellIndex] > -1)
		return 0;

	for (uint i = 0; i < (speed + safetyDistance); i++) {
		//sink
		if (outCellIndex <= -2000) {
			const uint index = blockDim.x*(gridDim.x*blockIdx.y + blockIdx.x) + threadIdx.x;
			curandState localState = state[index];
			float randomNumber = curand_uniform(&localState);
			state[index] = localState;
			if (i <= speed && !(randomNumber < sinkCarBlockedPossibilities[(outCellIndex + 2000)*-1]))
				return speed;
			return i;
		}
		//junction
		if (outCellIndex <= -1000)
			return i;

		//car in Cell
		if (roadCurrent[outCellIndex] > -1) {
			if (i <= safetyDistance) return 0;
			return i - safetyDistance;
		}

		//empty cell -> get next neighbor
		outCellIndex = neighbors[outCellIndex];
	}
	return speed;
}


inline __device__ uint getJunctionInCellsVectorIndex(uint* juncInCellIndices, uint size_juncInCells, uint cell) {
	for (uint i = 0; i < size_juncInCells; i++)
		if (juncInCellIndices[i] == cell)
			return i;
	//TODO real Error
	printf("no matching incoming cell to a junction found in: getJunctionInCellsVectorIndex()");
	return 65000;
}


__global__ void randomSetupKernel(curandState *state, uint size) {
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;

	const uint index = nx*(ny*z + y) + x;

	if (index >= size) return;

	curand_init((unsigned long long)clock() + index, index, 0, &state[index]);
}


__device__ real calcConcentration(uint oldSpeed, uint newSpeed)
{
	//printf("newIndex %d ", newIndex );
	if (oldSpeed == 0 && newSpeed > 0) //Start
		return 1.0f;
	else if (oldSpeed == 0 && newSpeed == 0) //Idle
		return 0.25f;
	else if (newSpeed == oldSpeed) //Drive
		return 0.5f;
	else if (newSpeed > oldSpeed) //Accelerate
		return 0.75f;
	else if (newSpeed < oldSpeed) //Brake
		return 0.45f;
	else
		printf("couldn't choose driving state in calcConcentration");
	return -1;
}


__device__ void putConcIntoArray(real * pConcArray, uint oldSpeed, uint newSpeed, uint newIndex)
{
	pConcArray[newIndex] = calcConcentration(oldSpeed, newSpeed);
}


__device__ void addConcToArray(real * pConcArray, uint oldSpeed, uint newSpeed, uint newIndex)
{
	pConcArray[newIndex] += calcConcentration(oldSpeed, newSpeed);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TrafficTimestep::combineJuncInCellIndices(std::vector<std::shared_ptr<Junction>> &junctions)
{
	for (auto& j : junctions) {
		if (juncInCellIndices.size() == 0) juncStartInIncells.push_back(0);
		else juncStartInIncells.push_back(castSizeT_Uint(juncInCellIndices.size()));

		for (uint i : j->getInCellIndices())
			this->juncInCellIndices.push_back(i);
	}
	size_juncInCells = castSizeT_Uint(juncInCellIndices.size());
}


void TrafficTimestep::combineJuncOutCellIndices(std::vector<std::shared_ptr<Junction>> &junctions)
{
	for (auto& j : junctions) {
		if (juncOutCellIndices.size() == 0) juncStartInOutcells.push_back(0);
		else juncStartInOutcells.push_back(castSizeT_Uint(juncOutCellIndices.size()));

		for (uint i : j->getOutCellIndices())
			this->juncOutCellIndices.push_back(i);
	}
	size_juncOutCells = castSizeT_Uint(juncOutCellIndices.size());
}


void TrafficTimestep::initJuncCarCanEnter()
{
	juncCarCanEnter.resize(size_juncInCells);
	thrust::fill(juncCarCanEnter.begin(), juncCarCanEnter.end(), true);
}

void TrafficTimestep::initJuncCarsOnJunction()
{
	juncCarsOnJunction.resize(size_juncInCells);
	thrust::fill(juncCarsOnJunction.begin(), juncCarsOnJunction.end(), -1);
}


void TrafficTimestep::combineSinkBlockedPossibilities(std::vector<std::shared_ptr<Sink>> &sinks)
{
	for (auto& s : sinks)
		this->sinkCarBlockedPossibilities.push_back(s->getPossibilityBeingBlocked());
}


void TrafficTimestep::combineSourcePossibilities(std::vector<std::shared_ptr<Source>> &sources) {
	for (auto& s : sources)
		this->sourcePossibilities.push_back(s->getPossibility());
}

void TrafficTimestep::combineSourceIndices(std::vector<std::shared_ptr<Source>> &sources)
{
	for (auto& s : sources)
		this->sourceIndices.push_back(s->getIndex());
}


void TrafficTimestep::initJuncAlreadyMoved()
{
	juncAlreadyMoved.resize(size_juncInCells);
	thrust::fill(juncAlreadyMoved.begin(), juncAlreadyMoved.end(), 0);
}

void TrafficTimestep::initJuncOldSpeeds()
{
	juncOldSpeeds.resize(size_juncInCells);
	thrust::fill(juncOldSpeeds.begin(), juncOldSpeeds.end(), 0);
}

void TrafficTimestep::combineUseTrafficLights(std::vector<std::shared_ptr<Junction>>& junctions)
{
	for (auto& j : junctions)
		juncTrafficLightSwitchTime.push_back(j->getTrafficLightSwitchTime());
}

void TrafficTimestep::initjuncOutCellIsOpen()
{
	juncOutCellIsOpen.resize(juncOutCellIndices.size());
	resetOutCellIsOpen();
}

void TrafficTimestep::combineJuncCarCanNotEnterThisOutCell(std::vector<std::shared_ptr<Junction>>& junctions)
{
	for (auto& j : junctions)
		for (int i : j->getCarCanNotEnterThisOutCell())
			this->juncCarCanNotEnterThisOutCell.push_back(i);

	if (juncCarCanNotEnterThisOutCell.size() < size_juncInCells)
		juncCarCanNotEnterThisOutCell.push_back(-2);
}

uint TrafficTimestep::getNumCarsOnJunctions()
{
	uint num = 0;
	for (int car : juncCarsOnJunction)
		if (car >= 0) ++num;
	return num;
}

void TrafficTimestep::resetOutCellIsOpen()
{
	thrust::fill(juncOutCellIsOpen.begin(), juncOutCellIsOpen.end(), true);
}


void TrafficTimestep::calculateTrafficTimestepKernelDimensions()
{
	unsigned int numberOfThreads = 64;
	int Grid = (size_roads / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid > 512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 gridRoad1(Grid1, Grid2);
	this->gridRoad = gridRoad1;
	dim3 threadsRoads1(numberOfThreads, 1, 1);
	this->threadsRoads = threadsRoads1;
}

void TrafficTimestep::calculateJunctionKernelDimensions()
{
	unsigned int numberOfThreads = 32;
	int Grid = (size_junctions / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid > 512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 gridJunctions1(Grid1, Grid2);
	this->gridJunctions = gridJunctions1;
	dim3 threadsJunctions1(numberOfThreads, 1, 1);
	this->threadsJunctions = threadsJunctions1;
}

void TrafficTimestep::calculateSourceKernelDimensions()
{
	unsigned int numberOfThreads = 32;
	int Grid = (size_sinks / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid > 512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 gridSources1(Grid1, Grid2);
	this->gridSources = gridSources1;
	dim3 threadsSources1(numberOfThreads, 1, 1);
	this->threadsSources = threadsSources1;
}
