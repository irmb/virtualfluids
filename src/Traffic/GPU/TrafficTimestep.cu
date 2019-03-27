#include "TrafficTimestep.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include <cmath>
#include <sstream>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>

#include <iomanip>

//#include "TrafficMovement.h"
#include "RoadNetwork/RoadNetworkData.h"
#include "Junction/Junction.h"
#include "Sink/Sink.h"
#include "Source/Source.h"
#include "Utilities/safe_casting.h"
#include "Utilities/invalidInput_error.h"

__global__   void trafficTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, int* oldSpeeds, uint* juncInCellIndices, bool* juncCarCanEnter,
	int* juncCarsOnJunction, uint* juncAlreadyMoved, uint* juncOldSpeeds, bool* sinkCarCanEnter,
	uint size_road, uint maxVelocity, uint maxAcceleration, uint safetyDistance, bool useSlowToStart, real slowStartPossibility, real dawdlePossibility, curandState *curandstate);

__global__	void sourceTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, uint* sourceIndices, bool* sinkCarCanEnter,
	float* sourcePossibilities, int* oldSpeeds, uint maxVelocity, uint safetyDistance, uint size_sources, curandState *curandstate);

__global__ void junctionTimestepKernel(int* juncCarsOnJunction, int* juncOutCellIndices, uint* juncStartInIncells, uint* juncStartInOutcells,
	uint* juncAlreadyMoved, int* juncCarCanNotEnterThisOutCell, bool* juncOutCellIsOpen, uint* juncOldSpeeds, bool* juncCarCanEnter,
	int* roadCurrent, int* roadNext, int* neighbors, int* oldSpeeds, bool* sinkCarCanEnter, uint safetyDistance,
	uint size_juncInCells, uint size_juncOutCells, uint size_junctions, curandState *curandstate);

__global__ void random_setup_kernel(curandState *state);

__device__ inline uint getJunctionInCellsVectorIndex(uint * juncInCellIndices, uint size_road, uint cell);

__device__ inline uint getGapAfterOutCell(int* roadCurrent, int* neighbors, bool* sinkCarCanEnter, int sourceIndex, uint speed, uint safetyDistance);



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TrafficTimestep::TrafficTimestep(std::shared_ptr<RoadNetworkData> road)
{
	//calculate sizes
	size_roads = road->roadLength;
	size_junctions = castSizeT_Uint(road->junctions.size());
	size_sources = castSizeT_Uint(road->sources.size());
	size_sinks = castSizeT_Uint(road->sinks.size());

	//set attributes
	maxVelocity = road->maxVelocity;
	safetyDistance = road->safetyDistance;
	dawdlePossibility = road->dawdlePossibility;
	useSlowToStart = road->useSlowToStart;
	slowStartPossibility = road->slowStartPossibility;
	maxAcceleration = road->maxAcceleration;

	//prepare road
	this->neighbors.resize(size_roads);
	this->neighbors = road->neighbors;
	roadCurrent.resize(size_roads);
	roadNext.resize(size_roads);
	oldSpeeds.resize(size_roads);

	//prepare junctions
	combineJuncInCellIndices(road->junctions);
	combineJuncOutCellIndices(road->junctions);
	combineJuncCarCanNotEnterThisOutCell(road->junctions);

	initJuncCarCanEnter();
	initJuncCarsOnJunction();
	initJuncAlreadyMoved();
	initJuncOldSpeeds();
	initjuncOutCellIsOpen();

	//prepare sinks
	sinkCarCanEnter.resize(size_sinks);

	//prepare sources
	combineSourcePossibilities(road->sources);
	combineSourceIndices(road->sources);

}

void TrafficTimestep::run(std::shared_ptr<RoadNetworkData> road)
{
	//copy road to device, reset next
	this->roadCurrent = *(road->pcurrent);
	std::fill(roadNext.begin(), roadNext.end(), -1);

	//reset Junction open outcells
	resetOutCellIsOpen();

	//copy sinks to device
	combineSinkCarCanEnterSink(road->sinks);

	callTrafficTimestepKernel();
	cudaDeviceSynchronize();
	getLastCudaError("trafficTimestepKernel execution failed");
	callJunctionKernel();
	cudaDeviceSynchronize();
	getLastCudaError("junctionTimestepKernel execution failed");
	callSourceKernel();
	cudaDeviceSynchronize();
	getLastCudaError("sourceTimestepKernel execution failed");

	thrust::copy(roadNext.begin(), roadNext.end(), (*(road->pnext)).begin());
	thrust::copy(roadCurrent.begin(), roadCurrent.end(), (*(road->pcurrent)).begin());
	thrust::copy(oldSpeeds.begin(), oldSpeeds.end(), road->oldSpeeds.begin());
}

void TrafficTimestep::callTrafficTimestepKernel()
{
	unsigned int numberOfThreads = 128;
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
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	//setUp random numbers
	curandState *d_state;
	cudaMalloc(&d_state, sizeof(curandState));
	random_setup_kernel << <1, 1 >> >(d_state);

	trafficTimestepKernel << < grid, threads >> > (
		roadCurrent.data().get(),
		roadNext.data().get(),
		neighbors.data().get(),
		oldSpeeds.data().get(),
		juncInCellIndices.data().get(),
		juncCarCanEnter.data().get(),
		juncCarsOnJunction.data().get(),
		juncAlreadyMoved.data().get(),
		juncOldSpeeds.data().get(),
		sinkCarCanEnter.data().get(),
		size_roads,
		maxVelocity,
		maxAcceleration,
		safetyDistance,
		useSlowToStart,
		slowStartPossibility,
		dawdlePossibility,
		d_state);
}

void TrafficTimestep::callSourceKernel()
{
	//calculate	grid and threads
	unsigned int numberOfThreads = 128;
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
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	//setUp random numbers
	curandState *d_state;
	cudaMalloc(&d_state, sizeof(curandState));
	random_setup_kernel << <1, 1 >> >(d_state);

	sourceTimestepKernel << < grid, threads >> > (
		roadCurrent.data().get(),
		roadNext.data().get(),
		neighbors.data().get(),
		sourceIndices.data().get(),
		sinkCarCanEnter.data().get(),
		sourcePossibilities.data().get(),
		oldSpeeds.data().get(),
		maxVelocity,
		safetyDistance,
		size_sources,
		d_state);
}


void TrafficTimestep::callJunctionKernel()
{
	//calculate	grid and threads
	unsigned int numberOfThreads = 128;
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
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	//setUp random numbers
	curandState *d_state;
	cudaMalloc(&d_state, sizeof(curandState));
	random_setup_kernel << <1, 1 >> >(d_state);
	
	junctionTimestepKernel << < grid, threads >> > (
		juncCarsOnJunction.data().get(),
		juncOutCellIndices.data().get(),
		juncStartInIncells.data().get(),
		juncStartInOutcells.data().get(),
		juncAlreadyMoved.data().get(),
		juncCarCanNotEnterThisOutCell.data().get(),
		juncOutCellIsOpen.data().get(),
		juncOldSpeeds.data().get(),
		juncCarCanEnter.data().get(),
		roadCurrent.data().get(),
		roadNext.data().get(),
		neighbors.data().get(),
		oldSpeeds.data().get(),
		sinkCarCanEnter.data().get(),
		safetyDistance,
		size_juncInCells,
		size_juncOutCells,
		size_junctions,
		d_state);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void trafficTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, int* oldSpeeds, uint* juncInCellIndices, bool* juncCarCanEnter, int* juncCarsOnJunction,
	uint* juncAlreadyMoved, uint* juncOldSpeeds, bool* sinkCarCanEnter,
	uint size_road, uint maxVelocity, uint maxAcceleration, uint safetyDistance, bool useSlowToStart, real slowStartPossibility, real dawdlePossibility, curandState *curandstate)
{
	//////////////////////////////////////////////////////////////////////////
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;

	const uint index = nx*(ny*z + y) + x;
	////////////////////////////////////////////////////////////////////////////////

	if (index >= size_road) return;
	////////////////////////////////////////////////////////////////////////////////

	if (roadCurrent[index] < 0) return;

	uint speed = roadCurrent[index];



	//// accelerate car ////////////////////////////////////////////////////////////////////
	if (speed < maxVelocity) {
		if (speed <= maxVelocity - maxAcceleration)
			speed += maxAcceleration;
		else
			speed = maxVelocity;
	}


	//// brake car /////////////////////////////////////////////////////////////////////////

	//getGapAfterCar
	uint gap;
	int neighbor = neighbors[index];
	uint currentCell = index;
	for (uint i = 0; i < (speed + safetyDistance); i++) {

		//sink
		if (neighbor <= -2000) {
			if (i <= speed && sinkCarCanEnter[(neighbor + 2000)*-1]) gap = speed;
			else gap = i;
			break;
		}

		//junction
		if (neighbor <= -1000) {
			if (juncCarCanEnter[getJunctionInCellsVectorIndex(juncInCellIndices, size_road, currentCell)] && i <= speed) gap = speed;
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



	//// dawdleCar ///////////////////////////////////////////////////////////////////////////
	float randomNumber = curand_uniform_double(curandstate);



	//Barlovic / SlowToStart
	if (useSlowToStart == true && roadCurrent[index] == 0) {
		if (randomNumber < slowStartPossibility) speed = 0;
	}
	else if (randomNumber < dawdlePossibility) //Standard NaSch
		if (speed >= maxAcceleration)
			speed -= maxAcceleration;
		else
			speed = 0;



	//// moveCar /////////////////////////////////////////////////////////////////////////////
	if (speed == 0) {
		(roadNext)[index] = 0;
		oldSpeeds[index] = roadCurrent[index];
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


	if (neighbor <= -1000 && neighbor > -2000) {
		//registerCar
		uint idx = getJunctionInCellsVectorIndex(juncInCellIndices, size_road, currentCell);
		juncCarsOnJunction[idx] = speed - 1; //all cars, which enter the junction have to slow down by one increment
		//getCarsOnJunction[index] = 0; //all cars, which enter the junction have to stop
		juncOldSpeeds[idx] = roadCurrent[index];
		juncCarCanEnter[idx] = false;
		juncAlreadyMoved[idx] = numberOfCellsMoved;
		return;
	}

	if (neighbor >= 0) {
		roadNext[neighbor] = speed;
		oldSpeeds[neighbor] = roadCurrent[index];
	}
}




__global__ void sourceTimestepKernel(int* roadCurrent, int* roadNext, int* neighbors, uint* sourceIndices, bool* sinkCarCanEnter, float* sourcePossibilities, int* oldSpeeds, 
	uint maxVelocity, uint safetyDistance, uint size_sources, curandState *curandstate)
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
	
	uint gap = getGapAfterOutCell(roadCurrent, neighbors, sinkCarCanEnter, sourceIndex, maxVelocity, safetyDistance);
	if (gap > 0) {
		//get car with random speed
		if (curand_uniform_double(curandstate) < sourcePossibilities[index]) {
			unsigned int speed = ceilf(curand_uniform(curandstate) * (maxVelocity + 1)) - 1;
			roadNext[sourceIndex] = speed;
			oldSpeeds[sourceIndex] = speed;
		}
	}
}


__global__ void junctionTimestepKernel(int* juncCarsOnJunction, int* juncOutCellIndices, uint* juncStartInIncells, uint* juncStartInOutcells,
	uint* juncAlreadyMoved, int* juncCarCanNotEnterThisOutCell, bool* juncOutCellIsOpen, uint* juncOldSpeeds, bool* juncCarCanEnter,
	int* roadCurrent, int* roadNext, int* neighbors, int* oldSpeeds, bool* sinkCarCanEnter, uint safetyDistance,
	uint size_juncInCells, uint size_juncOutCells, uint size_junctions, curandState* curandstate) {
	//////////////////////////////////////////////////////////////////////////
	const uint x = threadIdx.x;  // Globaler x-Index 
	const uint y = blockIdx.x;   // Globaler y-Index 
	const uint z = blockIdx.y;   // Globaler z-Index 

	const uint nx = blockDim.x;
	const uint ny = gridDim.x;

	const uint index = nx*(ny*z + y) + x;
	////////////////////////////////////////////////////////////////////////////////

	if (index >= size_junctions) return;
	////////////////////////////////////////////////////////////////////////////////

	uint inCellsSize = 0;
	uint firstInCellIndex = juncStartInIncells[index];
	if (index < size_junctions - 1) inCellsSize = juncStartInIncells[index + 1] - firstInCellIndex;
	else inCellsSize = size_juncInCells - firstInCellIndex;

	uint outCellSize = 0;
	uint firstOutCellIndex = juncStartInOutcells[index];
	if (index < size_junctions - 1) outCellSize = juncStartInOutcells[index + 1] - firstOutCellIndex;
	else outCellSize = size_juncOutCells - firstOutCellIndex;

	//loop through all cars
	for (uint inCellVectorIndex = firstInCellIndex; inCellVectorIndex < firstInCellIndex + inCellsSize; inCellVectorIndex++) {

		if (juncCarsOnJunction[inCellVectorIndex] >= 0) {

			//applyRules
			uint speed = juncCarsOnJunction[inCellVectorIndex];
			if (speed == 0 && juncAlreadyMoved[inCellVectorIndex] == 0) speed += 1;
			int remainingDist = speed - static_cast<int>(juncAlreadyMoved[inCellVectorIndex]);

			//printf("incell %d, speed %d, remainingDist %d, alreadyMoved %d", inCellVectorIndex, speed, remainingDist, juncAlreadyMoved[inCellVectorIndex]);
			if (remainingDist == 0) { //car can't leave the junction
				juncCarsOnJunction[inCellVectorIndex] = 0;
				juncAlreadyMoved[inCellVectorIndex] = 0;
				juncOldSpeeds[inCellVectorIndex] = 0;
				continue;
			}

			else {

				//calc numberOfPossibleOutCells
				uint numberOfPossibleOutCells = 0;

				for (uint outCellIndex = firstOutCellIndex; outCellIndex < firstOutCellIndex + outCellSize; outCellIndex++) {
					if (juncCarCanNotEnterThisOutCell[inCellVectorIndex] != juncOutCellIndices[outCellIndex] && juncOutCellIsOpen[outCellIndex] == true)
						numberOfPossibleOutCells++;
				}


				if (numberOfPossibleOutCells == 0)  //car can't leave the junction
				{
					juncCarsOnJunction[inCellVectorIndex] = 0;
					juncAlreadyMoved[inCellVectorIndex] = 0;
					juncOldSpeeds[inCellVectorIndex] = 0;
					continue;
				}

				//choose outcell
				int chosenCell = -1;
				uint random = ceilf(curand_uniform(curandstate) * numberOfPossibleOutCells);
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

				//brakeCar
				if (chosenCell < 0); //TODO cerr
				uint gap = getGapAfterOutCell(roadCurrent, neighbors, sinkCarCanEnter, chosenCell, speed, safetyDistance);
				if (gap < remainingDist) {
					if (gap > speed) gap = speed;
					speed = speed - remainingDist + gap;
					remainingDist = gap;
				}
				
				if (remainingDist <= 0) { //car can't leave the junction
					juncCarsOnJunction[inCellVectorIndex] = 0;
					juncAlreadyMoved[inCellVectorIndex] = 0;
					juncOldSpeeds[inCellVectorIndex] = 0;
					continue;
				} //moveCar

				if (remainingDist > 0) {

					if (remainingDist == 1) {
						roadNext[chosenCell] = speed;
						oldSpeeds[chosenCell] = juncOldSpeeds[inCellVectorIndex];
						juncCarsOnJunction[inCellVectorIndex] = -1;
						juncCarCanEnter[inCellVectorIndex] = true;
						return;
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
						oldSpeeds[neighbor] = juncOldSpeeds[inCellVectorIndex];
					}

					juncCarsOnJunction[inCellVectorIndex] = -1;
					juncCarCanEnter[inCellVectorIndex] = true;
				}
			}

		}
	}
	//TODO writeConc

}



inline __device__ uint getGapAfterOutCell(int* roadCurrent, int* neighbors, bool* sinkCarCanEnter, int sourceIndex, uint speed, uint safetyDistance)
{
	if (roadCurrent[sourceIndex] > -1)
		return 0;

	for (uint i = 0; i < (speed + safetyDistance); i++) {
		//sink
		if (sourceIndex <= -2000) {
			if (i <= speed && sinkCarCanEnter[(sourceIndex + 2000)*-1])
				return speed;
			return i;
		}
		//junction
		if (sourceIndex <= -1000)
			return i;

		//car in Cell
		if (roadCurrent[sourceIndex] > -1) {
			if (i <= safetyDistance) return 0;
			return i - safetyDistance;
		}

		//empty cell -> get next neighbor
		sourceIndex = neighbors[sourceIndex];
	}
	return speed;
}

__device__ uint getJunctionInCellsVectorIndex(uint* juncInCellIndices, uint size_road, uint cell) {
	for (uint i = 0; i < size_road; i++)
		if (juncInCellIndices[i] == cell)
			return i;
	// TODO Error: "no matching incoming cell to a junction found in: getJunctionInCellsVectorIndex()";
	return 65000;
}




__global__ void random_setup_kernel(curandState *state) {

	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	curand_init((unsigned long long)clock(), idx, 0, &state[idx]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void TrafficTimestep::combineJuncInCellIndices(std::vector<std::shared_ptr<Junction>> &junctions)
{
	for (auto& j : junctions) {
		juncStartInIncells.push_back(castSizeT_Uint(juncInCellIndices.size()));
		for (uint i : j->getInCellIndices())
			this->juncInCellIndices.push_back(i);
	}
	size_juncInCells = castSizeT_Uint(juncInCellIndices.size());
}

void TrafficTimestep::initJuncCarCanEnter()
{
	juncCarCanEnter.resize(size_juncInCells);
	std::fill(juncCarCanEnter.begin(), juncCarCanEnter.end(), true);
}

void TrafficTimestep::initJuncCarsOnJunction()
{
	juncCarsOnJunction.resize(size_juncInCells);
	std::fill(juncCarsOnJunction.begin(), juncCarsOnJunction.end(), -1);
}


void TrafficTimestep::combineSinkCarCanEnterSink(std::vector<std::shared_ptr<Sink>> &sinks)
{
	uint it = 0;
	for (auto& s : sinks) {
		this->sinkCarCanEnter[it] = s->carCanEnter();
		it++;
	}
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
	std::fill(juncAlreadyMoved.begin(), juncAlreadyMoved.end(), 0);
}

void TrafficTimestep::initJuncOldSpeeds()
{
	juncOldSpeeds.resize(size_juncInCells);
	std::fill(juncOldSpeeds.begin(), juncOldSpeeds.end(), 0);
}

void TrafficTimestep::initjuncOutCellIsOpen()
{
	juncOutCellIsOpen.resize(juncOutCellIndices.size());
	resetOutCellIsOpen();
}

void TrafficTimestep::combineJuncOutCellIndices(std::vector<std::shared_ptr<Junction>> &junctions)
{
	for (auto& j : junctions) {
		juncStartInOutcells.push_back(castSizeT_Uint(juncOutCellIndices.size()));
		for (uint i : j->getOutCellIndices())
			this->juncOutCellIndices.push_back(i);
	}
	size_juncOutCells = castSizeT_Uint(juncOutCellIndices.size());

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
	std::fill(juncOutCellIsOpen.begin(), juncOutCellIsOpen.end(), true);
}