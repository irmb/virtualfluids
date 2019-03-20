#include "TrafficMovement.h"

#include <iostream>
#include <stdexcept>

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Utilities/RandomHelper.h"
#include "Output/ConcentrationOutwriter.h"
#include "Output/CarDisplay.h"

TrafficMovement::TrafficMovement(std::unique_ptr<RoadNetworkData> road, const real dawdlePossibility)
{
	this->road = std::move(road);

	pcurrent = &this->road->current;
	pnext = &(this->road->next);

	checkCurrentForSafetyDistance();

	initDawdle(dawdlePossibility);
}


TrafficMovement::~TrafficMovement()
{
	pcurrent = NULL;
	pnext = NULL;
	pdummy = NULL;
}


void TrafficMovement::initDawdle(const real dawdlePossibility)
{
	try {
		if (dawdlePossibility >= 0 && dawdlePossibility < 1) {
			this->dawdlePossibility = dawdlePossibility;
		}
		else {
			throw invalidInput_error("The dawdlePossibility should be between 0 and 1.");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


void TrafficMovement::setSlowToStart(const real slowStartPossibility)
{
	try {
		if (slowStartPossibility >= 0 && slowStartPossibility < 1) {
			if (slowStartPossibility > 0) {
				this->slowStartPossibility = slowStartPossibility;
				useSlowToStart = true;
			}
		}
		else {
			throw invalidInput_error("The slowStartPossibility should be between 0 and 1.");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}

void TrafficMovement::setMaxAcceleration(uint maxAcceleration)
{
	this->maxAcceleration = maxAcceleration;
}

void TrafficMovement::setConcentrationOutwriter(std::unique_ptr<ConcentrationOutwriter> writer)
{
	this->concWriter = std::move(writer);
}


void TrafficMovement::setSaveResultsTrue(uint timeSteps)
{
	this->display = std::make_unique<CarDisplay>(&pcurrent, road->safetyDistance);
	if (display != nullptr) display->initResults(timeSteps);
}


const uint TrafficMovement::getNumberOfCars() const
{
	uint num = 0;
	for (auto& junc : road->junctions) {
		num += junc->getNumCarsOnJunction();
	}

	for (auto cell : *pcurrent) {
		if (cell >= 0) ++num;
	}
	return num;
}


const int TrafficMovement::getSpeedAtPosition(uint pos) const
{
	return (*pcurrent)[pos];
}


const uint TrafficMovement::getRoadLength() const
{
	return road->roadLength;
}


const uint TrafficMovement::getMaxVelocity() const
{
	return road->maxVelocity;
}


void TrafficMovement::loopTroughTimesteps(uint timeSteps)
{
	for (uint step = 1; step < timeSteps + 1; step++) {
		calculateTimestep(step);
	}
	dispResults();
}


void TrafficMovement::calculateTimestep(uint step)
{
	//if (display != nullptr)	display->dispCurrentRoad();
	//if (concWriter != nullptr) concWriter->dispConcentrations();

	VectorHelper::fillVector(*pnext, -1);
	if (concWriter != nullptr) concWriter->resetConcentrations();

	//apply rules on all cars  
	for (uint i = 0; i < road->roadLength; i++) {
		if ((*pcurrent)[i] > -1) {
			applyRules(i);
		}
	}

	calculateJunctionStep();

	calculateSourceStep();

	switchCurrentNext();

	if (display != nullptr)  display->putCurrentIntoResults(step);

	currentStep += 1;
}


void TrafficMovement::calculateSourceStep()
{
	uint sourceIndex;
	uint gap;
	for (auto &source : road->sources) {
		sourceIndex = source->getIndex();
		gap = getGapAfterOutCell(sourceIndex, road->maxVelocity);
		if (gap > 0) {
			uint speed = source->getSourceCar();
			(*pnext)[sourceIndex] = speed;
			writeConcentration(sourceIndex, speed);
		}
	}
}

void TrafficMovement::calculateJunctionStep()
{
	for (auto &junction : road->junctions) {
		junction->updateJunction();
		junction->calculateTimeStep(*this);
	}
}

void TrafficMovement::switchCurrentNext()
{
	pdummy = pcurrent;
	pcurrent = pnext;
	pnext = pdummy;
}

void TrafficMovement::applyRules(uint carIndex)
{
	uint speed = (*pcurrent)[carIndex];
	accelerateCar(speed);
	brakeCar(carIndex, speed);
	dawdleCar(carIndex, speed);
	moveCar(carIndex, speed);
}

void TrafficMovement::accelerateCar(uint & speed)
{
	if (speed <= road->maxVelocity - maxAcceleration) {
		speed += maxAcceleration;
	}
}

void TrafficMovement::brakeCar(uint carIndex, uint &speed)
{
	int neighbor = road->neighbors[carIndex];
	gap = getGapAfterCar(carIndex, speed, neighbor);
	if (speed > gap) {
		speed = gap;
	}
}

void TrafficMovement::dawdleCar(uint carIndex, uint & speed)
{
	randomNumber = distFloat(engine);

	//Barlovic / SlowToStart
	if (useSlowToStart == true && (*pcurrent)[carIndex] == 0) {
		if (randomNumber < slowStartPossibility) {
			speed = 0;
		}
		return;
	}

	//Standard NaSch
	if (randomNumber < dawdlePossibility) {
		if (speed >= maxAcceleration)
			speed -= maxAcceleration;
		else
			speed = 0;
	}
}

void TrafficMovement::moveCar(const uint carIndex, uint speed)
{
	if (speed == 0) {
		(*pnext)[carIndex] = 0;
		writeConcentration(carIndex, (*pcurrent)[carIndex]);
		return;
	}

	int neighbor = road->neighbors[carIndex];
	uint currentCell = carIndex;

	uint numberOfCellsMoved = iterateNeighborsInMove(currentCell, speed, neighbor);

	if (neighbor <= -1000 && neighbor > -2000) {
		getJunctionFromNeighbor(neighbor)->registerCar(currentCell, numberOfCellsMoved, speed, (*pcurrent)[carIndex]);
		return;
	}

	if (neighbor >= 0) {
		(*pnext)[neighbor] = speed;
		writeConcentration(neighbor, (*pcurrent)[carIndex]);
	}
}


void TrafficMovement::moveJunctionCar(uint outCellIndex, uint remainingDistance, uint speed, uint oldSpeed)
{
	if (remainingDistance == 1) {
		(*pnext)[outCellIndex] = speed;
		writeConcentration(outCellIndex, oldSpeed);
		return;
	}

	int neighbor = outCellIndex;

	uint numberOfCellsMoved = iterateNeighborsInMove(outCellIndex, remainingDistance, neighbor);

	try{
	if (neighbor <= -1000 && neighbor > -2000) {
		throw std::runtime_error("car entered two junctions in one timestep");
	}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}

	if (neighbor >= 0) {
		(*pnext)[neighbor] = speed;
		writeConcentration(neighbor, oldSpeed);
	}
}


uint TrafficMovement::iterateNeighborsInMove(uint & currentCell, uint speed, int & neighbor)
{
	uint numberOfCellsMoved = 1;

	for (uint i = 2; i <= speed; i++) {
		if (neighbor >= 0) {
			currentCell = neighbor;
			neighbor = road->neighbors[neighbor];
			++numberOfCellsMoved;
		}
		else
			break;
	}
	return numberOfCellsMoved;
}


std::shared_ptr<Junction>& TrafficMovement::getJunctionFromNeighbor(int neighbor)
{
	//calculate index in junctions vector for neighbor (-1000 to -1999)
	return road->junctions[((neighbor + 1000)*-1)];
}


std::shared_ptr<Sink>& TrafficMovement::getSinkFromNeighbor(int neighbor)
{
	//calculate index in junctions vector for neighbor (-2000 to -2999)
	int index = ((neighbor + 2000)*-1);
	return road->sinks[index];
}


uint TrafficMovement::getGapAfterCar(uint carIndex, uint speed, int neighbor)
{
	for (uint i = 0; i < (speed + road->safetyDistance); i++) {

		if (neighbor <= -2000)
			return getGapToSink(neighbor, i, speed);
		if (neighbor <= -1000)
			return getGapToJunction(neighbor, i, speed, carIndex);

		//car in Cell
		if ((*pcurrent)[neighbor] > -1)
			return adjustGapToSafetyDist(i);

		//empty cell -> get next neighbor, update currentCell
		carIndex = neighbor;
		neighbor = road->neighbors[neighbor];
	}
	return speed;
}


uint TrafficMovement::getGapAfterOutCell(uint outCellIndex, uint speed)
{
	if ((*pcurrent)[outCellIndex] > -1)
		return 0;
	 
	int neighbor = outCellIndex;

	for (uint i = 0; i < (speed + road->safetyDistance); i++) {

		if (neighbor <= -2000)
			return getGapToSink(neighbor, i, speed);
		if (neighbor <= -1000)
			return i;

		//car in Cell
		if ((*pcurrent)[neighbor] > -1)
			return adjustGapToSafetyDist(i);

		//empty cell -> get next neighbor
		neighbor = road->neighbors[neighbor];
	}
	return speed;
}


uint TrafficMovement::getGapToSink(int neighbor, uint i, uint speed)
{
	if (getSinkFromNeighbor(neighbor)->carCanEnter() && i <= speed)
		return speed;
	return i;
}


uint TrafficMovement::getGapToJunction(int neighbor, uint i, uint speed, uint currentCell)
{
	if (getJunctionFromNeighbor(neighbor)->acceptsCar(currentCell) && i <= speed)
		return speed;
	return i;
}


uint TrafficMovement::adjustGapToSafetyDist(uint gap)
{
	if (gap <= road->safetyDistance)
		return 0;
	else
		return gap - road->safetyDistance;
}

void TrafficMovement::writeConcentration(uint index, uint oldSpeed)
{
	if (concWriter != nullptr) {
		concWriter->calculateConcForSingleCar(index, oldSpeed, (*pnext)[index]);
	}
}


void TrafficMovement::writeConcentrationForJunction(uint inCellIndex, uint oldSpeed, uint speed)
{
	if (concWriter != nullptr) {
		concWriter->calculateConcForJunctionCar(inCellIndex, oldSpeed, speed);
	}
}


void TrafficMovement::dispResults() {
	if (display == nullptr)
		std::cout << "No results were saved." << std::endl;
	else
		display->dispResults(&road->neighbors, road->sinks, road->junctions, road->sources);
}


void TrafficMovement::checkCurrentForSafetyDistance()
{
	if (road->safetyDistance > 0) {
		uint neighbor;
		for (uint i = 0; i < road->roadLength; i++) {
			if ((*pcurrent)[i] > -1) {
				neighbor = road->neighbors[i];

				for (uint j = 1; j <= road->safetyDistance; j++) {
					if (neighbor <= -1000)
						break;
					if ((*pcurrent)[neighbor] > -1) {
						std::cerr << "timestep 0: safetyDistance was violated: carIndex: " << i << std::endl;
						break;
					}
					neighbor = road->neighbors[neighbor];
				}
			}
		}
	}
}


const std::vector<int> & TrafficMovement::getVehiclesForVTK()
{
	return road->currentWithLongVehicles;
}


void TrafficMovement::visualizeVehicleLengthForVTK()
{
	road->currentWithLongVehicles = *pcurrent;

	if (road->safetyDistance != 0) {
		for (uint i = 0; i < road->roadLength; i++) {
			if ((*pcurrent)[i] > -1) {
				checkSpeed((*pcurrent)[i]);
				int neighbor = road->neighbors[i];

				for (uint j = 1; j <= road->safetyDistance; j++) {

					if (neighbor <= -1000)
						break;
					if ((*pcurrent)[neighbor] > -1) {
						std::cerr << "safetyDistance was violated: timestep: " << currentStep << "\t carIndex: " << i << std::endl;
						break;
					}
					else
						(road->currentWithLongVehicles)[neighbor] = (*pcurrent)[i];
					neighbor = road->neighbors[neighbor];
				}
			}
		}
	}
}


void TrafficMovement::checkSpeed(uint speed)
{
	if (speed > road->maxVelocity) {
		std::cerr << "Speed is greater than allowed maxSpeed" << std::endl;
	}
}