#include "JunctionRandom.h"

#include <iostream>
#include <algorithm>

#include "TrafficMovement.h"

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Utilities/safe_casting.h"


JunctionRandom::JunctionRandom(const std::vector<uint> &inCellIndices, const std::vector<uint> &outCellIndices)
{
	data.inCellIndices = inCellIndices;
	data.outCellIndices = outCellIndices;

	uint inRoads = castSizeT_Uint(inCellIndices.size());

	data.carCanEnter.resize(inRoads);
	std::fill(data.carCanEnter.begin(), data.carCanEnter.end(), true);

	data.carsOnJunction.resize(inRoads);
	std::fill(data.carsOnJunction.begin(), data.carsOnJunction.end(), -1);

	data.alreadyMoved.resize(inRoads);
	std::fill(data.alreadyMoved.begin(), data.alreadyMoved.end(), 0);

	data.oldSpeeds.resize(inRoads);
}


void JunctionRandom::setCellIndexForNoUTurn(std::vector<int> carCanNotEnterThisOutCell)
{
	try {

		if (data.inCellIndices.size() != carCanNotEnterThisOutCell.size()) throw invalidInput_error("The Vector carCanNotEnterThisOutCell and inCellIndices have to be the same size.");
		data.carCanNotEnterThisOutCell = carCanNotEnterThisOutCell;

	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


bool JunctionRandom::acceptsCar(uint cellIndex)
{
	return data.carCanEnter[getInCellsVectorIndex(cellIndex)];
}


void JunctionRandom::registerCar(uint cellIndex, uint alreadyMoved, uint speed, uint oldSpeed)
{
	uint index = getInCellsVectorIndex(cellIndex);

	data.carsOnJunction[index] = speed - 1; //all cars, which enter the junction have to slow down by one increment
	//data.carsOnJunction[index] = 0; //all cars, which enter the junction have to stop
	data.oldSpeeds[index] = oldSpeed;
	data.carCanEnter[index] = false;
	data.alreadyMoved[index] = alreadyMoved;
}


uint JunctionRandom::getInCellsVectorIndex(uint cellIndex)
{
	try {
		auto it = find(data.inCellIndices.begin(), data.inCellIndices.end(), cellIndex);

		if (it != data.inCellIndices.end())
			return static_cast <uint> (distance(data.inCellIndices.begin(), it));

		throw std::runtime_error("The passed cell is not an incoming cell to this junction.");
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


void JunctionRandom::updateJunction()
{
	data.possibleOutCells = data.outCellIndices;
}


void JunctionRandom::calculateTimeStep(TrafficMovement& road)
{
	uint index = 0;
	for (int carSpeed : data.carsOnJunction) {
		if (carSpeed >= 0) { //check if there is a car on the junction
			applyRules(carSpeed, index, road);
			data.alreadyMoved[index] = 0;
		}
		++index;
	}
	writeConcentrations(road);
}


void JunctionRandom::applyRules(int & carSpeed, int index, TrafficMovement& road)
{
	if (carSpeed == 0 && data.alreadyMoved[index] == 0)
		carSpeed += 1;

	int remainingDistance = static_cast<uint>(carSpeed) - data.alreadyMoved[index];
	if (remainingDistance > 0) {
		int outCell = chooseOutCell(index);
		if (outCell >= 0) {
			brakeCar(outCell, carSpeed, remainingDistance, road);
			if (remainingDistance > 0) {
				moveCar(outCell, carSpeed, index, remainingDistance, road);
				return;
			}
		}
	}
	data.carsOnJunction[index] = 0;				//cars, which can't cross the junctionin one timestep, because they already moved to many cells, loose their speed.
	//data.getCarsOnJunction[index] = carSpeed;	//cars, which can't cross the junction in one timestep, because they already moved to many cells, keep their speed.
}


void JunctionRandom::brakeCar(uint outCellIndex, int &speed, int &remainingDistance, TrafficMovement& road)
{
	int gap = road.getGapAfterOutCell(outCellIndex, remainingDistance);
	if (gap < remainingDistance) {
		if (gap > speed) gap = speed;
		speed = speed - remainingDistance + gap;
		remainingDistance = gap;
	}
}


void JunctionRandom::moveCar(uint outCell, int carSpeed, int index, int remainingDistance, TrafficMovement& road)
{
	road.moveJunctionCar(outCell, remainingDistance, carSpeed, data.oldSpeeds[index]);
	data.carsOnJunction[index] = -1;
	data.carCanEnter[index] = true;
}


int JunctionRandom::chooseOutCell(int index)
{
	std::vector<uint> outCellsTemp;

	if (data.carCanNotEnterThisOutCell.size() > 0 && data.carCanNotEnterThisOutCell[index] >= 0) {
		for (uint cell : data.possibleOutCells) {
			if (cell != data.carCanNotEnterThisOutCell[index])
				outCellsTemp.push_back(cell);
		}
	}
	else
		outCellsTemp = data.possibleOutCells;

	int random = generateRandomOutCellIndex(castSizeT_Uint(outCellsTemp.size()));
	if (random < 0) return -1;

	int outCell = outCellsTemp[random];
	data.possibleOutCells.erase(std::remove(data.possibleOutCells.begin(), data.possibleOutCells.end(), outCell), data.possibleOutCells.end());
	return outCell;
}


int JunctionRandom::generateRandomOutCellIndex(uint outCellsTempSize)
{
	int random = 1000;

	switch (outCellsTempSize) {
	case 0:
		random = -1;
		break;
	case 1:
		random = 0;
		break;
	case 2:
		random = data.distInt2(data.engine);
		break;
	case 3:
		random = data.distInt3(data.engine);
		break;
	case 4:
		random = data.distInt4(data.engine);
		break;
	default:
		std::cerr << "default in JunctionRandom::chooseOutCell()" << std::endl;
		break;
	}

	return random;
}

void JunctionRandom::writeConcentrations(TrafficMovement& road)
{
	int i = 0;
	for (int carSpeed : data.carsOnJunction) {
		if (carSpeed >= 0) {
			road.writeConcentrationForJunction(data.inCellIndices[i], data.oldSpeeds[i], data.carsOnJunction[i]);
			data.oldSpeeds[i] = data.carsOnJunction[i];
		}
		++i;
	}
}

const std::vector<uint>& JunctionRandom::getInCellIndices() const
{
	return data.inCellIndices;
}

const std::vector<uint>& JunctionRandom::getOutCellIndices() const
{
	return data.outCellIndices;
}

const std::vector<bool>& JunctionRandom::getCarCanEnter() const
{
	return data.carCanEnter;
}

const std::vector<int>& JunctionRandom::getCarsOnJunction() const
{
	return data.carsOnJunction;
}

const std::vector<uint>& JunctionRandom::getAlreadyMoved() const
{
	return data.alreadyMoved;
}

const std::vector<uint>& JunctionRandom::getOldSpeeds() const
{
	return data.oldSpeeds;
}

const std::vector<int>& JunctionRandom::getCarCanNotEnterThisOutCell() const
{
	return data.carCanNotEnterThisOutCell;
}

void JunctionRandom::dispJunction(const uint index, const uint roadLength) const
{
	if (find(data.inCellIndices.begin(), data.inCellIndices.end(), (roadLength - index - 1)) != data.inCellIndices.end()) {
		std::cout << std::setw(4) << "in";
	}
	else if (find(data.outCellIndices.begin(), data.outCellIndices.end(), (roadLength - index - 1)) != data.outCellIndices.end()) {
		std::cout << std::setw(4) << "out";
	}
	else {
		std::cout << std::setw(4) << " ";
	}
}

uint JunctionRandom::getNumCarsOnJunction() const
{
	uint num = 0;
	for (auto car : data.carsOnJunction)
		if (car >= 0)
			++num;
	return num;
}

void JunctionRandom::checkOutCellIndices(const uint roadLength) const
{
	try {
		for (uint cell : data.outCellIndices)
			if (cell >= roadLength) throw invalidInput_error("The indices of incoming cells to a junction are greater than the roadLength.");
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


