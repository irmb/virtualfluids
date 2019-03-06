#include "JunctionRandom.h"



JunctionRandom::JunctionRandom(const std::vector<uint> &inCellIndices, const std::vector<uint> &outCellIndices)
{
	data.inCellIndices = inCellIndices;
	data.outCellIndices = outCellIndices;

	uint inRoads = inCellIndices.size();

	data.carCanEnter.resize(inRoads);
	std::fill(data.carCanEnter.begin(), data.carCanEnter.end(), true);

	data.carsOnJunction.resize(inRoads);
	std::fill(data.carsOnJunction.begin(), data.carsOnJunction.end(), -1);

	data.alreadyMoved.resize(inRoads);
	std::fill(data.alreadyMoved.begin(), data.alreadyMoved.end(), 0);
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


void JunctionRandom::registerCar(uint cellIndex, uint numberOfCellsAlreadyMoved, uint speed)
{
	uint index = getInCellsVectorIndex(cellIndex);
	data.carsOnJunction[index] = speed;
	data.carCanEnter[index] = false;
	data.alreadyMoved[index] = numberOfCellsAlreadyMoved;
}


uint JunctionRandom::getInCellsVectorIndex(uint cellIndex)
{
	try {
		auto it = find(data.inCellIndices.begin(), data.inCellIndices.end(), cellIndex);

		if (it != data.inCellIndices.end())
		{
			return distance(data.inCellIndices.begin(), it);
		}

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
	index = 0;
	for (int carSpeed : data.carsOnJunction) {

		if (carSpeed >= 0) { //check if there is a car on the junction

			if (carSpeed == 0)
				std::cerr << "speed on junction was 0" << std::endl;
			else
				applyRules(carSpeed, index, road);
			data.alreadyMoved[index] = 0;

		}
		++index;
	}
}


void JunctionRandom::applyRules(int & carSpeed, const int & index, TrafficMovement& road)
{
	remainingDistance = static_cast<uint>(carSpeed) - data.alreadyMoved[index];

	if (remainingDistance > 0) {
		outCell = chooseOutCell(index);
		if (outCell >= 0) {
			breakCar(outCell, carSpeed, remainingDistance, road);

			if (remainingDistance > 0)
				moveCar(outCell, carSpeed, index, road);
			else
				data.carsOnJunction[index] = 1;
			return;
		}
	}
	data.carsOnJunction[index] = 1;
}


void JunctionRandom::breakCar(uint outCellIndex, int &speed, uint &remainingDistance, TrafficMovement& road)
{
	gap = road.getGapAfterOutCell(outCellIndex, remainingDistance);
	if (gap < remainingDistance) {
		speed = speed - remainingDistance + gap;
		remainingDistance = gap;
	}
}


void JunctionRandom::moveCar(uint outCell, int & carSpeed, const int & index, TrafficMovement& road)
{
	road.moveJunctionCar(outCell, remainingDistance, carSpeed);
	data.carsOnJunction[index] = -1;
	data.carCanEnter[index] = true;
}


int JunctionRandom::chooseOutCell(const int & index)
{
	std::vector<uint> outCellsTemp;

	if (data.carCanNotEnterThisOutCell.size() > 0 && data.carCanNotEnterThisOutCell[index]) {
		for (uint cell : data.possibleOutCells) {
			if (cell != data.carCanNotEnterThisOutCell[index])
				outCellsTemp.push_back(cell);
		}
	}
	else
		outCellsTemp = data.possibleOutCells;


	switch (outCellsTemp.size()) {
	case 0:
		return -1;
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
	}

	outCell = outCellsTemp[random];
	data.possibleOutCells.erase(std::remove(data.possibleOutCells.begin(), data.possibleOutCells.end(), outCell), data.possibleOutCells.end());
	return outCell;
}


const std::vector<uint>& JunctionRandom::getInCellIndices() const
{
	return data.inCellIndices;
}

void JunctionRandom::dispJunction(const uint index, uint roadLength) const
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

const uint JunctionRandom::getNumCarsOnJunction() const
{
	uint num = 0;
	for (auto car : data.carsOnJunction) {
		if (car >= 0)
			++num;
	}
	return num;
}


void JunctionRandom::checkOutCellIndices(uint roadLength)
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