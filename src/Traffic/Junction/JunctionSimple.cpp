#include "JunctionSimple.h"



JunctionSimple::JunctionSimple(const vector<unsigned int> &inCellIndices, const vector<unsigned int> &outCellIndices)
{
	data.inCellIndices = inCellIndices;
	data.outCellIndices = outCellIndices;

	unsigned int inRoads = inCellIndices.size();

	data.carCanEnter.resize(inRoads);
	fill(data.carCanEnter.begin(), data.carCanEnter.end(), true);

	data.carsOnJunction.resize(inRoads);
	fill(data.carsOnJunction.begin(), data.carsOnJunction.end(), -1);

	data.alreadyMoved.resize(inRoads);
	fill(data.alreadyMoved.begin(), data.alreadyMoved.end(), 0);
}


void JunctionSimple::setCellIndexForNoUTurn(std::vector<int> carCanNotEnterThisOutCell)
{
	try {

		if (data.inCellIndices.size() != carCanNotEnterThisOutCell.size()) throw invalidInput_error("The Vector carCanNotEnterThisOutCell and inCellIndices have to be the same size.");
		data.carCanNotEnterThisOutCell = carCanNotEnterThisOutCell;

	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}


bool JunctionSimple::acceptsCar(unsigned int cellIndex)
{
	return data.carCanEnter[getInCellsVectorIndex(cellIndex)];
}


void JunctionSimple::registerCar(unsigned int cellIndex, unsigned int numberOfCellsAlreadyMoved, unsigned int speed)
{
	unsigned int index = getInCellsVectorIndex(cellIndex);
	data.carsOnJunction[index] = speed;
	data.carCanEnter[index] = false;
	data.alreadyMoved[index] = numberOfCellsAlreadyMoved;
}


unsigned int JunctionSimple::getInCellsVectorIndex(unsigned int cellIndex)
{
	try {
		auto it = find(data.inCellIndices.begin(), data.inCellIndices.end(), cellIndex);

		if (it != data.inCellIndices.end())
		{
			return distance(data.inCellIndices.begin(), it);
		}

		throw runtime_error("The passed cell is not an incoming cell to this junction.");
	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}


void JunctionSimple::updateJunction()
{
	data.freeOutCells = data.outCellIndices;
}


void JunctionSimple::calculateTimeStep(OneWayRoadSSJ& road)
{
	index = 0;
	for (int carSpeed : data.carsOnJunction) {

		if (carSpeed >= 0) { //check if there is a car on the junction

			if (carSpeed == 0)
				cerr << "speed on junction was 0" << endl;
			else
				applyRules(carSpeed, index, road);
			data.alreadyMoved[index] = 0;

		}
		++index;
	}
}


void JunctionSimple::applyRules(int & carSpeed, const int & index, OneWayRoadSSJ& road)
{
	remainingDistance = static_cast<unsigned int>(carSpeed) - data.alreadyMoved[index];

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


void JunctionSimple::breakCar(unsigned int outCellIndex, int &speed, unsigned int &remainingDistance, OneWayRoadSSJ& road)
{
	gap = road.getGapAfterOutCell(outCellIndex, remainingDistance);
	if (gap < remainingDistance) {
		speed = speed - remainingDistance + gap;
		remainingDistance = gap;
	}
}


void JunctionSimple::moveCar(unsigned int outCell, int & carSpeed, const int & index, OneWayRoadSSJ& road)
{
	road.moveJunctionCar(outCell, remainingDistance, carSpeed);
	data.carsOnJunction[index] = -1;
	data.carCanEnter[index] = true;
}


int JunctionSimple::chooseOutCell(const int & index)
{
	vector<unsigned int> freeOutCellsTemp;

	if (data.carCanNotEnterThisOutCell.size() > 0 && data.carCanNotEnterThisOutCell[index]) {
		for (unsigned int cell : data.freeOutCells) {
			if (cell != data.carCanNotEnterThisOutCell[index])
				freeOutCellsTemp.push_back(cell);
		}
	}
	else
		freeOutCellsTemp = data.freeOutCells;


	switch (freeOutCellsTemp.size()) {
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

	outCell = freeOutCellsTemp[random];
	data.freeOutCells.erase(std::remove(data.freeOutCells.begin(), data.freeOutCells.end(), outCell), data.freeOutCells.end());
	return outCell;
}


const vector<unsigned int>& JunctionSimple::getInCellIndices() const
{
	return data.inCellIndices;
}

void JunctionSimple::dispJunction(const unsigned int index, unsigned int roadLength) const
{
	if (find(data.inCellIndices.begin(), data.inCellIndices.end(), (roadLength - index - 1)) != data.inCellIndices.end()) {
		cout << setw(4) << "in";
	}
	else if (find(data.outCellIndices.begin(), data.outCellIndices.end(), (roadLength - index - 1)) != data.outCellIndices.end()) {
		cout << setw(4) << "out";
	}
	else {
		cout << setw(4) << " ";
	}
}

const unsigned int JunctionSimple::getNumCarsOnJunction() const
{
	unsigned int num = 0;
	for (auto car : data.carsOnJunction) {
		if (car >= 0)
			++num;
	}
	return num;
}


void JunctionSimple::checkOutCellIndices(unsigned int roadLength)
{
	try {
		for (unsigned int cell : data.outCellIndices)
			if (cell >= roadLength) throw invalidInput_error("The indices of incoming cells to a junction are greater than the roadLength.");
	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}