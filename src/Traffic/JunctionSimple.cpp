#include "JunctionSimple.h"



JunctionSimple::JunctionSimple(const vector<unsigned int> &inCellIndices, const vector<unsigned int> &outCellIndices) :
	engine{ random_device{}() }
{
	data.inCellIndices = inCellIndices;
	data.outCellIndices = outCellIndices;

	unsigned int inRoads = inCellIndices.size();

	carCanEnter.resize(inRoads);
	fill(carCanEnter.begin(), carCanEnter.end(), true);

	carsOnJunction.resize(inRoads);
	fill(carsOnJunction.begin(), carsOnJunction.end(), -1);

	alreadyMoved.resize(inRoads);
	fill(alreadyMoved.begin(), alreadyMoved.end(), 0);
}


JunctionSimple::~JunctionSimple()
{
}


bool JunctionSimple::acceptsCar(unsigned int cellIndex)
{
	return carCanEnter[getInCellsVectorIndex(cellIndex)];
}


void JunctionSimple::registerCar(unsigned int cellIndex, unsigned int numberOfCellsAlreadyMoved, unsigned int speed)
{
	unsigned int index = getInCellsVectorIndex(cellIndex);
	carsOnJunction[index] = speed;
	carCanEnter[index] = false;
	alreadyMoved[index] = numberOfCellsAlreadyMoved;
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
	freeOutCells = data.outCellIndices;
}


void JunctionSimple::calculateTimeStep(OneWayRoadSSJ& road)
{
	index = 0;
	for (int carSpeed : carsOnJunction) {

		if (carSpeed >= 0) { //check if there is a car on the junction

			if (carSpeed == 0)
				cerr << "speed on junction was 0" << endl;
			else
				applyRules(carSpeed, index, road);
			alreadyMoved[index] = 0;

		}
		++ index;
	}
}


void JunctionSimple::applyRules(int & carSpeed,const int & index, OneWayRoadSSJ& road)
{
	remainingDistance = static_cast<unsigned int>(carSpeed) - alreadyMoved[index];

	if (remainingDistance > 0 && freeOutCells.size() > 0) {
		outCell = chooseOutCell();
		breakCar(outCell, carSpeed, remainingDistance, road);

		if (remainingDistance > 0)
			moveCar(carSpeed, index, road);
		else
			carsOnJunction[index] = 1;
	}
	else
		carsOnJunction[index] = 1;
}


void JunctionSimple::breakCar(unsigned int outCellIndex, int &speed, unsigned int &remainingDistance, OneWayRoadSSJ& road)
{
	gap = road.getGapAfterOutCell(outCellIndex, remainingDistance);
	if (gap < remainingDistance) {
		speed = speed - remainingDistance + gap;
		remainingDistance = gap;
	}

}


void JunctionSimple::moveCar(int & carSpeed, const int & index, OneWayRoadSSJ& road)
{
	road.moveJunctionCar(outCell, remainingDistance, carSpeed);
	carsOnJunction[index] = -1;
	carCanEnter[index] = true;
}


unsigned int JunctionSimple::chooseOutCell()
{
	uniform_int_distribution<unsigned int> distFloat(0, freeOutCells.size() - 1);
	random = distFloat(engine);

	outCell = freeOutCells[random];
	freeOutCells.erase(freeOutCells.begin() + random);
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
	for (auto car : carsOnJunction) {
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