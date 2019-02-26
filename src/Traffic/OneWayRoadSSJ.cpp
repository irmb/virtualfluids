#include "OneWayRoadSSJ.h"

OneWayRoadSSJ::OneWayRoadSSJ(unique_ptr<RoadNetworkData> road, const float dawdlePossibility) 
{
	this->road = move(road);

	pcurrent = &this->road->current;
	pnext = &(this->road->next);

	initDawdle(dawdlePossibility);
	initResults();
}

OneWayRoadSSJ::OneWayRoadSSJ()
{
}

OneWayRoadSSJ::~OneWayRoadSSJ()
{
}


void OneWayRoadSSJ::calculateTimestep(unsigned int step)
{
//	dispCurrentRoad();
	writeConcentration();

	VectorHelper::fillVector(*pnext, -1);
	vector<unsigned int> cars = findCarIndicesInCurrent();

	for (auto &car : cars)
		applyRules(car);

	for (auto &junction : road->junctions) {
		junction->updateJunction();
		junction->calculateTimeStep(*this);
	}

	calculateSourceStep();

	switchCurrentNext();

	putCurrentIntoResults(step);
	currentStep += 1;
}


void OneWayRoadSSJ::calculateSourceStep()
{
	unsigned int sourceIndex;
	unsigned int maxSpeed;
	for (auto &source : road->sources) {
		sourceIndex = source->getIndex();
		maxSpeed = getGapAfterOutCell(sourceIndex, road->maxVelocity);
		if (maxSpeed > 0)
			(*pnext)[sourceIndex] = source->getSourceCar();
	}
}


void OneWayRoadSSJ::moveCar(unsigned int carIndex, unsigned int speed)
{
	if (speed == 0) {
		(*pnext)[carIndex] = 0;
		return;
	}

	int neighbor = road->neighbors[carIndex];
	unsigned int currentCell = carIndex;

	unsigned int numberOfCellsMoved = iterateNeighborsInMove(currentCell, speed, neighbor);

	if (neighbor <= -1000 && neighbor > -2000) {
		getJunctionFromNeighbor(neighbor)->registerCar(currentCell, numberOfCellsMoved, speed); 
		return;
	}

	if (neighbor >= 0)
		(*pnext)[neighbor] = speed;
}

unsigned int OneWayRoadSSJ::iterateNeighborsInMove(unsigned int &currentCell, unsigned int speed, int &neighbor)
{
	unsigned int numberOfCellsMoved = 1;

	for (unsigned int i = 2; i <= speed; i++) {
		if (neighbor >= 0) {
			currentCell = neighbor;
			numberOfCellsMoved += 1;
			neighbor = road->neighbors[neighbor];
		}
		else
			break;
	}

	return numberOfCellsMoved;
}

unsigned int OneWayRoadSSJ::getGapAfterCar(unsigned int carIndex, unsigned int speed, int neighbor)
{
	unsigned int currentCell = carIndex;

	for (unsigned int i = 0; i < (speed + road->safetyDistance); i++) {

		if (neighbor <= -2000)
			return getGapToSink(neighbor, i, speed);
		if (neighbor <= -1000)
			return getGapToJunction(neighbor, i, speed, currentCell);

		//car in Cell
		if ((*pcurrent)[neighbor] > -1)
			return adjustGapConsideringSafetyDistance(i);

		//empty cell -> get next neighbor, update currentCell
		currentCell = neighbor;
		neighbor = road->neighbors[neighbor];
	}
	return speed;
}



unsigned int OneWayRoadSSJ::getGapToJunction(int neighbor, unsigned int i, unsigned int speed, unsigned int currentCell)
{
	if (getJunctionFromNeighbor(neighbor)->acceptsCar(currentCell) && i <= speed)
		return speed;
	return i;
}

unsigned int OneWayRoadSSJ::getGapToSink(int neighbor, unsigned int i, unsigned int speed)
{
	if (getSinkFromNeighbor(neighbor)->carCanEnter())
		return speed;
	return adjustGapConsideringSafetyDistance(i);
}


unsigned int OneWayRoadSSJ::getGapAfterOutCell(unsigned int outCellIndex, unsigned int speed)
{
	if ((*pcurrent)[outCellIndex] > -1)
		return 0;
	gap=getGapAfterCar(outCellIndex, speed, outCellIndex);
	return  gap;
}


void OneWayRoadSSJ::moveJunctionCar(unsigned int outCellIndex, unsigned int remainingDistance, unsigned int speed)
{
	if (remainingDistance == 1) {
		(*pnext)[outCellIndex] = speed;
		return;
	}

	int neighbor = outCellIndex;
	unsigned int currentCell = outCellIndex;

	unsigned int numberOfCellsMoved = iterateNeighborsInMove(currentCell, remainingDistance, neighbor);

	if (neighbor <= -1000 && neighbor > -2000) {
		(*pnext)[currentCell] = speed; // car can't enter more than one junction in one timestep 
		return;
	}

	if (neighbor >= 0)
		(*pnext)[neighbor] = speed;
}


shared_ptr<Junction>&  OneWayRoadSSJ::getJunctionFromNeighbor(int neighbor)
{
	//calculate index in junctions vector for neighbor (-1000 to -1999)
	return road->junctions[((neighbor + 1000)*-1)];
}

shared_ptr<Sink>& OneWayRoadSSJ::getSinkFromNeighbor(int neighbor)
{
	//calculate index in junctions vector for neighbor (-2000 to -2999)
	int index = ((neighbor + 2000)*-1);
	return road->sinks[index];
}


const unsigned int OneWayRoadSSJ::getNumberOfCars() const
{
	unsigned int num = 0;
	for (auto& junc : road->junctions) {
		num += junc->getNumCarsOnJunction();
	}

	return OneWayRoad::getNumberOfCars() + num;
}


void OneWayRoadSSJ::dispResults()
{
	if (!saveResults) {
		std::cout << "No results were saved" << std::endl;
		return;
	}

	visualizeSafetyDistanceForConsole();
	reverse(results.begin(), results.end());

	for (unsigned int i = 0; i < results.size(); i++) {

		dispJunctionsAtCell(i);
		dispSinksAtCell(i);
		dispSourcesAtCell(i);

		for (unsigned int j = 0; j < results[i].size(); j++) {
			VectorHelper::makeVectorOutputColourful(results[i][j]);
			cout << setw(4) << results[i][j];
		}

		cout << endl;
	}
	cout << endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // set output default white 7;
}

void OneWayRoadSSJ::visualizeSafetyDistanceForConsole() {
	if (road->safetyDistance != 0) {
		int neighbor;
		for (unsigned int step = 0; step <= timeSteps; step++) {
			for (unsigned int i = 0; i < road->roadLength; i++) {
				if (results[i][step] > -1) {
					neighbor = road->neighbors[i];
					for (unsigned int j = 1; j <= road->safetyDistance; j++) {
						//junction or sink
						if (neighbor <= -1000) {
							break;
						}
						if (results[neighbor][step] > -1) {
							cerr << "safetyDistance was violated: timestep: " << step << "\t carIndex: " << i << endl;
							break;
						}
						else
							results[neighbor][step] = -5;
						neighbor = road->neighbors[neighbor];
					}
				}
			}
		}
	}

}

void OneWayRoadSSJ::dispJunctionsAtCell(unsigned int index)  const
{
	for (auto& junc : road->junctions) {
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // set output default white 7;
		junc->dispJunction(index, road->roadLength);
	}
}

void OneWayRoadSSJ::dispSinksAtCell(unsigned int index)  const
{
	for (auto& sink : road->sinks) {
		if (sink->getIndex() == road->roadLength - index - 1) {
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12); //set output bright green 10, bright red 12
			cout << setw(4) << 1 - (sink->getPossibilityBeingBlocked());
			return;
		}
		cout << setw(4) << " ";
	}
}


void OneWayRoadSSJ::dispSourcesAtCell(unsigned int index)  const
{
	for (auto& source : road->sources) {
		if (source->getIndex() == road->roadLength - index - 1) {
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 10); //set output bright green 10, bright red 12
			cout << setw(4) << source->getPossibility();
			return;
		}
		cout << setw(4) << " ";
	}
}