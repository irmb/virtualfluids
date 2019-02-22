#include "OneWayRoad.h"
#include <stdexcept>

OneWayRoad::OneWayRoad()
{
}


OneWayRoad::OneWayRoad(unique_ptr<RoadNetworkData> road, const float dawdlePossibility)
{
	this->road = move(road);

	pcurrent = &this->road->current;
	pnext = &(this->road->next);

	initDawdle(dawdlePossibility);
	initResults();
}


OneWayRoad::~OneWayRoad()
{
	pcurrent = NULL;
	pnext = NULL;
	pdummy = NULL;
}


void OneWayRoad::initResults()
{
	if (saveResults) {
		results.resize(road->roadLength, vector<int>(1));
		VectorHelper::fillVector(results, -10);
	}
}


void OneWayRoad::initDawdle(const float dawdlePossibility)
{
	try {
		if (dawdlePossibility >= 0 && dawdlePossibility < 1) {
			this->dawdlePossibility = dawdlePossibility;
		}
		else {
			throw invalidInput_error("The dawdlePossibility should be between 0 and 1.");
		}
	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}


void OneWayRoad::setSlowToStart(const float slowStartPossibility)
{
	try {
		if (slowStartPossibility >= 0 && slowStartPossibility < 1) {

			this->slowStartPossibility = slowStartPossibility;
			useSlowToStart = true;

		}
		else {
			throw invalidInput_error("The slowStartPossibility should be between 0 and 1.");
		}
	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	}
}


void OneWayRoad::setConcentrationOutwriter(unique_ptr<ConcentrationOutwriter> writer)
{
	this->concWriter = move(writer);
}


void OneWayRoad::setSaveResults(bool saveResults)
{
	saveResults = true;
}


const unsigned int OneWayRoad::getNumberOfCars() const
{
	return findCarIndicesInCurrent().size();
}

const int OneWayRoad::getSpeedAtPosition(unsigned int pos) const
{
	return (*pcurrent)[pos];
}

unsigned int OneWayRoad::getRoadLength() const
{
	return road->roadLength;
}

unsigned int OneWayRoad::getMaxVelocity() const
{
	return road->maxVelocity;
}

void OneWayRoad::initCalculation(int timeSteps)
{
	this->timeSteps = timeSteps;

	initResultsForCalculation();

	putCurrentIntoResults(0);
}

void OneWayRoad::initResultsForCalculation()
{
	if (!saveResults) return;

	for (unsigned int i = 0; i < road->roadLength; i++) {
		results[i].resize(timeSteps + 1);
	}
	VectorHelper::fillVector(results, -1);
}


void OneWayRoad::loopTroughTimesteps()
{
	for (int step = 1; step < timeSteps + 1; step++) {
		calculateTimestep(step);
	}
}


void OneWayRoad::calculateTimestep(unsigned int step)
{
	//	dispCurrentRoad();
	writeConcentration();

	VectorHelper::fillVector(*pnext, -1);
	vector<unsigned int> cars = findCarIndicesInCurrent();

	for (auto &car : cars) {
		applyRules(car);
	}

	switchCurrentNext();
	putCurrentIntoResults(step);
	currentStep += 1;
}

void OneWayRoad::switchCurrentNext()
{
	pdummy = pcurrent;
	pcurrent = pnext;
	pnext = pdummy;
}

void OneWayRoad::applyRules(unsigned int carIndex)
{
	unsigned int speed = (*pcurrent)[carIndex];
	accelerateCar(speed);
	breakCar(carIndex, speed);
	dawdleCar(carIndex, speed);
	moveCar(carIndex, speed);
}

void OneWayRoad::accelerateCar(unsigned int & speed)
{
	if (speed < road->maxVelocity) {
		speed += 1;
	}
}

void OneWayRoad::breakCar(unsigned int carIndex, unsigned int &speed)
{
	int neighbor = road->neighbors[carIndex];
	gap = getGapAfterCar(carIndex, speed, neighbor);
	if (speed > gap) {
		speed = gap;
	}
}

void OneWayRoad::dawdleCar(unsigned int carIndex, unsigned int & speed)
{
	randomNumber = distFloat(engine);

	//Barlovic / SlowToStart
	if ((*pcurrent)[carIndex] == 0 && useSlowToStart == true) {
		if (randomNumber < slowStartPossibility) {
			speed = 0;
		}
		return;
	}

	//Standard NaSch
	if (randomNumber < dawdlePossibility) {
		if ((speed) > 0) {
			speed -= 1;
		}
		else {
			speed = 0;
		}
	}
}

void OneWayRoad::moveCar(unsigned int carIndex, unsigned int speed)
{
	int neighbor = road->neighbors[carIndex];

	if (speed == 0) {
		(*pnext)[carIndex] = 0;
	}
	else {
		for (unsigned int i = 2; i <= speed; i++) {
			neighbor = road->neighbors[neighbor];
		}
		(*pnext)[neighbor] = speed;
	}
}


unsigned int OneWayRoad::getGapAfterCar(unsigned int carIndex, unsigned int speed, int neighbor)
{
	for (unsigned int i = 0; i < (speed + road->safetyDistance); i++) {
		if ((*pcurrent)[neighbor] > -1) {
			return adjustGapConsideringSafetyDistance(i);
		}
		neighbor = road->neighbors[neighbor];
	}
	return speed;
}


unsigned int OneWayRoad::adjustGapConsideringSafetyDistance(unsigned int speed)
{
	if (speed <= road->safetyDistance)
		return 0;
	else
		return speed - road->safetyDistance;
}


vector<unsigned int> OneWayRoad::findCarIndicesInCurrent() const
{
	vector<unsigned int> cars;
	for (unsigned int i = 0; i < road->roadLength; i++) {
		if ((*pcurrent)[i] > -1) {
			cars.push_back(i);
		}
	}
	return cars;
}

void OneWayRoad::putCurrentIntoResults(unsigned int step)
{
	if (!saveResults) return;

	for (unsigned int i = 0; i < road->roadLength; i++) {
		results[i][step] = (*pcurrent)[i];
	}
}

void OneWayRoad::writeConcentration()
{
	if (concWriter != nullptr) {
		concWriter->writeToArray(*pcurrent);
	}
}


void OneWayRoad::writeResultsToFile() const
{
	try {


		fstream outFile("results.txt", fstream::out | fstream::trunc);
		if (outFile.is_open())
		{
			for (unsigned int i = 0; i < results.size(); i++) {
				for (unsigned int j = 0; j < results[i].size() - 1; j++)
					outFile << results[i][j] << " ";

				outFile << results[i][results[i].size() - 1];
				outFile << endl;
			}
			cout << "Finished writing data to file" << endl;
		}


		else
			throw runtime_error("Couldn't open file");

	}
	catch (const exception& e) {
		cerr << e.what() << endl;
	}
	catch (...) {
		cerr << "unknown exception while writing to file" << endl;
	}
}



void OneWayRoad::dispCurrentRoad() const
{
	cout << "current: ( step: " << currentStep << " )" << endl;
	VectorHelper::dispVectorColour(*pcurrent);
}

void OneWayRoad::dispResults()
{
	if (!saveResults) {
		std::cout << "No results were saved" << std::endl;
		return;
	}

	cout << "results:" << endl;

	visualizeSafetyDistanceForConsole();
	reverse(results.begin(), results.end());
	VectorHelper::dispVectorColour(results);
}

void OneWayRoad::visualizeSafetyDistanceForConsole()
{
	if (road->safetyDistance != 0) {
		for (unsigned int step = 0; step <= timeSteps; step++) {
			for (unsigned int i = 0; i < road->roadLength; i++) {
				if (results[i][step] > -1) {
					int neighbor = road->neighbors[i];
					for (unsigned int j = 1; j <= road->safetyDistance; j++) {
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


const std::vector<int>& OneWayRoad::getVehiclesForVTK()
{
	return road->currentWithLongVehicles;
}

void OneWayRoad::visualizeVehicleLengthForVTK()
{
	road->currentWithLongVehicles = *pcurrent;

	if (road->safetyDistance != 0) {
		for (unsigned int i = 0; i < road->roadLength; i++) {
			if ((*pcurrent)[i] > -1) {
				int neighbor = road->neighbors[i];
				for (unsigned int j = 1; j <= road->safetyDistance; j++) {
					if ((*pcurrent)[neighbor] > -1) {
						cerr << "safetyDistance was violated: timestep: " << currentStep << "\t carIndex: " << i << endl;
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