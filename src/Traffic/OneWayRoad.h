#pragma once
#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <random> 
#include <memory>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output

#include "invalidInput_error.h"
#include "VectorHelper.h"
#include "RandomHelper.h"
#include "RoadNetworkData.h"
#include "ConcentrationOutwriter.h"

using namespace std;

class VF_PUBLIC OneWayRoad //Periodic BC
{
protected:
	shared_ptr<RoadNetworkData> road;

	float dawdlePossibility;
	unsigned int timeSteps;
	unsigned int currentStep = 0;			//only needed for dispCurrent()

	bool useSlowToStart = false;
	float slowStartPossibility;

	mt19937 engine = RandomHelper::make_engine();
	uniform_real_distribution<float> distFloat{ 0.0, 1.0 };

	vector<int> *pcurrent;
	vector<int> *pnext;
	vector<int> *pdummy;

	vector<vector<int> > results;		//saves the results of the calculation; x-axis = timesteps, y axis = positions

	unique_ptr<ConcentrationOutwriter> concWriter = nullptr;

public:
	OneWayRoad(shared_ptr<RoadNetworkData> road, const float dawdlePossibility);
	OneWayRoad();
	OneWayRoad(const OneWayRoad&) = delete;
	virtual ~OneWayRoad();

	void setSlowToStart(const float slowStartPossibility);
	void setConcentrationOutwriter(unique_ptr<ConcentrationOutwriter> writer);

	virtual void calculateResults(int timeSteps);

	virtual const unsigned int getNumberOfCars() const;			//only use for testing
	const int getSpeedAtPosition(unsigned int pos) const;       //only use for testing
	unsigned int getRoadLength() const;
	unsigned int getMaxVelocity() const;

	void dispCurrentRoad() const;
	virtual void dispResults();
	void writeResultsToFile() const;


protected:
	void initDawdle(const float dawdlePossibility);
	void initResults();

	void initResultsForCalculation();

	virtual void calculateTimestep(unsigned int step);
	void putCurrentIntoResults(unsigned int step);
	void writeConcentration();
	void switchCurrentNext();

	vector<unsigned int> findCarIndicesInCurrent() const;
	virtual unsigned int getGapAfterCar(unsigned int carIndex, unsigned int speed, int neighbor);
	unsigned int adjustGapConsideringSafetyDistance(unsigned int speed);

	void applyRules(unsigned int carIndex);
	void accelerateCar(unsigned int &speed);
	void breakCar(unsigned int carIndex, unsigned int &speed);
	void dawdleCar(unsigned int carIndex, unsigned int &speed);
	virtual void moveCar(unsigned int carIndex, unsigned int speed);

	virtual void visualizeSafetyDistance();

protected:
	//temporary varables for calculation
	unsigned int gap;
	double randomNumber;
};

