#pragma once
#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <random> 
#include <memory>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Utilities/RandomHelper.h"
#include "RoadNetwork/RoadNetworkData.h"
#include "Output/ConcentrationOutwriter.h"

using namespace std;

class VF_PUBLIC OneWayRoad //Periodic BC
{
protected:
	unique_ptr<RoadNetworkData> road;

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

	bool saveResults = false;
	vector<vector<int> > results;		//saves the results of the calculation; x-axis = timesteps, y axis = positions

	unique_ptr<ConcentrationOutwriter> concWriter = nullptr;

public:
	OneWayRoad(unique_ptr<RoadNetworkData> road, const float dawdlePossibility);
	OneWayRoad();
	OneWayRoad(const OneWayRoad&) = delete;
	virtual ~OneWayRoad();

	void setSlowToStart(const float slowStartPossibility);
	void setConcentrationOutwriter(unique_ptr<ConcentrationOutwriter> writer);
	void setSaveResults(bool saveResults);

	virtual void initCalculation(int timeSteps);
	virtual void calculateTimestep(unsigned int step);
	void loopTroughTimesteps();

	virtual const unsigned int getNumberOfCars() const;			//only use for testing
	const int getSpeedAtPosition(unsigned int pos) const;       //only use for testing
	unsigned int getRoadLength() const;
	unsigned int getMaxVelocity() const;

	void dispCurrentRoad() const;
	virtual void dispResults();
	void writeResultsToFile() const;

	virtual void visualizeVehicleLengthForVTK();
	virtual const std::vector<int>& getVehiclesForVTK();

protected:
	void initDawdle(const float dawdlePossibility);
	void initResults();

	void initResultsForCalculation();


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
	
	void checkCurrentForSafetyDistance();
	virtual void visualizeSafetyDistanceForConsole();

protected:
	//temporary varables for calculation
	unsigned int gap;
	double randomNumber;
};

