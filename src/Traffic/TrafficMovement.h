#pragma once
#include <VirtualFluidsDefinitions.h>

#include <iostream>

#include <vector>
#include <random> 
#include <memory>
#include <stdexcept>

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Utilities/RandomHelper.h"
#include "RoadNetwork/RoadNetworkData.h"
#include "Output/ConcentrationOutwriter.h"
#include "Output/CarDisplay.h"


class VF_PUBLIC TrafficMovement //Periodic BC
{
public:
	TrafficMovement(std::unique_ptr<RoadNetworkData> road, const float dawdlePossibility);
	//TrafficMovement() {};
	TrafficMovement(const TrafficMovement&) = delete;
	~TrafficMovement();

	//setUp
	void setSlowToStart(const float slowStartPossibility);
	void setConcentrationOutwriter(std::unique_ptr<ConcentrationOutwriter> writer);
	void setSaveResultsTrue();
	void initCalculation(uint timeSteps);

	//timpestep
	void loopTroughTimesteps();
	void calculateTimestep(uint step);

	//get
	const uint getRoadLength() const;
	const uint getMaxVelocity() const;
	const uint getNumberOfCars() const;			//only use for testing
	const int getSpeedAtPosition(uint pos) const;       //only use for testing

	//methods used by junctions and sinks
	uint getGapAfterOutCell(uint outCellIndex, uint speed);
	void moveJunctionCar(uint outCellIndex, uint remainingDistance, uint speed);

	//disp
	void dispResults();

	//vtk
	void visualizeVehicleLengthForVTK();
	const std::vector<int>& getVehiclesForVTK();

	//pollution
	const std::vector<float>& getConcentrations();


private:
	//init
	void initDawdle(const float dawdlePossibility);
	void checkCurrentForSafetyDistance();

	//calculate timestep
	std::vector<uint> findCarIndicesInCurrent() const;
	void calculateSourceStep();
	void switchCurrentNext();

	//gap
	uint getGapAfterCar(uint carIndex, uint speed, int neighbor);
	uint getGapToSink(int neighbor, uint i, uint speed);
	uint getGapToJunction(int neighbor, uint i, uint speed, uint currentCell);
	uint adjustGapToSafetyDist(uint speed);

	//getVectorIndex
	std::shared_ptr<Junction>& getJunctionFromNeighbor(int neighbor);
	std::shared_ptr<Sink>& getSinkFromNeighbor(int neighbor);

	//apply rules
	void applyRules(uint carIndex);
	void accelerateCar(uint &speed);
	void breakCar(uint carIndex, uint &speed);
	void dawdleCar(uint carIndex, uint &speed);
	void moveCar(uint carIndex, uint speed);
	uint iterateNeighborsInMove(uint &currentCell, uint speed, int &neighbor);

	void writeConcentration(uint index);


private:
	std::unique_ptr<RoadNetworkData> road;
	std::unique_ptr<ConcentrationOutwriter> concWriter = nullptr;
	std::unique_ptr<CarDisplay> display = nullptr;

	std::vector<int> *pcurrent;
	std::vector<int> *pnext;
	std::vector<int> *pdummy; 

	float dawdlePossibility;

	bool useSlowToStart = false;
	float slowStartPossibility;

	uint timeSteps;
	uint currentStep;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };

private:
	//temporary varables for calculation
	uint gap;
	double randomNumber;
};

