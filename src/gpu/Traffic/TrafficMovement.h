#pragma once

#include <vector>
#include <random> 
#include <memory>

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include "RoadNetwork/RoadNetworkData.h"
#include "Utilities/RandomHelper.h"

#include "Output/ConcentrationOutwriter.h"
#include "Output/CarDisplay.h"

#include "Traffic_export.h"

class TrafficLogger;
//class ConcentrationOutwriter;
//class CarDisplay;
class TrafficTimestep;

class TRAFFIC_EXPORT TrafficMovement
{
public:
	TrafficMovement(std::shared_ptr<RoadNetworkData> road, const real dawdlePossibility);
	//TrafficMovement() {};
	TrafficMovement(const TrafficMovement&) = delete;
	~TrafficMovement();

	//setUp
	void setSlowToStart(const real slowStartPossibility);
	void setMaxAcceleration(uint maxAcceleration);
	void setConcentrationOutwriter(uint roadlength, real* concArrayStart = 0);
	void setSaveResultsTrue(uint timeSteps);
	void setUseGPU(real * pConcArray = nullptr, int* naschVelocity = nullptr);
	void setUseLogger();

	//timpestep
	void loopTroughTimesteps(uint numberOfTimesteps);
	void calculateTimestep(uint step);

	//get
	uint getRoadLength() const;
	uint getMaxVelocity() const;
	uint getNumberOfCars() const;			//only use for testing
	int getSpeedAtPosition(uint pos) const;       //only use for testing
	real getDawdlePossibility();
	bool getUseSlowToStart();
	real getSlowToStartPossibility();
	uint getMaxAcceleration();


	//methods used by junctions and sources
	uint getGapAfterOutCell(uint outCellIndex, uint speed);
	void moveJunctionCar(uint outCellIndex, uint remainingDistance, uint speed, uint oldSpeed);
	void writeConcentrationForJunction(uint inCellIndex, uint oldSpeed, uint speed);

	//vtk
	void visualizeVehicleLengthForVTK();
	const std::vector<int> & getVehiclesForVTK();

	//for debugging
	void checkSpeed(uint speed);

private:
	//init
	void initDawdle(const real dawdlePossibility);
	void checkCurrentForSafetyDistance();

	//calculate timestep
	void calculateSourceStep();
	void calculateJunctionStep();
	void switchCurrentNext();

	//gap
	uint getGapAfterCar(uint carIndex, uint speed, int neighbor);
	uint getGapToSink(int neighbor, uint i, uint speed);
	uint getGapToJunction(int neighbor, uint i, uint speed, uint currentCell);
	uint adjustGapToSafetyDist(uint gap);

	//getVectorIndex
	std::shared_ptr<Junction>& getJunctionFromNeighbor(int neighbor);
	std::shared_ptr<Sink>& getSinkFromNeighbor(int neighbor);

	//apply rules
	void applyRules(uint carIndex);
	void accelerateCar(uint &speed);
	void brakeCar(uint carIndex, uint &speed);
	void dawdleCar(uint carIndex, uint &speed);
	void moveCar(const uint carIndex, uint speed);
	uint iterateNeighborsInMove(uint &currentCell, uint speed, int &neighbor);

	//disp
	void dispResults();
	void dispCurrentConcFromGPU();

	//pollution
	void writeConcentration(uint index, uint oldSpeed);

	//gpu
	void copyDevToHost();

private:
	std::shared_ptr<RoadNetworkData> road;
	std::unique_ptr<ConcentrationOutwriter> concWriter = nullptr;
	std::unique_ptr<CarDisplay> display = nullptr;

	bool useGPU = false;
	std::unique_ptr<TrafficTimestep> gpuCalculation;
	bool copiedDevToHost = false;

	bool useLogger = false;

	uint currentStep = 0;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<real> distFloat{ 0.0, 1.0 };

private:
	//temporary variables for calculation
	uint gap;
	float randomNumber;
};

