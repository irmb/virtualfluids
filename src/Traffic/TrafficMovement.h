#pragma once

#include <vector>
#include <random> 
#include <memory>

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include "RoadNetwork/RoadNetworkData.h"

class ConcentrationOutwriter;
class CarDisplay;

class VF_PUBLIC TrafficMovement
{
public:
	TrafficMovement(std::unique_ptr<RoadNetworkData> road, const real dawdlePossibility);
	//TrafficMovement() {};
	TrafficMovement(const TrafficMovement&) = delete;
	~TrafficMovement();

	//setUp
	void setSlowToStart(const real slowStartPossibility);
	void setConcentrationOutwriter(std::unique_ptr<ConcentrationOutwriter> writer);
	void setSaveResultsTrue(uint timeSteps);

	//timpestep
	void loopTroughTimesteps(uint numberOfTimesteps);
	void calculateTimestep(uint step);

	//get
	const uint getRoadLength() const;
	const uint getMaxVelocity() const;
	const uint getNumberOfCars() const;			//only use for testing
	const int getSpeedAtPosition(uint pos) const;       //only use for testing

	//methods used by junctions and sources
	uint getGapAfterOutCell(uint outCellIndex, uint speed);
	void moveJunctionCar(uint outCellIndex, uint remainingDistance, uint speed, uint oldSpeed);
	void writeConcentrationForJunction(uint inCellIndex, uint oldSpeed, uint speed);

	//vtk
	void visualizeVehicleLengthForVTK();
	const std::vector<int> & getVehiclesForVTK();


private:
	//init
	void initDawdle(const real dawdlePossibility);
	void checkCurrentForSafetyDistance();

	//calculate timestep
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
	void moveCar(const uint carIndex, uint speed);
	uint iterateNeighborsInMove(uint &currentCell, uint speed, int &neighbor);

	//disp
	void dispResults();

	//pollution
	void writeConcentration(uint index, uint oldSpeed);


private:
	std::unique_ptr<RoadNetworkData> road;
	std::unique_ptr<ConcentrationOutwriter> concWriter = nullptr;
	std::unique_ptr<CarDisplay> display = nullptr;

	std::vector<int> *pcurrent;
	std::vector<int> *pnext;
	std::vector<int> *pdummy;

	real dawdlePossibility;

	bool useSlowToStart = false;
	real slowStartPossibility;

	uint currentStep;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<real> distFloat{ 0.0, 1.0 };

private:
	//temporary variables for calculation
	uint gap;
	double randomNumber;
};

