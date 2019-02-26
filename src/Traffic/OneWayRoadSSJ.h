#pragma once
#include <VirtualFluidsDefinitions.h>

#include "OneWayRoad.h"

class VF_PUBLIC OneWayRoadSSJ :
	public OneWayRoad
{
	
public:
	OneWayRoadSSJ(unique_ptr<RoadNetworkData> road, const float dawdlePossibility);
	OneWayRoadSSJ();
	virtual ~OneWayRoadSSJ();

	virtual void calculateTimestep(unsigned int step);

	unsigned int getGapAfterOutCell(unsigned int outCellIndex, unsigned int speed);
	void moveJunctionCar(unsigned int outCellIndex, unsigned int remainingDistance, unsigned int speed);

	virtual void dispResults();

	virtual const unsigned int getNumberOfCars() const; //only use for testing

private:

	void calculateSourceStep();

	virtual unsigned int getGapAfterCar(unsigned int carIndex, unsigned int speed, int neighbor);
	unsigned int getGapToSink(int neighbor, unsigned int i, unsigned int speed);
	unsigned int getGapToJunction(int neighbor, unsigned int i, unsigned int speed, unsigned int currentCell);

	virtual void moveCar(unsigned int carIndex, unsigned int speed);
	unsigned int iterateNeighborsInMove(unsigned int &currentCell, unsigned int speed, int &neighbor);
	shared_ptr<Junction>& getJunctionFromNeighbor(int neighbor);
	shared_ptr<Sink>& getSinkFromNeighbor(int neighbor);

	virtual void visualizeSafetyDistanceForConsole();
	void dispJunctionsAtCell(unsigned int index) const;
	void dispSinksAtCell(unsigned int index) const;
	void dispSourcesAtCell(unsigned int index) const;
};
