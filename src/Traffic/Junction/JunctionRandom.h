#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <iostream>
#include <random> 
#include <algorithm>

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Junction.h"
#include "TrafficMovement.h"


class VF_PUBLIC JunctionRandom :
	public Junction
{

private:
	JunctionData data;

public:
	JunctionRandom(const std::vector<uint> &inCellIndices, const std::vector<uint> &outCellIndices);
	~JunctionRandom() {};

	virtual void setCellIndexForNoUTurn(std::vector<int> carCanNotEnterThisOutCell);

	virtual bool acceptsCar(uint cellIndex); //determines if a car can enter the junction
	virtual void registerCar(uint cellIndex, uint numberOfCellsAlreadyMoved,  uint speed); //registers all cars entering the junction
	virtual void calculateTimeStep(TrafficMovement& road);
	virtual void updateJunction();

	virtual const std::vector<uint>& getInCellIndices() const;
	
	virtual void dispJunction(const uint index, uint roadLength) const;
	virtual const uint getNumCarsOnJunction() const; 

	virtual void checkOutCellIndices(uint roadLength); 

private:
	uint getInCellsVectorIndex(uint cellIndex);

	void applyRules(int &carSpeed,const int &index, TrafficMovement& road);
	void breakCar(uint outCellIndex, int &speed, uint &remainingDistance, TrafficMovement& road);
	void moveCar(uint outCell, int & carSpeed, const int & index, TrafficMovement& road);
	int chooseOutCell(const int & index);


private:
	//variables for temporaray calculations
	uint remainingDistance;
	int outCell;
	uint random;
	uint gap;
	int index;
};

