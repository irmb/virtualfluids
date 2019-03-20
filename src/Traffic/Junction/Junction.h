#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>

#include "JunctionData.h"

class TrafficMovement;

class VF_PUBLIC Junction
{
public:
	virtual void checkOutCellIndices(const uint roadLength) const = 0;

	virtual void setCellIndexForNoUTurn(std::vector<int> carCanNotEnterThisOutCell) = 0;

	virtual bool acceptsCar(uint cellIndex) = 0; //determines if a car can enter the junction
	virtual void registerCar(uint cellIndex, uint numberOfCellsAlreadyMoved, uint speed, uint oldSpeed) = 0; //registers all cars entering the junction
	virtual void calculateTimeStep(TrafficMovement &road) = 0;
	virtual void updateJunction() = 0;

	virtual const std::vector<uint>& getInCellIndices() const  = 0;

	virtual void dispJunction(const uint index, const uint roadLength) const = 0;
	virtual const uint getNumCarsOnJunction() const = 0;
};