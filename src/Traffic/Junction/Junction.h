#pragma once
#include <vector>

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include "JunctionData.h"

class TrafficMovement;

class VF_PUBLIC Junction
{
public:
	virtual void checkOutCellIndices(const uint roadLength) const = 0;

	virtual void setCellIndicesForNoUTurn(std::vector<int> carCanNotEnterThisOutCell) = 0;

	virtual bool acceptsCar(uint cellIndex) = 0; //determines if a car can enter the junction
	virtual void registerCar(uint cellIndex, uint numberOfCellsAlreadyMoved, uint speed, uint oldSpeed) = 0; //registers all cars entering the junction
	virtual void calculateTimeStep(TrafficMovement &road, uint currentTimestep) = 0;

	virtual const std::vector<uint>& getInCellIndices()const = 0;
	virtual const std::vector<uint>& getOutCellIndices() const = 0;
	virtual const std::vector<bool>& getCarCanEnter() const = 0;
	virtual const std::vector<int>& getCarsOnJunction()const = 0;
	virtual const std::vector<uint>& getAlreadyMoved()const = 0;
	virtual const std::vector<uint>& getOldSpeeds()const = 0;
	virtual const std::vector<int>& getCarCanNotEnterThisOutCell() const = 0;
	virtual uint getTrafficLightSwitchTime()const = 0;

	virtual void dispJunction(const uint index, const uint roadLength) const = 0;
	virtual uint getNumCarsOnJunction() const = 0;
};