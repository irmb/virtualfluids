#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>

#include "JunctionData.h"

class OneWayRoadSSJ;

class VF_PUBLIC Junction
{
public:
	virtual ~Junction() {};

	virtual void checkOutCellIndices(unsigned int roadLength) = 0;

	virtual bool acceptsCar(unsigned int cellIndex) = 0; //determines if a car can enter the junction
	virtual void registerCar(unsigned int cellIndex, unsigned int numberOfCellsAlreadyMoved, unsigned int speed) = 0; //registers all cars entering the junction
	virtual void calculateTimeStep(OneWayRoadSSJ& road) = 0;
	virtual void updateJunction() = 0;

	virtual const std::vector<unsigned int>& getInCellIndices() const  = 0;

	virtual void dispJunction(const unsigned int index, unsigned int roadLength) const = 0;
	virtual const unsigned int getNumCarsOnJunction() const = 0;
};