#pragma once
#include <VirtualFluidsDefinitions.h>

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcBySpeedAndAcceleration :
	public ConcentrationOutwriter
{
public:
	ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed, uint maxAcceleration);
	~ConcBySpeedAndAcceleration() {};

	virtual void calculateConcForSingleCar(uint index, uint speed, uint acceleration);

private:
	float maxSpeed = 0;
	float maxAcceleration;
};

