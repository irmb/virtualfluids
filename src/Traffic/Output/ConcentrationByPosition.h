#pragma once
#include <VirtualFluidsDefinitions.h>

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcentrationByPosition:
	public ConcentrationOutwriter
{
public:
	ConcentrationByPosition(uint roadlength, uint maxSpeed = 0, uint maxAcceleration = 0);
	~ConcentrationByPosition() {};

	virtual void calculateConcForSingleCar(uint index, uint speed = 0, uint acceleration = 0);
	virtual void setMaxSpeedAndAcceleration(uint maxSpeed, uint maxAcceleration) {};


};

