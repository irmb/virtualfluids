#pragma once
#include <VirtualFluidsDefinitions.h>

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcentrationByPosition:
	public ConcentrationOutwriter
{
public:
	ConcentrationByPosition(uint roadlength, uint maxSpeed = 0);
	ConcentrationByPosition(uint roadlength, float* concArrayStart, uint concArraySize, uint maxSpeed = 0);
	~ConcentrationByPosition() {};

	virtual void calculateConcForSingleCar(uint index, uint oldSpeed = 0, uint speed = 0);
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed = 0, uint speed = 0);

	virtual void setMaxSpeedAndAcceleration(uint maxSpeed, uint maxAcceleration) {};


};

