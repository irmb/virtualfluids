#pragma once
#include <VirtualFluidsDefinitions.h>

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcBySpeedAndAcceleration :
	public ConcentrationOutwriter
{
public:
	ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed);
	ConcBySpeedAndAcceleration(uint roadlength, float* concArrayStart, uint concArraySize, uint maxSpeed);	
	~ConcBySpeedAndAcceleration() {};

	virtual void calculateConcForSingleCar(uint index, uint oldSpeed, uint speed);
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed, uint speed);


private:
	float chooseConc(uint oldSpeed, uint speed);

private:
	float maxSpeed = 0;
	float maxAcceleration;
};

