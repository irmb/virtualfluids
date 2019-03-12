#pragma once

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcBySpeedAndAcceleration :
	public ConcentrationOutwriter
{
public:
	ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed);
	ConcBySpeedAndAcceleration(uint roadlength, real* concArrayStart, uint maxSpeed=0);	
	~ConcBySpeedAndAcceleration() {};

	virtual void calculateConcForSingleCar(uint index, uint oldSpeed, uint speed);
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed, uint speed);


private:
	real chooseConc(uint oldSpeed, uint speed);

private:
	real maxSpeed = 0;
};

