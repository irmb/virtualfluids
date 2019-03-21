#pragma once

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcentrationByPosition:
	public ConcentrationOutwriter
{
public:
	ConcentrationByPosition(uint roadlength, real* concArrayStart = nullptr, uint maxSpeed = 0);
	~ConcentrationByPosition() {};

	virtual void calculateConcForSingleCar(uint index, uint oldSpeed = 0, uint speed = 0);
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed = 0, uint speed = 0);
};

