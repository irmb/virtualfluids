#pragma once

#include "ConcentrationOutwriter.h"

class VIRTUALFLUIDS_GPU_EXPORT ConcBySpeedAndAcceleration :
	public ConcentrationOutwriter
{
public:
	ConcBySpeedAndAcceleration(uint roadlength, real* concArrayStart = 0);
	~ConcBySpeedAndAcceleration() {};

	virtual void calculateConcForSingleCar(uint index, uint oldSpeed, uint speed);
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed, uint speed);
	virtual void calculateConcForAllCars(const std::vector<int> oldSpeeds, const std::vector<int> newSpeeds);


private:
	real chooseConc(uint oldSpeed, uint speed);

};

