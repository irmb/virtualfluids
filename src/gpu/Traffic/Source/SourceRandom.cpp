#include "SourceRandom.h"

#include <iostream>

#include "Utilities/invalidInput_error.h"

SourceRandom::SourceRandom(const uint sourceIndex, const real sourcePossibility, uint maxVelocity) 
{
	data.sourceIndex = sourceIndex;
	data.maxVelocity = maxVelocity;

	try {
		if (sourcePossibility >= 0 && sourcePossibility <= 1) {
			data.sourcePossibility = sourcePossibility;
			std::uniform_int_distribution<uint> distInt2{ 0, maxVelocity };
			distInt = distInt2;
		}
		else {
			throw invalidInput_error("possibility of a car leaving the sink should be between 0 and 1");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	};
}



SourceRandom::~SourceRandom()
{
}

uint SourceRandom::getIndex() const
{
	return data.sourceIndex;
}

real SourceRandom::getPossibility() const
{
	return data.sourcePossibility;
}


uint SourceRandom::getSourceCar()
{
	randomNumber = distFloat(engine);
	if (randomNumber < data.sourcePossibility) {
		return randomSpeed();
	}
	return -1;
}


uint SourceRandom::randomSpeed()
{
	return distInt(engine);
}
