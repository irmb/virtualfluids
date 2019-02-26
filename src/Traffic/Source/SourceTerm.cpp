#pragma once

#include "SourceTerm.h"


SourceTerm::SourceTerm(const unsigned int sourceIndex, const float sourcePossibility, unsigned int maxVelocity) 
{
	data.sourceIndex = sourceIndex;
	data.maxVelocity = maxVelocity;

	try {
		if (sourcePossibility >= 0 && sourcePossibility <= 1) {
			data.sourcePossibility = sourcePossibility;
			uniform_int_distribution<unsigned int> distInt2{ 0, maxVelocity };
			distInt = distInt2;
		}
		else {
			throw invalidInput_error("possibility of a car leaving the sink should be between 0 and 1");
		}
	}
	catch (const exception& e) {
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	};
}



SourceTerm::~SourceTerm()
{
}

unsigned int SourceTerm::getIndex() const
{
	return data.sourceIndex;
}

float SourceTerm::getPossibility() const
{
	return data.sourcePossibility;
}


unsigned int SourceTerm::getSourceCar()
{
	randomNumber = distFloat(engine);
	if (randomNumber < data.sourcePossibility) {
		return randomSpeed();
	}
	return -1;
}


unsigned int SourceTerm::randomSpeed()
{
	return distInt(engine);
}
