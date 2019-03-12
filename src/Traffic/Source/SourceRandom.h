#pragma once
#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <random>
#include "Source.h"
#include "Utilities\RandomHelper.h"
#include "Utilities\invalidInput_error.h"


class VF_PUBLIC SourceRandom:
	public Source
{
private:
	SourceData data;

	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
	std::uniform_int_distribution<uint> distInt{ 0, 1 };

public:
	SourceRandom(const uint sourceIndex, const real sourcePossibility, uint maxVelocity);
	~SourceRandom();

	virtual uint getIndex() const;
	virtual real getPossibility() const;
	virtual uint getSourceCar();

private:
	uint randomSpeed();

private:
	//variables for temporaray calculations
	real randomNumber;
};

