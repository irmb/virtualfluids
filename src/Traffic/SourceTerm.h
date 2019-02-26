#pragma once
#include <VirtualFluidsDefinitions.h>

#include <iostream>
#include <random>
#include "Source.h"
#include "OneWayRoadSSJ.h"

using namespace std;

class VF_PUBLIC SourceTerm:
	public Source
{
private:
	SourceData data;

	mt19937 engine = RandomHelper::make_engine();
	uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
	uniform_int_distribution<unsigned int> distInt{ 0, 1 };

public:
	SourceTerm(const unsigned int sourceIndex, const float sourcePossibility, unsigned int maxVelocity);
	~SourceTerm();

	virtual unsigned int getIndex() const;
	virtual float getPossibility() const;
	virtual unsigned int getSourceCar();

private:
	unsigned int randomSpeed();

private:
	//variables for temporaray calculations
	float randomNumber;
};

