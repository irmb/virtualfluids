#pragma once

#include <vector>
#include <iostream>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include "DrivingStates.h"

class VF_PUBLIC ConcentrationOutwriter
{
public:
	virtual void resetConcentrations();
	virtual void calculateConcForSingleCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;

	void dispConcentrations();

protected:
	void putConcIntoArrayOrVector(uint index, float conc);
	void addConcToArrayOrVector(uint index, float conc);

protected:
	std::vector<float> concentration;
	bool useLBMConcArray = false;
	float* concArrayStart;
	uint concArraySize;


private:
	void dispSingleConcentration(float conc);
	
};

