#pragma once

#include <vector>
#include <iostream>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

class VF_PUBLIC ConcentrationOutwriter
{
public:
	virtual void resetConcentrations();
	virtual void calculateConcForSingleCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;

	void dispConcentrations();

protected:
	void putConcIntoArrayOrVector(uint index, real conc);
	void addConcToArrayOrVector(uint index, real conc);

protected:
	std::vector<real> concentration;
	bool useLBMConcArray = false;
	real* concArrayStart;
	uint roadLength;

private:
	void dispSingleConcentration(real conc);
	
};

