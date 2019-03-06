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
	virtual void resetConcVector();
	virtual void calculateConcForSingleCar(uint index, uint speed = 0, uint acceleration = 0) = 0;

	void dispConcentration();
	const std::vector<float> & getConcentrations();

protected:
	std::vector<float> concentration;
};

