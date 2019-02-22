#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <iostream>

#include <iomanip>	//formatting output streams

class VF_PUBLIC ConcentrationOutwriter
{
public:
	virtual void writeToArray(const std::vector<int> & currentCarDistribution) = 0;
};

