#ifndef BASIC_SIMULATION_PARAMETER_STRUCT_H
#define BASIC_SIMULATION_PARAMETER_STRUCT_H

#include <vector>

struct BasicSimulationParameterStruct
{
	std::vector<unsigned int> devices;
	unsigned int numberOfTimeSteps;
};

#endif