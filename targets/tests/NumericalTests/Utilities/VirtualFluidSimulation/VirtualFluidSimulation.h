#ifndef VIRTUAL_FLUID_SIMULATION_H
#define VIRTUAL_FLUID_SIMULATION_H

#include "VirtualFluids_GPU/LBM/LB.h"

#include <string>
#include <memory>

class Parameter;
class GridProvider;
class DataWriter;
class Calculator;
class TestResults;

class VirtualFluidSimulation
{
public:
	virtual std::shared_ptr<Parameter> getParameter() = 0;
	virtual std::shared_ptr<GridProvider> getGrid() = 0;
	virtual std::shared_ptr<DataWriter> getDataWriter() = 0;
	virtual std::shared_ptr<Calculator> getCalculator() = 0;
private:

};
#endif