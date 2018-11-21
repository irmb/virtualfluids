#ifndef VIRTUAL_FLUID_SIMULATION_H
#define VIRTUAL_FLUID_SIMULATION_H

#include "VirtualFluids_GPU/LBM/LB.h"

class VirtualFluidSimulation
{
public:
	virtual void run() = 0;
private:

};
#endif