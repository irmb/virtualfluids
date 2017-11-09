#ifndef INIT_LATTICE_H
#define INIT_LATTICE_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Temperature/FindTemperature.h"

extern "C" void initLattice(Parameter* para);

#endif
