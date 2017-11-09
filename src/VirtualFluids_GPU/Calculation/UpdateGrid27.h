#ifndef UPDATEGRID27_H
#define UPDATEGRID27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"

extern "C" void updateGrid27(Parameter* para, Communicator* comm, int level, int max_level, unsigned int t);

#endif
