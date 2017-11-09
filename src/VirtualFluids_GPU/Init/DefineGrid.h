#ifndef DEFINE_GRID_H
#define DEFINE_GRID_H

#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"

extern "C" void defineGrid(Parameter* para, Communicator* comm);

#endif
