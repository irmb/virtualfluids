#ifndef DEFINE_BCS_H
#define DEFINE_BCS_H

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Parameter/Parameter.h"

extern "C" void findQ27(Parameter* para);

extern "C" void findBC27(Parameter* para);

extern "C" void findPressQShip(Parameter* para);

#endif
