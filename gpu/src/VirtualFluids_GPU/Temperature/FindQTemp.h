#ifndef FIND_Q_TEMP_H
#define FIND_Q_TEMP_H

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Parameter/Parameter.h"


extern "C" void findTempPress(Parameter* para);

extern "C" void findKforTempPress(Parameter* para);

extern "C" void findTempVel(Parameter* para);

extern "C" void findKforTempVel(Parameter* para);

extern "C" void findTemp(Parameter* para);

extern "C" void findKforTemp(Parameter* para);


#endif
