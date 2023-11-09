#ifndef FIND_Q_TEMP_H
#define FIND_Q_TEMP_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"


void findTempPress(Parameter* para);

void findKforTempPress(Parameter* para);

void findTempVel(Parameter* para);

void findKforTempVel(Parameter* para);

void findTemp(Parameter* para);

void findKforTemp(Parameter* para);


#endif
