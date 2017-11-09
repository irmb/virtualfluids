#ifndef FIND_TEMPERATURE_H
#define FIND_TEMPERATURE_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Temperature/FindQTemp.h"

extern "C" void initTemperatur(Parameter* para, int lev);

extern "C" void findTempSim(Parameter* para);
							
extern "C" void findTempVelSim(Parameter* para);
								
extern "C" void findTempPressSim(Parameter* para);

#endif
