#ifndef CalcMedian_H
#define CalcMedian_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

extern "C" void allocMedian(Parameter* para);
extern "C" void calcMedian(Parameter* para, unsigned int tdiff);
extern "C" void resetMedian(Parameter* para);

#endif
