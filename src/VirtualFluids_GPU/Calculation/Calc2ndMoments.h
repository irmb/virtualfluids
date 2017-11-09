#ifndef Calc2ndMoments_H
#define Calc2ndMoments_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

//2nd
extern "C" void alloc2ndMoments(Parameter* para);
extern "C" void init2ndMoments(Parameter* para);
extern "C" void calc2ndMoments(Parameter* para);

//3rd
extern "C" void alloc3rdMoments(Parameter* para);
extern "C" void init3rdMoments(Parameter* para);
extern "C" void calc3rdMoments(Parameter* para);

//higher order
extern "C" void allocHigherOrderMoments(Parameter* para);
extern "C" void initHigherOrderMoments(Parameter* para);
extern "C" void calcHigherOrderMoments(Parameter* para);

#endif
