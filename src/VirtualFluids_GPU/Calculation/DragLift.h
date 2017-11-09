#ifndef DragLift_H
#define DragLift_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Utilities/StringUtil.hpp"
#include "Parameter/Parameter.h"

extern "C" void calcDragLift(Parameter* para, int lev);
extern "C" void allocDragLift(Parameter* para);
extern "C" void printDragLift(Parameter* para, int timestep);

#endif
