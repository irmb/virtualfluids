#ifndef Particles_H
#define Particles_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Utilities/StringUtil.hpp"
#include "Parameter/Parameter.h"

//extern "C" void calcDragLift(Parameter* para, int lev);
extern "C" void allocParticles(Parameter* para);
extern "C" void initParticles(Parameter* para);
extern "C" void propagateParticles(Parameter* para, unsigned int t);
extern "C" void copyAndPrintParticles(Parameter* para, unsigned int t, bool isInit);

extern "C" void rearrangeGeometry(Parameter* para);

#endif
