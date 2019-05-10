//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#ifndef CalcMedian_H
#define CalcMedian_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

extern "C" void allocMedian(Parameter* para);
extern "C" void calcMedian(Parameter* para, uint tdiff);
extern "C" void resetMedian(Parameter* para);

//Advection-Diffusion
extern "C" void allocMedianAD(Parameter* para);
extern "C" void calcMedianAD(Parameter* para, uint tdiff);
extern "C" void resetMedianAD(Parameter* para);

#endif
