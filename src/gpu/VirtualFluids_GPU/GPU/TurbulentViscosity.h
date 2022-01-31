#ifndef TURBULENT_VISCOSITY_H_
#define TURBULENT_VISCOSITY_H_

#include "Core/DataTypes.h"

class Parameter;

extern "C" void calcTurbulentViscosityAMD(Parameter* para, int level);

#endif //TURBULENT_VISCOSITY_H_