//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#ifndef Outflow_H
#define Outflow_H

#include "LBM/LB.h"

struct LBMSimulationParameter;

void OutflowNonReflecting(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition);

void OutflowNonReflectingPressureCorrection(LBMSimulationParameter* parameterDevice, QforBoundaryConditions* boundaryCondition);

#endif
