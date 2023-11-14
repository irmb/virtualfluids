#ifndef SHEAR_WAVE_PARAMETER_STRUCT_H
#define SHEAR_WAVE_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities/Structs/BasicSimulationParameterStruct.h"

struct ShearWaveParameterStruct
{
    std::shared_ptr<BasicSimulationParameterStruct> basicSimulationParameter;

    double ux;
    double uz;
    int basicTimeStepLength;
    double l0;
    double rho0;
    std::string vtkFilePath;
    std::vector<std::string> dataToCalcTests;
};
#endif