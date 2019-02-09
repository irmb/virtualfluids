#ifndef TAYLOR_GREEN_VORTEX_UX_PARAMETER_STRUCT_H
#define TAYLOR_GREEN_VORTEX_UX_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities/Structs/BasicSimulationParameterStruct.h"

struct TaylorGreenVortexUxParameterStruct
{
	std::shared_ptr<BasicSimulationParameterStruct> basicSimulationParameter;

	double ux;
	double amplitude;
	int basicTimeStepLength;
	double l0;
	double rho0;
	std::string vtkFilePath;
};
#endif