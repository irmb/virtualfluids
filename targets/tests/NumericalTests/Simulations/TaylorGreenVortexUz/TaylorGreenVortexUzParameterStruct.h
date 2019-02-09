#ifndef TAYLOR_GREEN_VORTEX_UZ_PARAMETER_STRUCT_H
#define TAYLOR_GREEN_VORTEX_UZ_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities/Structs/BasicSimulationParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

struct TaylorGreenVortexUzParameterStruct
{
	std::shared_ptr<BasicSimulationParameterStruct> basicSimulationParameter;

	double uz;
	double amplitude;
	int basicTimeStepLength;
	double l0;
	double rho0;
	std::string vtkFilePath;
};
#endif