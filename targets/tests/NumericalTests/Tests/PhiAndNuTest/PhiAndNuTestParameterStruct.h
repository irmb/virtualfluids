#ifndef PHI_AND_NU_TEST_PARAMETER_STRUCT_H
#define PHI_AND_NU_TEST_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities\Structs\BasicTestParameterStruct.h"

struct PhiAndNuTestParameterStruct
{
	std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

	double minOrderOfAccuracy;
	unsigned int startTimeStepCalculation;
	unsigned int endTimeStepCalculation;
};
#endif 