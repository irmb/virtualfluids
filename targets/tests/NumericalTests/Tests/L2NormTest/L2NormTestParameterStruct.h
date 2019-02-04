#ifndef L2_NORM_TEST_PARAMETER_STRUCT_H
#define L2_NORM_TEST_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities\Structs\BasicTestParameterStruct.h"

struct L2NormTestParameterStruct
{
	std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

	double maxDiff;
	unsigned int basicTimeStep;
	unsigned int divergentTimeStep;
};

#endif