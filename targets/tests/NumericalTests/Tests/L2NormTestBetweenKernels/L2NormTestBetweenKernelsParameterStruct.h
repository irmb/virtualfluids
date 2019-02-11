#ifndef L2_NORM_TEST_BETWEEN_KERNELS_PARAMETER_STRUCT_H
#define L2_NORM_TEST_BETWEEN_KERNELS_PARAMETER_STRUCT_H

#include <memory>
#include <string>
#include <vector>

#include "Utilities/Structs/BasicTestParameterStruct.h"

struct L2NormTestBetweenKernelsParameterStruct
{
	std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

	std::string basicKernel;

	std::vector<std::string> kernelsToTest;
	std::vector<int> timeSteps;

	std::string normalizeWith;
};

#endif 