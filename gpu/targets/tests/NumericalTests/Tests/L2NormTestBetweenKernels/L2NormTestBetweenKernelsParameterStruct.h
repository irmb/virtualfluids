#ifndef L2_NORM_TEST_BETWEEN_KERNELS_PARAMETER_STRUCT_H
#define L2_NORM_TEST_BETWEEN_KERNELS_PARAMETER_STRUCT_H

#include <memory>
#include <string>
#include <vector>

#include "Utilities/Structs/BasicTestParameterStruct.h"

#include "VirtualFluids_GPU/Kernel//Utilities/KernelType.h"

struct L2NormTestBetweenKernelsParameterStruct
{
	std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

	KernelType basicKernel;

	std::vector<KernelType> kernelsToTest;
	std::vector<int> timeSteps;

	std::vector<std::string> normalizeData;
};

#endif 