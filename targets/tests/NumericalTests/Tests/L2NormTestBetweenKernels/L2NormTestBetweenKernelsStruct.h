#ifndef L2_NORM_TEST_BETWEEN_KERNELS_STRUCT_H
#define L2_NORM_TEST_BETWEEN_KERNELS_STRUCT_H

#include <memory>
#include <string>
#include <vector>

class L2NormTestBetweenKernels;
class L2NormBetweenKernelPostProcessingStrategy;
class L2NormBetweenKernelsInformation;

struct L2NormTestBetweenKernelsStruct
{
	std::shared_ptr<L2NormBetweenKernelsInformation> logFileInfo;

	std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests;
};

#endif 