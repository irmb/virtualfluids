#ifndef L2_NORM_TEST_STRUCT_H
#define L2_NORM_TEST_STRUCT_H

#include <memory>
#include <vector>

class L2NormTest;
class L2NormInformation;
class L2NormPostProcessingStrategy;


struct L2NormTestStruct
{
	std::shared_ptr<L2NormInformation> logFileInfo;

	std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > postProcessingStrategies;
	std::vector<std::shared_ptr<L2NormTest> > tests;
};

#endif