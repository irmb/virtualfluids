#ifndef PHI_AND_NY_TEST_STRUCT_H
#define PHI_AND_NY_TEST_STRUCT_H

#include <memory>
#include <vector>

class PhiAndNyInformation;
class PhiAndNyTestPostProcessingStrategy;
class PhiAndNyTest;

struct PhiAndNyTestStruct
{
	std::shared_ptr<PhiAndNyInformation> logFileInfo;
	std::vector<std::shared_ptr<PhiAndNyTest> > tests;
};
#endif 