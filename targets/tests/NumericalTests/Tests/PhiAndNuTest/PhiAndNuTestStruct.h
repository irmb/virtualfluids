#ifndef PHI_AND_NU_TEST_STRUCT_H
#define PHI_AND_NU_TEST_STRUCT_H

#include <memory>
#include <vector>

class PhiAndNuInformation;
class PhiAndNuTestPostProcessingStrategy;
class PhiAndNuTest;

struct PhiAndNuTestStruct
{
	std::shared_ptr<PhiAndNuInformation> logFileInfo;
	std::vector<std::shared_ptr<PhiAndNuTest> > tests;
};
#endif 