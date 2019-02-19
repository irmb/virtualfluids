#ifndef PHI_TEST_STRUCT_H
#define PHI_TEST_STRUCT_H

#include <memory>
#include <vector>

class PhiTestLogFileInformation;
class PhiTest;

struct PhiTestStruct
{
	std::shared_ptr<PhiTestLogFileInformation> logFileInfo;
	std::vector<std::shared_ptr<PhiTest> > tests;
};
#endif 