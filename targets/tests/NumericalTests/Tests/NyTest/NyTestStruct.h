#ifndef NY_TEST_STRUCT_H
#define NY_TEST_STRUCT_H

#include <memory>
#include <vector>

class NyTestLogFileInformation;
class NyTest;

struct NyTestStruct
{
	std::shared_ptr<NyTestLogFileInformation> logFileInfo;
	std::vector<std::shared_ptr<NyTest> > tests;
};
#endif 