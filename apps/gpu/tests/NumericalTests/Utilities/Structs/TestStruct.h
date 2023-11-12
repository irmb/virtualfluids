#ifndef TEST_STRUCT_H
#define TEST_STRUCT_H

#include <memory>
#include <string>
#include <vector>

class TestLogFileInformation;
class Test;

struct TestStruct
{
    std::shared_ptr<TestLogFileInformation> logFileInfo;
    std::vector<std::shared_ptr<Test> > tests;
    std::string testName;
};
#endif 