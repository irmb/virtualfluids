#ifndef BASIC_TEST_LOGFILE_INFORMATION_H
#define BASIC_TEST_LOGFILE_INFORMATION_H

#include "../LogFileInformationImp.h" 

#include <iostream>
#include <memory>
#include <vector>

class BasicTestLogFileInformation : public LogFileInformationImp
{
public:
    static std::shared_ptr<BasicTestLogFileInformation> getNewInstance();

    std::string getOutput();
    void addTest(std::string testName, bool testRun);

private:
    BasicTestLogFileInformation();

    void buildOutput();

    bool outputBuild;
    std::vector<std::string> testName;
    std::vector<bool> testRun;
};
#endif