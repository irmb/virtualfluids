#ifndef NY_TEST_LOGFILE_INFORMATION_H
#define NY_TEST_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"
#include "Utilities/Test/TestStatus.h"

#include <memory>
#include <vector>

class NyTest;
struct NyTestParameterStruct;

class NyTestLogFileInformation : public TestLogFileInformation
{
public:
    static std::shared_ptr<NyTestLogFileInformation> getNewInstance(std::shared_ptr<NyTestParameterStruct> testPara);

    std::string getOutput();
    void addTestGroup(std::vector<std::shared_ptr<NyTest> > tests);

private:
    NyTestLogFileInformation() {};
    NyTestLogFileInformation(std::shared_ptr<NyTestParameterStruct> testPara);

    void fillMyData(std::vector<std::shared_ptr<NyTest> > testGroup);

    std::vector<std::vector<std::shared_ptr<NyTest> > > testGroups;
    unsigned int startTimeStepCalculation, endTimeStepCalculation;
    std::vector<int> lx;
    std::vector<int> lxForErase;
    std::vector<double> ny, nyDiff;
    std::vector<double> orderOfAccuracyNyDiff;
    std::vector<std::string> dataToCalc;
    std::vector<TestStatus> status;
};
#endif