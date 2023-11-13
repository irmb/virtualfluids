#ifndef PHI_TEST_LOGFILE_INFORMATION_H
#define PHI_TEST_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"
#include "Utilities/Test/TestStatus.h"

#include <memory>
#include <vector>

class PhiTest;
struct PhiTestParameterStruct;

class PhiTestLogFileInformation : public TestLogFileInformation
{
public:
    static std::shared_ptr<PhiTestLogFileInformation> getNewInstance(std::shared_ptr<PhiTestParameterStruct> testPara);

    std::string getOutput();
    void addTestGroup(std::vector<std::shared_ptr<PhiTest> > tests);

private:
    PhiTestLogFileInformation() {};
    PhiTestLogFileInformation(std::shared_ptr<PhiTestParameterStruct> testPara);

    void fillMyData(std::vector<std::shared_ptr<PhiTest> > testGroup);

    std::vector<std::vector<std::shared_ptr<PhiTest> > > testGroups;
    unsigned int startTimeStepCalculation, endTimeStepCalculation;
    std::vector<int> lx;
    std::vector<int> lxForErase;
    std::vector<double> phiDiff;
    std::vector<double> orderOfAccuracy;
    std::vector<std::string> dataToCalc;
    std::vector<TestStatus> status;
};
#endif