#ifndef PHIANDNY_LOGFILE_INFORMATION_H
#define PHIANDNY_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"

#include <memory>
#include <vector>

class PhiAndNyTest;
struct PhiAndNyTestParameterStruct;

class PhiAndNyInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr<PhiAndNyInformation> getNewInstance(std::shared_ptr<PhiAndNyTestParameterStruct> testPara);

	std::string getOutput();
	void addTestGroup(std::vector<std::shared_ptr<PhiAndNyTest> > tests);

private:
	PhiAndNyInformation() {};
	PhiAndNyInformation(std::shared_ptr<PhiAndNyTestParameterStruct> testPara);

	void fillMyData(std::vector<std::shared_ptr<PhiAndNyTest> > testGroup);

	std::vector<std::vector<std::shared_ptr<PhiAndNyTest> > > testGroups;
	unsigned int startTimeStepCalculation, endTimeStepCalculation;
	std::vector<int> lx;
	std::vector<int> lxForErase;
	std::vector<double> phiDiff;
	std::vector<double> ny, nyDiff;
	std::vector<double> orderOfAccuracyPhiDiff;
	std::vector<double> orderOfAccuracyNyDiff;
	std::vector<std::string> dataToCalc;
};
#endif