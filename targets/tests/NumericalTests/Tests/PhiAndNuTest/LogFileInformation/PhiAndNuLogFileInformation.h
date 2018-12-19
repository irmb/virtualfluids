#ifndef PHIANDNU_LOGFILE_INFORMATION_H
#define PHIANDNU_LOGFILE_INFORMATION_H

#include "Utilities\LogFileInformation\TestLogFileInformation\TestLogFileInformation.h"

#include <memory>
#include <vector>

class PhiAndNuTest;

class PhiAndNuInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr< PhiAndNuInformation> getNewInstance(unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation);

	std::string getOutput();
	void addTestGroup(std::vector< std::shared_ptr< PhiAndNuTest>> tests);

private:
	void fillMyData(std::vector< std::shared_ptr< PhiAndNuTest>> testGroup);

	PhiAndNuInformation() {};
	PhiAndNuInformation(unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation);

	std::vector< std::vector< std::shared_ptr< PhiAndNuTest>>> testGroups;
	unsigned int startTimeStepCalculation, endTimeStepCalculation;

	std::vector<int> lx;
	std::vector<int> lxForErase;
	std::vector<double> phiDiff;
	std::vector<double> nu, nuDiff;
	std::vector<double> orderOfAccuracyPhiDiff;
	std::vector<double> orderOfAccuracyNuDiff;
	std::vector<std::string> dataToCalc;
};
#endif