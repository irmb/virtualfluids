#ifndef PHIANDNU_LOGFILE_INFORMATION_H
#define PHIANDNU_LOGFILE_INFORMATION_H

#include "Utilities\LogFileInformation\TestLogFileInformation\TestLogFileInformation.h"

#include <memory>
#include <vector>

class PhiAndNuTest;

class PhiAndNuInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr< PhiAndNuInformation> getNewInstance(std::vector< std::shared_ptr< PhiAndNuTest>> tests, unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation);

	std::string getOutput();

private:
	PhiAndNuInformation() {};
	PhiAndNuInformation(std::vector< std::shared_ptr< PhiAndNuTest>> tests, unsigned int startTimeStepCalculation, unsigned int endTimeStepCalculation);

	std::vector< std::shared_ptr< PhiAndNuTest>> tests;
	unsigned int startTimeStepCalculation, endTimeStepCalculation;
};
#endif