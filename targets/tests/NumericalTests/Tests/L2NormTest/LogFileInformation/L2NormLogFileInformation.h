#ifndef L2NORM_LOGFILE_INFORMATION_H
#define L2NORM_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"

#include <memory>
#include <vector>

class L2NormTest;
struct L2NormTestParameterStruct;

class L2NormInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr<L2NormInformation> getNewInstance(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests);

	std::string getOutput();

private:
	L2NormInformation() {};
	L2NormInformation(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests);

	std::vector<std::shared_ptr<L2NormTest> > tests;

	unsigned int basicTimeStep, divergentTimeStep;
	std::vector<std::string> dataToCalc;
	std::vector<std::string> normalizeData;
};
#endif