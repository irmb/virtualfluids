#ifndef L2NORM_BETWEEN_KERNELS_LOGFILE_INFORMATION_H
#define L2NORM_BETWEEN_KERNELS_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"

#include <memory>
#include <vector>

class L2NormTestBetweenKernels;
struct L2NormTestBetweenKernelsParameterStruct;

class L2NormBetweenKernelsInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr<L2NormBetweenKernelsInformation> getNewInstance(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::string> dataToCalcTests);

	std::string getOutput();

private:
	L2NormBetweenKernelsInformation() {};
	L2NormBetweenKernelsInformation(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::string> dataToCalcTests);

	void deleteLastCharInOss();

	std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests;

	std::string basicKernel;
	std::vector<int> timeSteps;
	std::vector<std::string> dataToCalc;
	std::vector<std::string> normalizeData;
};
#endif