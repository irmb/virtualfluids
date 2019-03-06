#include "L2NormLogFileInformationBetweenKernels.h"

#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernels.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"

#include <iomanip>
#include <sstream>

std::shared_ptr<L2NormBetweenKernelsInformation> L2NormBetweenKernelsInformation::getNewInstance(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
{
	return std::shared_ptr<L2NormBetweenKernelsInformation>(new L2NormBetweenKernelsInformation(tests, testPara, dataToCalcTests));
}

std::string L2NormBetweenKernelsInformation::getOutput()
{
	std::ostringstream headName;
	headName << "L2Norm Test Between Kernels";
	makeCenterHead(headName.str());

	oss << "BasicKernel_L2Norm_BK=" << basicKernel << std::endl;
	oss << "DataToCalculate_L2Norm_BK=\"";
	for (int i = 0; i < dataToCalc.size(); i++)
		oss << dataToCalc.at(i) << " ";
	deleteLastCharInOss();
	oss << "\"" << std::endl;
	oss << "TimeSteps_L2Norm_BK=\""; 
	for (int i = 0; i < timeSteps.size(); i++)
		oss << timeSteps.at(i) << " ";
	deleteLastCharInOss();
	oss << "\""<< std::endl << std::endl;
	oss << "NormalizeWith_L2Norm_BK=\"";
	for (int i = 0; i < normalizeData.size(); i++)
		oss << normalizeData.at(i) << " ";
	deleteLastCharInOss();
	oss << "\"" << std::endl << std::endl;

	std::ostringstream failMessage;
	failMessage << "FailTests_L2Norm_BK=\"";
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)->getTestStatus() == passed || tests.at(i)->getTestStatus() == failed)
			oss << tests.at(i)->getLogFileOutput();
		if (tests.at(i)->getTestStatus() == error || tests.at(i)->getTestStatus() == simulationCrashed)
			failMessage << tests.at(i)->getErrorLogFileOutput() << " ";
	}
	std::string fail = failMessage.str();
	if (fail.back() == ' ')
		fail = fail.substr(0, fail.size() - 1);
	failMessage.str(std::string());
	failMessage << fail << "\"";
	oss << failMessage.str() << std::endl << std::endl;

	return oss.str();
}

L2NormBetweenKernelsInformation::L2NormBetweenKernelsInformation(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::string> dataToCalcTests)
	: tests(tests), dataToCalc(dataToCalcTests)
{
	basicKernel = testPara->basicKernel;
	timeSteps = testPara->timeSteps;
	normalizeData = testPara->normalizeData;
}

void L2NormBetweenKernelsInformation::deleteLastCharInOss()
{
	std::string myString = oss.str().substr(0, oss.str().size() - 1);
	oss.str("");
	oss.clear();
	oss << myString;
}
