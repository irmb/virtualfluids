#include "L2NormLogFileInformationBetweenKernels.h"

#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernels.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"

#include <iomanip>
#include <sstream>

std::shared_ptr<L2NormBetweenKernelsInformation> L2NormBetweenKernelsInformation::getNewInstance(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara)
{
	return std::shared_ptr<L2NormBetweenKernelsInformation>(new L2NormBetweenKernelsInformation(tests, testPara));
}

std::string L2NormBetweenKernelsInformation::getOutput()
{
	std::ostringstream headName;
	headName << tests.at(0)->getSimulationName() << "L2Norm Test Between Kernels";
	makeCenterHead(headName.str());

	oss << "BasicKernel=" << basicKernel << std::endl;
	oss << "DataToCalculate_L2Norm_BK=\"";
	for (int i = 0; i < dataToCalc.size(); i++)
		oss << dataToCalc.at(i) << " ";
	deleteLastCharInOss();
	oss << "\"" << std::endl;
	oss << "TimeSteps=\""; 
	for (int i = 0; i < timeSteps.size(); i++)
		oss << timeSteps.at(i) << " ";
	deleteLastCharInOss();
	oss << "\""<< std::endl << std::endl;

	for (int i = 0; i < tests.size(); i++)
		oss << tests.at(i)->getLogFileOutput();

	return oss.str();
}

L2NormBetweenKernelsInformation::L2NormBetweenKernelsInformation(std::vector<std::shared_ptr<L2NormTestBetweenKernels> > tests, std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara)
	: tests(tests)
{
	basicKernel = testPara->basicKernel;
	timeSteps = testPara->timeSteps;
	dataToCalc = testPara->basicTestParameter->dataToCalc;
}

void L2NormBetweenKernelsInformation::deleteLastCharInOss()
{
	std::string myString = oss.str().substr(0, oss.str().size() - 1);
	oss.str("");
	oss.clear();
	oss << myString;
}
