#include "L2NormLogFileInformation.h"

#include "Tests\L2NormTest\L2NormTest.h"

#include <iomanip>
#include <sstream>

std::shared_ptr<L2NormInformation> L2NormInformation::getNewInstance(std::vector<std::shared_ptr<L2NormTest>> tests, unsigned int basicTimeStep, unsigned int divergentTimeStep)
{
	return std::shared_ptr<L2NormInformation>(new L2NormInformation(tests, basicTimeStep, divergentTimeStep));
}

std::string L2NormInformation::getOutput()
{
	std::ostringstream headName;
	headName << tests.at(0)->getSimulationName() << " L2Norm Test";
	makeCenterHead(headName.str());

	oss << "BasicTimeStep_L2Norm=" << basicTimeStep << std::endl;
	oss << "DivergentTimeStep_L2Norm=" << divergentTimeStep << std::endl;
	oss << std::endl;

	for (int i = 0; i < tests.size(); i++)
		oss << tests.at(i)->getLogFileOutput();

	return oss.str();
}

L2NormInformation::L2NormInformation(std::vector<std::shared_ptr<L2NormTest>> tests, unsigned int basicTimeStep, unsigned int divergentTimeStep) : tests(tests), basicTimeStep(basicTimeStep), divergentTimeStep(divergentTimeStep)
{

}