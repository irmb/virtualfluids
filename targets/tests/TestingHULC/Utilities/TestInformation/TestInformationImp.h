#ifndef TESTINFORMATIONIMP_H
#define TESTINFORMATIONIMP_H

#include "TestInformation.h"

#include <sstream>
#include <memory>

class TestInformationImp : public TestInformation
{
public:
	static std::shared_ptr<TestInformation> getNewInstance(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity, bool tgv, double u0TGV, double AmplitudeTGV, bool sw, double u0SW, double v0SW);
	std::string getInformation();

protected:
	TestInformationImp(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity, bool tgv, double u0TGV, double AmplitudeTGV, bool sw, double u0SW, double v0SW);
	TestInformationImp() {};

private:
	void makeHastTags();
	void makeCenterHead(std::string head);

	std::ostringstream oss;

};
#endif // !TESTINFORMATIONIMP_H
