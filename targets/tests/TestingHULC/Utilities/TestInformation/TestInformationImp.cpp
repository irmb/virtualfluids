#include "TestInformationImp.h"

#include <iomanip>

std::shared_ptr<TestInformation> TestInformationImp::getNewInstance(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity, bool tgv, double u0TGV, double AmplitudeTGV, bool sw, double u0SW, double v0SW)
{
	return std::shared_ptr<TestInformation>(new TestInformationImp(numberOfTimeSteps, basisTimeStepLength, startStepCalculation, viscosity, tgv, u0TGV, AmplitudeTGV, sw, u0SW, v0SW));
}

TestInformationImp::TestInformationImp(int numberOfTimeSteps, int basisTimeStepLength, int startStepCalculation, double viscosity, bool tgv, double u0TGV, double amplitudeTGV, bool sw, double u0SW, double v0SW)
{
	makeCenterHead("Basic Information");
	oss << "NumberOfTimeSteps: " << numberOfTimeSteps << std::endl;
	oss << "BasisTimeStepLength: " << basisTimeStepLength << std::endl;
	oss << "StartStepCalculation: " << startStepCalculation << std::endl;
	oss << "Viscosity: " << viscosity << std::endl;
	oss << std::endl;
	
	if (tgv) {
		makeCenterHead("TaylorGreenVortex Information");
		oss << "u0: " << u0TGV << std::endl;
		oss << "Amplitude: " << amplitudeTGV << std::endl;
		oss << std::endl;
	}
	if (sw) {
		makeCenterHead("ShearWave Information");
		oss << "u0: " << u0SW << std::endl;
		oss << "v0: " << v0SW << std::endl;
		oss << std::endl;
	}
}

std::string TestInformationImp::getInformation()
{
	return oss.str();
}

void TestInformationImp::makeHastTags()
{
	oss << "#################################################" << std::endl;
}

void TestInformationImp::makeCenterHead(std::string output)
{
	makeHastTags();
	oss << "#" << std::setfill(' ') << std::right << std::setw(24 + output.length() / 2) << output << std::setw(24 - output.length() / 2) << "#" << std::endl;
	makeHastTags();
}
