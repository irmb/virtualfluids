#include "ShearWaveLogFileInformation.h"

std::shared_ptr<ShearWaveInformation> ShearWaveInformation::getNewInstance(double u0, double v0, std::vector< bool> tests, std::vector< real> l, int l0)
{
	return std::shared_ptr<ShearWaveInformation>(new ShearWaveInformation(u0, v0, tests, l, l0));
}

std::string ShearWaveInformation::getOutput()
{
	makeCenterHead("ShearWave Information");
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)) {
			oss << "Lx=" << l.at(i) << std::endl;
			oss << "l0=" << l0 << std::endl;
			oss << "u0=" << u0 / (l.at(i) / l0) << std::endl;
			oss << "v0=" << v0 / (l.at(i) / l0) << std::endl;
			oss << std::endl;
		}
	}
	return oss.str();
}

std::string ShearWaveInformation::getFilePathExtensionOne()
{
	std::ostringstream oss;
	oss << "ShearWave\\";
	return oss.str();
}

std::string ShearWaveInformation::getFilePathExtensionTwo()
{
	std::ostringstream oss;
	oss << "u0_" << u0 << "_v0_" << v0 << "\\";
	return oss.str();
}

ShearWaveInformation::ShearWaveInformation(double u0, double v0, std::vector< bool> tests, std::vector< real> l, int l0) : u0(u0), v0(v0), tests(tests), l(l), l0(l0)
{

}