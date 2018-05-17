#include "ShearWaveLogFileInformation.h"

std::shared_ptr<LogFileInformation> ShearWaveInformation::getNewInstance(double u0, double v0, std::vector< bool> tests, std::vector< real> l)
{
	return std::shared_ptr<LogFileInformation>(new ShearWaveInformation(u0, v0, tests, l));
}

std::string ShearWaveInformation::getOutput()
{
	makeCenterHead("ShearWave Information");
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)) {
			oss << "Lx:" << l.at(i) << std::endl;
			oss << "u0: " << u0 / (l.at(i) / l0) << std::endl;
			oss << "v0: " << v0 / (l.at(i) / l0) << std::endl;
			oss << std::endl;
		}
	}
	return oss.str();
}

ShearWaveInformation::ShearWaveInformation(double u0, double v0, std::vector< bool> tests, std::vector< real> l) : u0(u0), v0(v0), tests(tests), l(l)
{
	l0 = 32;
}