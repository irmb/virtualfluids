#include "TaylorGreenLogFileInformation.h"

std::shared_ptr<LogFileInformation> TaylorGreenInformation::getNewInstance(double u0, double amplitude, std::vector< bool> tests, std::vector< real> l)
{
	return std::shared_ptr<LogFileInformation>(new TaylorGreenInformation(u0, amplitude, tests, l));
}

std::string TaylorGreenInformation::getOutput()
{
	makeCenterHead("TaylorGreenVortex Information");
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)) {
			oss << "Lx:" << l.at(i) << std::endl;
			oss << "u0: " << u0 / (l.at(i) / l0) << std::endl;
			oss << "Amplitude: " << amplitude / (l.at(i) / l0) << std::endl;
			oss << std::endl;
		}
	}
	
	return oss.str();
}

TaylorGreenInformation::TaylorGreenInformation(double u0, double amplitude, std::vector< bool> tests, std::vector< real> l) : u0(u0), amplitude(amplitude), tests(tests), l(l)
{
	l0 = 32;
}
