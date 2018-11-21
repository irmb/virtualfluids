#include "TaylorGreenLogFileInformation.h"

std::shared_ptr<TaylorGreenInformation> TaylorGreenInformation::getNewInstance(double u0, double amplitude, std::vector< bool> tests, std::vector<double> l, int l0)
{
	return std::shared_ptr<TaylorGreenInformation>(new TaylorGreenInformation(u0, amplitude, tests, l, l0));
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

std::string TaylorGreenInformation::getFilePathExtension()
{
	std::ostringstream oss;
	oss <<"TaylorGreenVortex\\u0_ " << u0 << "_Amplitude_ " << amplitude;
	return oss.str();
}

TaylorGreenInformation::TaylorGreenInformation(double u0, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0) : u0(u0), amplitude(amplitude), tests(tests), l(l), l0(l0)
{
}
