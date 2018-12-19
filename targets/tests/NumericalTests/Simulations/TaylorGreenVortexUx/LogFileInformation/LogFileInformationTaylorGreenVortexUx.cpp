#include "LogFileInformationTaylorGreenVortexUx.h"

std::shared_ptr<LogFileInformationTaylorGreenUx> LogFileInformationTaylorGreenUx::getNewInstance(double ux, double amplitude, std::vector< bool> tests, std::vector<double> l, int l0)
{
	return std::shared_ptr<LogFileInformationTaylorGreenUx>(new LogFileInformationTaylorGreenUx(ux, amplitude, tests, l, l0));
}

std::string LogFileInformationTaylorGreenUx::getOutput()
{
	makeCenterHead("TaylorGreenVortex U0 Information");
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)) {
			oss << "Lx=" << l.at(i) << std::endl;
			oss << "ux=" << ux / (l.at(i) / l0) << std::endl;
			oss << "Amplitude= " << amplitude / (l.at(i) / l0) << std::endl;
			oss << "l0=" << l0 << std::endl;
			oss << std::endl;
		}
	}
	
	return oss.str();
}

std::string LogFileInformationTaylorGreenUx::getFilePathExtensionTwo()
{
	std::ostringstream oss;
	oss << "ux_ " << ux << "_Amplitude_ " << amplitude << "\\";
	return oss.str();
}

std::string LogFileInformationTaylorGreenUx::getFilePathExtensionOne()
{
	std::ostringstream oss;
	oss <<"TaylorGreenVortexUx\\";
	return oss.str();
}

LogFileInformationTaylorGreenUx::LogFileInformationTaylorGreenUx(double ux, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0) : ux(ux), amplitude(amplitude), tests(tests), l(l), l0(l0)
{
}
