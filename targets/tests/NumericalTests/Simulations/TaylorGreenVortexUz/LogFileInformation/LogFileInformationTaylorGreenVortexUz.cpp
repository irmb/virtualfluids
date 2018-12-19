#include "LogFileInformationTaylorGreenVortexUz.h"

std::shared_ptr<LogFileInformationTaylorGreenUz> LogFileInformationTaylorGreenUz::getNewInstance(double uz, double amplitude, std::vector< bool> tests, std::vector<double> l, int l0)
{
	return std::shared_ptr<LogFileInformationTaylorGreenUz>(new LogFileInformationTaylorGreenUz(uz, amplitude, tests, l, l0));
}

std::string LogFileInformationTaylorGreenUz::getOutput()
{
	makeCenterHead("TaylorGreenVortex V0 Information");
	for (int i = 0; i < tests.size(); i++) {
		if (tests.at(i)) {
			oss << "Lx=" << l.at(i) << std::endl;
			oss << "l0=" << l0 << std::endl;
			oss << "uz=" << uz / (l.at(i) / l0) << std::endl;
			oss << "Amplitude=" << amplitude / (l.at(i) / l0) << std::endl;
			oss << std::endl;
		}
	}
	
	return oss.str();
}

std::string LogFileInformationTaylorGreenUz::getFilePathExtensionOne()
{
	std::ostringstream oss;
	oss << "TaylorGreenVortexUz\\";
	return oss.str();
}

std::string LogFileInformationTaylorGreenUz::getFilePathExtensionTwo()
{
	std::ostringstream oss;
	oss << "uz_ " << uz << "_Amplitude_ " << amplitude << "\\";
	return oss.str();
}

LogFileInformationTaylorGreenUz::LogFileInformationTaylorGreenUz(double uz, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0) : uz(uz), amplitude(amplitude), tests(tests), l(l), l0(l0)
{
}
