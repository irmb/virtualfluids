#include "TaylorGreenInformation.h"

std::shared_ptr<LogFileInformation> TaylorGreenInformation::getNewInstance(double u0, double amplitude)
{
	return std::shared_ptr<LogFileInformation>(new TaylorGreenInformation(u0, amplitude));
}

std::string TaylorGreenInformation::getOutput()
{
	makeCenterHead("TaylorGreenVortex Information");
	oss << "u0: " << u0 << std::endl;
	oss << "Amplitude: " << amplitude << std::endl;
	oss << std::endl;
	return oss.str();
}

TaylorGreenInformation::TaylorGreenInformation(double u0, double amplitude) : u0(u0), amplitude(amplitude)
{
}
