#include "ShearWaveLogFileInformation.h"

std::shared_ptr<LogFileInformation> ShearWaveInformation::getNewInstance(double u0, double v0)
{
	return std::shared_ptr<LogFileInformation>(new ShearWaveInformation(u0,v0));
}

std::string ShearWaveInformation::getOutput()
{
	makeCenterHead("ShearWave Information");
	oss << "u0: " << u0 << std::endl;
	oss << "v0: " << v0 << std::endl;
	oss << std::endl;
	return oss.str();
}

ShearWaveInformation::ShearWaveInformation(double u0, double v0) : u0(u0), v0(v0)
{
}