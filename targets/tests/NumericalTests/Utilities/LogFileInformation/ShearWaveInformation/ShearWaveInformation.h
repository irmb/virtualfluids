#ifndef SHEAR_WAVE_INFORMATION_H
#define SHEAR_WAVE_INFORMATION_H

#include "../LogFileInformationImp.h"

#include <memory>

class ShearWaveInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(double u0, double v0);
	std::string getOutput();

private:
	ShearWaveInformation() {};
	ShearWaveInformation(double u0, double v0);

	double u0;
	double v0;
};
#endif