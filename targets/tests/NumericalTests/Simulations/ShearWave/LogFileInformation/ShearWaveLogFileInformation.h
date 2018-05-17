#ifndef SHEAR_WAVE_INFORMATION_H
#define SHEAR_WAVE_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"

#include "LBM\LB.h"

#include <memory>
#include <vector>

class ShearWaveInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(double u0, double v0, std::vector< bool> tests, std::vector< real> l);
	std::string getOutput();

private:
	ShearWaveInformation() {};
	ShearWaveInformation(double u0, double v0, std::vector< bool> tests, std::vector< real> l);

	double u0;
	double v0;
	std::vector< bool> tests;
	std::vector< real> l;
	int l0;
};
#endif