#ifndef SHEAR_WAVE_INFORMATION_H
#define SHEAR_WAVE_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"

#include "LBM\LB.h"

#include <memory>
#include <vector>

class ShearWaveInformation : public LogFileInformationImp, public SimulationLogFileInformation
{
public:
	static std::shared_ptr<ShearWaveInformation> getNewInstance(double u0, double v0, std::vector< bool> tests, std::vector< real> l, int l0);

	std::string getOutput();
	std::string getFilePathExtension();

private:
	ShearWaveInformation() {};
	ShearWaveInformation(double u0, double v0, std::vector< bool> tests, std::vector< real> l, int l0);

	double u0;
	double v0;
	std::vector< bool> tests;
	std::vector< real> l;
	int l0;
};
#endif