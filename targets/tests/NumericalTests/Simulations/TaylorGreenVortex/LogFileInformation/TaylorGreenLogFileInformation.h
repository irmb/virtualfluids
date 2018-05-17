#ifndef TAYLOR_GREEN_INFORMATION_H
#define TAYLOR_GREEN_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"

#include "LBM\LB.h"

#include <memory>
#include <vector>

class TaylorGreenInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(double u0, double amplitude, std::vector< bool> tests, std::vector< real> l);
	std::string getOutput();

private:
	TaylorGreenInformation() {};
	TaylorGreenInformation(double u0, double amplitude, std::vector< bool> tests, std::vector< real> l);

	double u0;
	double amplitude;
	std::vector< bool> tests;
	std::vector< real> l;
	int l0;
};
#endif 