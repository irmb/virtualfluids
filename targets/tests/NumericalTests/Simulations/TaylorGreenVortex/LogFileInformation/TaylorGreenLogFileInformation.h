#ifndef TAYLOR_GREEN_INFORMATION_H
#define TAYLOR_GREEN_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"

#include <memory>

class TaylorGreenInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(double u0, double amplitude);
	std::string getOutput();

private:
	TaylorGreenInformation() {};
	TaylorGreenInformation(double u0, double amplitude);

	double u0;
	double amplitude;
};
#endif 