#ifndef LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H
#define LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"

#include <memory>
#include <vector>

class LogFileInformationTaylorGreenUx : public LogFileInformationImp, public SimulationLogFileInformation
{
public:
	static std::shared_ptr<LogFileInformationTaylorGreenUx> getNewInstance(double ux, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);
	
	std::string getOutput();
	std::string getFilePathExtensionOne();
	std::string getFilePathExtensionTwo();

private:
	LogFileInformationTaylorGreenUx() {};
	LogFileInformationTaylorGreenUx(double ux, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);

	double ux;
	double amplitude;
	std::vector< bool> tests;
	std::vector< double> l;
	int l0;
};
#endif 