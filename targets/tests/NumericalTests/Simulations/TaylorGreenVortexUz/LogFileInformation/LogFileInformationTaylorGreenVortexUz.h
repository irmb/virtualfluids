#ifndef LOGFILE_INFORMATION_TAYLOR_GREEN_UZ_H
#define LOGFILE_INFORMATION_TAYLOR_GREEN_UZ_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"

#include <memory>
#include <vector>

class LogFileInformationTaylorGreenUz : public LogFileInformationImp, public SimulationLogFileInformation
{
public:
	static std::shared_ptr<LogFileInformationTaylorGreenUz> getNewInstance(double uz, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);
	
	std::string getOutput();
	std::string getFilePathExtensionOne();
	std::string getFilePathExtensionTwo();

private:
	LogFileInformationTaylorGreenUz() {};
	LogFileInformationTaylorGreenUz(double uz, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);

	double uz;
	double amplitude;
	std::vector< bool> tests;
	std::vector< double> l;
	int l0;
};
#endif 