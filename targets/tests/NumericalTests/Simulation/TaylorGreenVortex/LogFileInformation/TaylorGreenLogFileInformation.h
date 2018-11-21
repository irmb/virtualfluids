#ifndef TAYLOR_GREEN_INFORMATION_H
#define TAYLOR_GREEN_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"

#include <memory>
#include <vector>

class TaylorGreenInformation : public LogFileInformationImp, public SimulationLogFileInformation
{
public:
	static std::shared_ptr<TaylorGreenInformation> getNewInstance(double u0, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);
	
	std::string getOutput();
	std::string getFilePathExtension();

private:
	TaylorGreenInformation() {};
	TaylorGreenInformation(double u0, double amplitude, std::vector< bool> tests, std::vector< double> l, int l0);

	double u0;
	double amplitude;
	std::vector< bool> tests;
	std::vector< double> l;
	int l0;
};
#endif 