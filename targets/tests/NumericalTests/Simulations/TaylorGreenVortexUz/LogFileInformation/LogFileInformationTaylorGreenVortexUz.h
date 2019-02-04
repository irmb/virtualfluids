#ifndef LOGFILE_INFORMATION_TAYLOR_GREEN_UZ_H
#define LOGFILE_INFORMATION_TAYLOR_GREEN_UZ_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities\LogFileInformation\SimulationLogFileInformation\SimulationLogFileInformation.h"

#include <memory>
#include <vector>

struct TaylorGreenVortexUzParameterStruct;
struct GridInformationStruct;

class LogFileInformationTaylorGreenUz : public SimulationLogFileInformation
{
public:
	static std::shared_ptr<LogFileInformationTaylorGreenUz> getNewInstance(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
	
	std::string getOutput();
	std::string getFilePathExtensionOne();
	std::string getFilePathExtensionTwo();

private:
	LogFileInformationTaylorGreenUz() {};
	LogFileInformationTaylorGreenUz(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

	double uz;
	double amplitude;
	std::vector<double> lz;
	int l0;
};
#endif 