#ifndef LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H
#define LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities/LogFileInformation/SimulationLogFileInformation/SimulationLogFileInformation.h"

#include <memory>
#include <vector>

struct TaylorGreenVortexUxParameterStruct;
struct GridInformationStruct;

class LogFileInformationTaylorGreenUx : public SimulationLogFileInformation
{
public:
	static std::shared_ptr<LogFileInformationTaylorGreenUx> getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
	
	std::string getOutput();
	std::vector<std::string> getFilePathExtension();
	

private:
	LogFileInformationTaylorGreenUx() {};
	LogFileInformationTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

	double ux;
	double amplitude;
	std::vector<double> lx;
	int l0;
};
#endif 