#include "LogFileInformationTaylorGreenVortexUx.h"

#include "Simulations\TaylorGreenVortexUx\TaylorGreenVortexUxParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

std::shared_ptr<LogFileInformationTaylorGreenUx> LogFileInformationTaylorGreenUx::getNewInstance(std::shared_ptr< TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector< std::shared_ptr< GridInformationStruct> > gridInfoStruct)
{
	return std::shared_ptr<LogFileInformationTaylorGreenUx>(new LogFileInformationTaylorGreenUx(simParaStruct, gridInfoStruct));
}

std::string LogFileInformationTaylorGreenUx::getOutput()
{
	makeCenterHead("TaylorGreenVortex U0 Information");
	for (int i = 0; i < lx.size(); i++) {
		oss << "Lx=" << lx.at(i) << std::endl;
		oss << "ux=" << ux / (lx.at(i) / l0) << std::endl;
		oss << "Amplitude= " << amplitude / (lx.at(i) / l0) << std::endl;
		oss << "l0=" << l0 << std::endl;
		oss << std::endl;
	}
	
	return oss.str();
}

std::string LogFileInformationTaylorGreenUx::getFilePathExtensionTwo()
{
	std::ostringstream oss;
	oss << "ux_ " << ux << "_Amplitude_ " << amplitude << "\\";
	return oss.str();
}

std::string LogFileInformationTaylorGreenUx::getFilePathExtensionOne()
{
	std::ostringstream oss;
	oss <<"TaylorGreenVortexUx\\";
	return oss.str();
}

LogFileInformationTaylorGreenUx::LogFileInformationTaylorGreenUx(std::shared_ptr< TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector< std::shared_ptr< GridInformationStruct> > gridInfoStruct)
{
	this->ux = simParaStruct->ux;
	this->amplitude = simParaStruct->amplitude;
	this->l0 = simParaStruct->l0;

	for(int i = 0; i < gridInfoStruct.size(); i++)
		lx.push_back(gridInfoStruct.at(i)->lx);
}
