#include "ShearWaveLogFileInformation.h"

#include "Simulations\ShearWave\ShearWaveParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

std::shared_ptr<ShearWaveInformation> ShearWaveInformation::getNewInstance(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
	return std::shared_ptr<ShearWaveInformation>(new ShearWaveInformation(simParaStruct, gridInfoStruct));
}

std::string ShearWaveInformation::getOutput()
{
	makeCenterHead("ShearWave Information");
	for (int i = 0; i < lx.size(); i++) {
			oss << "Lx=" << lx.at(i) << std::endl;
			oss << "l0=" << l0 << std::endl;
			oss << "ux=" << ux / (lx.at(i) / l0) << std::endl;
			oss << "uz=" << uz / (lx.at(i) / l0) << std::endl;
			oss << std::endl;
	}
	return oss.str();
}

std::string ShearWaveInformation::getFilePathExtensionOne()
{
	std::ostringstream oss;
	oss << "ShearWave\\";
	return oss.str();
}

std::string ShearWaveInformation::getFilePathExtensionTwo()
{
	std::ostringstream oss;
	oss << "ux_" << ux << "_uz_" << uz << "\\";
	return oss.str();
}

ShearWaveInformation::ShearWaveInformation(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
	this->ux = simParaStruct->ux;
	this->uz = simParaStruct->uz;
	this->l0 = simParaStruct->l0;

	for (int i = 0; i < gridInfoStruct.size(); i++)
		lx.push_back(gridInfoStruct.at(i)->lx);
}