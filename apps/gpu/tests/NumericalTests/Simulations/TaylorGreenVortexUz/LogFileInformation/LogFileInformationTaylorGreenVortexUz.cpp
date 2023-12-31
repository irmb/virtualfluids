#include "LogFileInformationTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"

std::shared_ptr<LogFileInformationTaylorGreenUz> LogFileInformationTaylorGreenUz::getNewInstance(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    return std::shared_ptr<LogFileInformationTaylorGreenUz>(new LogFileInformationTaylorGreenUz(simParaStruct, gridInfoStruct));
}

std::string LogFileInformationTaylorGreenUz::getOutput()
{
    makeCenterHead("TaylorGreenVortex V0 Information");
    oss << "SimulationName=TaylorGreenVortexUz" << std::endl;
    oss << "Lx=\"";
    for (int i = 0; i < lx.size(); i++) {
        oss << lx.at(i);
        if (i < lx.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl << std::endl;
    }

    for (int i = 0; i < lz.size(); i++) {
        oss << "l0_" << lx.at(i) << "=" << l0 << std::endl;
        oss << "uz_" << lx.at(i) << "=" << uz / (lz.at(i) / l0) << std::endl;
        oss << "Amplitude_" << lx.at(i) << "=" << amplitude / (lz.at(i) / l0) << std::endl;
        oss << std::endl;
    }
    
    return oss.str();
}

std::vector<std::string> LogFileInformationTaylorGreenUz::getFilePathExtension()
{
    std::vector<std::string> myFilePath;
    myFilePath.push_back("TaylorGreenVortexUz");
    std::ostringstream oss;
    oss << "uz_" << uz << "_Amplitude_" << amplitude;
    myFilePath.push_back(oss.str());
    return myFilePath;
}

LogFileInformationTaylorGreenUz::LogFileInformationTaylorGreenUz(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    this->uz = simParaStruct->uz;
    this->amplitude = simParaStruct->amplitude;
    this->l0 = simParaStruct->l0;

    for (int i = 0; i < gridInfoStruct.size(); i++)
        lz.push_back(gridInfoStruct.at(i)->lz);
    for (int i = 0; i < gridInfoStruct.size(); i++)
        lx.push_back(gridInfoStruct.at(i)->lx);
}
