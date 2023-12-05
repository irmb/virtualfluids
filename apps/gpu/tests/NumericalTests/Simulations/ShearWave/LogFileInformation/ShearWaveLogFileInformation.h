#ifndef SHEAR_WAVE_INFORMATION_H
#define SHEAR_WAVE_INFORMATION_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities/LogFileInformation/SimulationLogFileInformation/SimulationLogFileInformation.h"

#include "Calculation/Calculation.h"

#include <memory>
#include <vector>

struct ShearWaveParameterStruct;
struct GridInformationStruct;

class ShearWaveInformation : public SimulationLogFileInformation
{
public:
    static std::shared_ptr<ShearWaveInformation> getNewInstance(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

    std::string getOutput();
    std::vector<std::string> getFilePathExtension();

private:
    ShearWaveInformation() {};
    ShearWaveInformation(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

    double ux;
    double uz;
    std::vector<real> lx;
    int l0;
};
#endif