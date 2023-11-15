#include "ShearWaveSimulationInfo.h"

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#include <sstream>

std::shared_ptr<ShearWaveSimulationInfo> ShearWaveSimulationInfo::getNewInstance(int simID, std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
{
    return std::shared_ptr<ShearWaveSimulationInfo>(new ShearWaveSimulationInfo(simID, kernel,viscosity, simParaStruct, gridInfoStruct, numberOfSimulations));
}

ShearWaveSimulationInfo::ShearWaveSimulationInfo(int simID, std::string kernel, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct, int numberOfSimulations)
    : SimulationInfoImp(simID, kernel, viscosity, gridInfoStruct->lx, numberOfSimulations, "ShearWave", simParaStruct->dataToCalcTests)
{
    std::ostringstream oss;
    oss << " ux: " << simParaStruct->ux / (gridInfoStruct->lx / simParaStruct->l0) << " uz: " << simParaStruct->uz / (gridInfoStruct->lx / simParaStruct->l0);
    this->simulationParameterString = oss.str();

}