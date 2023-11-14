#ifndef CONFIG_DATA_STRUCT_H
#define CONFIG_DATA_STRUCT_H

#include <memory>
#include <string>
#include <vector>

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"
#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"
#include "Tests/L2NormTest/L2NormTestParameterStruct.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"
#include "Tests/NyTest/NyTestParameterStruct.h"
#include "Tests/PhiTest/PhiTestParameterStruct.h"
#include "Utilities/Structs/BasicSimulationParameterStruct.h"
#include "Utilities/Structs/VectorWriterInformationStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"
#include "Utilities/Structs/LogFileParameterStruct.h"

struct ConfigDataStruct
{
    std::vector<double> viscosity;
    std::vector<std::string> kernelsToTest;

    std::vector<std::shared_ptr<TaylorGreenVortexUxParameterStruct> > taylorGreenVortexUxParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > taylorGreenVortexUxGridInformation;

    std::vector<std::shared_ptr<TaylorGreenVortexUzParameterStruct> > taylorGreenVortexUzParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > taylorGreenVortexUzGridInformation;

    std::vector<std::shared_ptr<ShearWaveParameterStruct> > shearWaveParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > shearWaveGridInformation;

    


    bool writeAnalyticalToVTK;
    unsigned int ySliceForCalculation;
    
    std::string logFilePath;

    int numberOfSimulations;

    std::shared_ptr<PhiTestParameterStruct> phiTestParameter;
    std::shared_ptr<NyTestParameterStruct> nyTestParameter;
    std::shared_ptr<L2NormTestParameterStruct> l2NormTestParameter;
    std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> l2NormTestBetweenKernelsParameter;

    std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo;

    std::shared_ptr<LogFileParameterStruct> logFilePara;
};

#endif 