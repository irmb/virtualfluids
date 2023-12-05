#ifndef SIMULATION_PARAMETER_IMP_H
#define SIMULATION_PARAMETER_IMP_H

#include "SimulationParameter.h"

#include "Calculation/Calculation.h"
struct GridInformationStruct;
struct BasicSimulationParameterStruct;

class SimulationParameterImp : public SimulationParameter
{
public:
    double getViscosity();
    std::string getGridPath();
    std::string getFilePath();
    unsigned int getNumberOfGridLevels();
    unsigned int getEndTime();
    unsigned int getTimeStepLength();
    unsigned int getLx();
    unsigned int getLz();
    unsigned int getL0();
    std::vector<unsigned int> getDevices();
    double getMaxVelocity();
    
    std::shared_ptr<KernelConfiguration> getKernelConfiguration();

protected:
    SimulationParameterImp() {};
    SimulationParameterImp(std::string kernelName, double viscosity, std::shared_ptr<BasicSimulationParameterStruct> basicSimPara, std::shared_ptr<GridInformationStruct> gridInfo);

    void generateFileDirectionInMyStystem(std::string filePath);

    unsigned int timeStepLength;
    std::string filePath;
    double maxVelocity;
    real lx, l0, lz;

private:
    real viscosity;
    unsigned int numberOfTimeSteps, basisTimeStepLength;
    std::string gridPath;
    std::vector<unsigned int> devices;
    unsigned int maxLevel, numberOfGridLevels;
    std::shared_ptr<KernelConfiguration> kernelConfig;
};

#endif
