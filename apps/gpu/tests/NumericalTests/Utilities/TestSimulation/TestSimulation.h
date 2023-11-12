#ifndef TEST_SIMULATION_H
#define TEST_SIMULATION_H

#include <memory>
#include <string>

class InitialCondition;
class DataWriter;
class SimulationParameter;
class SimulationResults;
class SimulationObserver;
class TimeTracking;
class Parameter;

class TestSimulation
{
public:
    virtual ~TestSimulation() = default;
    virtual void run() = 0;
    virtual void makeSimulationHeadOutput() = 0;
    virtual void startPostProcessing() = 0;

    virtual std::shared_ptr<SimulationParameter> getSimulationParameter() = 0;
    virtual std::shared_ptr<TimeTracking> getTimeTracking() = 0;
    virtual void setParameter(std::shared_ptr<Parameter> para) = 0;
};
#endif