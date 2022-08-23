#include "TestImp.h"
#include <algorithm>

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/NumericalTestSimulation/NumericalTestSimulation.h"
#include "Utilities/PostProcessingStrategy/PostProcessingStrategy.h"
#include "Utilities/SimulationInfo/SimulationInfo.h"

void TestImp::run()
{
    for (size_t i = 0; i < simulations.size(); i++) {
        auto sim = simulations.at(i);
        auto simInfo = simInfos.at(i);

        // NOTE: Simulations can be in this vector multiple times
        // Therefore, we skip the simulation if it has been run already
        if (simulationRun.at(i))
            continue;
        sim->run();
    }
}

void TestImp::update()
{
    for (size_t i = 0; i < simulations.size(); i++) {
        if (simulationRun.at(i) == false) {
            switch (simulations.at(i)->getSimulationStatus()) {
                case executed:
                    simulationRun.at(i) = true;
                    postProStrategies.at(i)->evaluate();
                    break;
                case crashed:
                    simulationRun.at(i) = true;
                    testStatus = simulationCrashed;
                    break;
                case initialized:
                    simulationRun.at(i) = false;
                    break;
                default:
                    break;
            }
        }
    }
    if (CheckAllSimulationRun()) {
        if (testStatus != simulationCrashed)
            evaluate();
        else
            makeConsoleOutput();
    }
}

void TestImp::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo,
                            std::shared_ptr<PostProcessingStrategy> postProStrategy)
{
    simulations.push_back(sim);
    simInfos.push_back(simInfo);
    postProStrategies.push_back(postProStrategy);
    simulationRun.push_back(false);
}

TestStatus TestImp::getTestStatus()
{
    return testStatus;
}

void TestImp::makeConsoleOutput()
{
    switch (testStatus) {
        case passed:
            colorOutput->makeTestOutput(buildTestOutput(), testStatus);
            break;
        case failed:
            colorOutput->makeTestOutput(buildTestOutput(), testStatus);
            break;
        case test_error:
            colorOutput->makeTestOutput(buildErrorTestOutput(), testStatus);
            break;
        case simulationCrashed:
            colorOutput->makeTestOutput(buildSimulationFailedTestOutput(), testStatus);
            break;
        default:
            break;
    }
}

TestImp::TestImp(std::shared_ptr<ColorConsoleOutput> colorOutput) : colorOutput(colorOutput)
{
    simulationRun.resize(0);
    simulations.resize(0);
    simInfos.resize(0);
}

bool TestImp::CheckAllSimulationRun()
{
    return std::all_of(simulationRun.begin(), simulationRun.end(), [](bool run) { return run; });
}

std::vector<std::string> TestImp::buildSimulationFailedTestOutput()
{
    std::vector<std::string> output = buildBasicTestOutput();
    std::ostringstream oss;

    oss << "Simulation crashed!";
    output.push_back(oss.str());
    oss.str(std::string());

    return output;
}
