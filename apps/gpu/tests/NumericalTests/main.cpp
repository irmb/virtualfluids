#include <gmock/gmock.h>
#include <mpi.h>

#include "Utilities/ConfigFileReaderNT/ConfigFileReaderNT.h"
#include "Utilities/LogFileQueue/LogFileQueue.h"
#include "Utilities/NumericalTestFactory/NumericalTestFactoryImp.h"
#include "Utilities/TestQueue/TestQueue.h"
#include "Utilities/VirtualFluidSimulationFactory/VirtualFluidSimulationFactory.h"

// validation
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include <fstream>
#include <iostream>

static TestSuiteResult startNumericalTests(const std::string &configFile)
{
    auto configData = vf::gpu::tests::readConfigFile(configFile);

    std::shared_ptr<NumericalTestFactoryImp> numericalTestFactory = NumericalTestFactoryImp::getNewInstance(configData);

    std::shared_ptr<TestQueue> testQueue = numericalTestFactory->getTestQueue();
    std::shared_ptr<LogFileQueue> logFileQueue = numericalTestFactory->getLogFileQueue();

    auto result = testQueue->run();
    logFileQueue->writeLogFiles();

    return result;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    auto tests_passed = TestSuiteResult::FAILED;

    if (argc > 1)
        tests_passed = startNumericalTests(argv[1]);
    else
        std::cout << "Configuration file must be set!: lbmgm <config file>" << std::endl << std::flush;

    MPI_Finalize();

    return tests_passed;
}
