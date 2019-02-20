#include <mpi.h>
#include <gmock/gmock.h>

#include "Utilities/ConfigFileReader/ConfigFileReader.h"
#include "Utilities/LogFileQueue/LogFileQueue.h"
#include "Utilities/NumericalTestFactory/NumericalTestFactoryImp.h"
#include "Utilities/TestQueue/TestQueue.h"
#include "Utilities/VirtualFluidSimulation/VirtualFluidSimulation.h"
#include "Utilities/VirtualFluidSimulationFactory/VirtualFluidSimulationFactoryImp.h"

static void startNumericalTests(const std::string &configFile)
{
	std::shared_ptr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance(configFile);

	std::shared_ptr<NumericalTestFactoryImp> numericalTestFactory = NumericalTestFactoryImp::getNewInstance(configReader->getConfigData());

	std::vector<std::shared_ptr<TestSimulation> > testSim = numericalTestFactory->getTestSimulations();
	std::shared_ptr<TestQueue> testQueue = numericalTestFactory->getTestQueue();
	std::shared_ptr<LogFileQueue> logFileQueue = numericalTestFactory->getLogFileQueue();

	std::shared_ptr<VirtualFluidSimulationFactory> factory = VirtualFluidSimulationFactoryImp::getNewInstance();
	std::vector<std::shared_ptr<VirtualFluidSimulation> > vfSimulations = factory->makeVirtualFluidSimulations(testSim);

	for (int i = 0; i < vfSimulations.size(); i++)
		vfSimulations.at(i)->run();

	testQueue->makeFinalOutput();
	logFileQueue->writeLogFiles();
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	if (argc > 1)
		startNumericalTests(argv[1]);
	else
		std::cout << "Configuration file must be set!: lbmgm <config file>" << std::endl << std::flush;

    MPI_Finalize();

	return 0;
}
