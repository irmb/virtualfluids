#include <gmock/gmock.h>
#include "mpi.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"

#include "Utilities/ConfigFileReader/ConfigFileReader.h"
#include "Utilities\VirtualFluidSimulation\VirtualFluidSimulation.h"
#include "Utilities\VirtualFluidSimulationFactory\VirtualFluidSimulationFactoryImp.h"
#include "Utilities/Calculator/Calculator.h"
#include "Utilities/TestInformation/TestInformation.h"

static void startNumericalTests(const std::string &configFile)
{
	std::shared_ptr< ConfigFileReader> configReader = ConfigFileReader::getNewInstance();
	configReader->readConfigFile(configFile);

	std::vector< std::shared_ptr< SimulationParameter> > simPara = configReader->getSimulationParameter();
	std::shared_ptr< TestInformation> testInfo = configReader->getTestInformation();

	std::shared_ptr< VirtualFluidSimulationFactory> factory = VirtualFluidSimulationFactoryImp::getNewInstance();
	std::vector< std::shared_ptr< VirtualFluidSimulation> > vfSimulations = factory->makeVirtualFluidSimulations(simPara);

	for (int i = 0; i < vfSimulations.size(); i++)
	{
		testInfo->makeSimulationHeadOutput(i);
		testInfo->setSimulationStartTime(i);
		Simulation sim;
		sim.init(vfSimulations.at(i)->getParameter(), vfSimulations.at(i)->getGrid(), vfSimulations.at(i)->getDataWriter());
		sim.run();
		testInfo->setSimulationEndTime(i);

		vfSimulations.at(i)->getCalculator()->calcAndCopyToTestResults();
	}

	testInfo->makeFinalTestOutput();

	testInfo->writeLogFile();
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