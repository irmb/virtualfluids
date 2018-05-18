#include <gmock/gmock.h>
#include "mpi.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"

#include "Utilities/ConfigFileReader/ConfigFileReader.h"
#include "Utilities/TestCondition/TestCondition.h"
#include "Utilities/TestConditionFactory/TestConditionFactoryImp.h"
#include "Utilities/Calculator/Calculator.h"
#include "Utilities/TestResults/TestResults.h"
#include "Utilities/TestInformation/TestInformation.h"

static void startNumericalTests(const std::string &configFile)
{
	std::shared_ptr< ConfigFileReader> configReader = ConfigFileReader::getNewInstance();
	configReader->readConfigFile(configFile);

	std::vector< std::shared_ptr< TestParameter> > testPara = configReader->getTestParameter();
	std::shared_ptr< TestInformation> testInfo = configReader->getTestInformation();

	std::shared_ptr< TestConditionFactory> factory = TestConditionFactoryImp::getNewInstance();
	std::vector< std::shared_ptr< TestCondition> > testConditions = factory->makeTestConditions(testPara);

	for (int i = 0; i < testConditions.size(); i++)
	{
		testInfo->makeSimulationHeadOutput(i);
		testInfo->setSimulationStartTime(i);
		Simulation sim;
		sim.init(testConditions.at(i)->getParameter(), testConditions.at(i)->getGrid(), testConditions.at(i)->getDataWriter());
		sim.run();
		testInfo->setSimulationEndTime(i);

		testConditions.at(i)->getCalculator()->calcAndCopyToTestResults();
		testConditions.at(i)->getTestResults()->evaluate();
	}

	testInfo->makeFinalTestOutput();
	testInfo->writeLogFile();
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	if (argc > 1)
		startNumericalTests(argv[1]);

    MPI_Finalize();

	return 0;
}