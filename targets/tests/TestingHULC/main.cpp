#include <gmock/gmock.h>
#include "mpi.h"

#include "VirtualFluids_GPU\LBM\Simulation.h"

#include "Utilities\Reader\Reader.h"
#include "Utilities\EvaluationParameter\EvaluationParameter.h"
#include "Utilities\TestCondition\TestCondition.h"
#include "Utilities\Calculator\Calculator.h"
#include "Utilities\LogFileWriter\LogFileWriter.h"
#include "Utilities\TestConditionFactory\TestConditionFactoryImp.h"

#include "Tests\DataCollector\DataCollector.h"
#include "Tests\DataQueue\DataQueue.h"

#include "Tests\OrderOfAccuracy\OrderOfAccuracy.h"

#include "Tests\TestCout\TestCout.h"

#include <iostream>


//muss nicht unbedingt
#include "Utilities\Results\Results.h"

using std::shared_ptr;


const int numberOfTests = 5;

DataQueue nuTGV[numberOfTests - 1];
DataQueue phiTGV[numberOfTests - 1];
DataQueue nuSW[numberOfTests - 1];
DataQueue phiSW[numberOfTests - 1];

static void testHULC(const std::string &configFile)
{
	std::shared_ptr< Reader > configReader = Reader::getNewInstance(configFile);

	std::vector< std::shared_ptr< EvaluationParameter > > evaPara = configReader->makeEvaluationParameter();
	std::shared_ptr<TestInformation> testInfo = configReader->makeTestInformation();
	std::vector<std::shared_ptr<TestParameter> > testPara = configReader->makeTestParameter();

	std::shared_ptr<TestConditionFactory> factory = TestConditionFactoryImp::getNewInstance(testPara);
	std::vector<std::shared_ptr<TestCondition>> testConditions = factory->makeTestConditions();

	DataCollector tgvCollector = DataCollector(nuTGV, phiTGV, numberOfTests - 1, "TaylorGreenVortex");
	DataCollector swCollector = DataCollector(nuSW, phiSW, numberOfTests - 1, "ShearWave");

	for (int i = 0; i < testConditions.size(); i++)
	{
		evaPara.at(i)->setStartTime();
		TEST_HEAD(evaPara.at(i)->getTestName(), evaPara.at(i)->getLx());
		Simulation sim;
		sim.init(testConditions.at(i)->getParameter(), testConditions.at(i)->getGrid(), testConditions.at(i)->getDataWriter());
		sim.run();
		evaPara.at(i)->setEndTime();

		std::shared_ptr<Calulator> calc = std::shared_ptr<Calulator>(new Calulator(testConditions.at(i)->getSimulationResults(), evaPara.at(i)));

		double nu = calc->calcNu();
		double nudiff = calc->calcNuDiff(nu);
		double phidiff = calc->calcPhiDiff();

		tgvCollector.addNuDiffAndPhi(nudiff, phidiff, evaPara.at(i));
		swCollector.addNuDiffAndPhi(nudiff, phidiff, evaPara.at(i));
	}

	std::shared_ptr<LogFileWriter> logFile = std::shared_ptr<LogFileWriter>(new LogFileWriter(evaPara, testInfo));
	logFile->makeDataQueueOutput(nuTGV, numberOfTests - 1);
	logFile->makeDataQueueOutput(phiTGV, numberOfTests - 1);
	logFile->makeDataQueueOutput(nuSW, numberOfTests - 1);
	logFile->makeDataQueueOutput(phiSW, numberOfTests - 1);
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	if (argc > 1)
		testHULC(argv[1]);

	::testing::InitGoogleTest(&argc, argv);
    MPI_Finalize();

	return RUN_ALL_TESTS();
}

INSTANTIATE_TEST_CASE_P(TaylorGreenVortexNu, OrderOfAccuracy, ValuesIn(nuTGV));
INSTANTIATE_TEST_CASE_P(TaylorGreenVortexPhi, OrderOfAccuracy, ValuesIn(phiTGV));

INSTANTIATE_TEST_CASE_P(ShearWaveNu, OrderOfAccuracy, ValuesIn(nuSW));
INSTANTIATE_TEST_CASE_P(ShearWavePhi, OrderOfAccuracy, ValuesIn(phiSW));