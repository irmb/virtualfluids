#include <mpi.h>
#include <gmock/gmock.h>

#include "Utilities/ConfigFileReaderNT/ConfigFileReaderNT.h"
#include "Utilities/LogFileQueue/LogFileQueue.h"
#include "Utilities/NumericalTestFactory/NumericalTestFactoryImp.h"
#include "Utilities/TestQueue/TestQueue.h"
#include "Utilities/VirtualFluidSimulation/VirtualFluidSimulation.h"
#include "Utilities/VirtualFluidSimulationFactory/VirtualFluidSimulationFactoryImp.h"

//validation
#include <fstream>
#include <iostream>
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"

static void validateTestSuite()
{
	const int timeSteps = 10;
	const int begin = 11;
	const int end = 20;
	const int l0 = 32;
	const double viscosity = 0.001;
	std::string kernelName = "Cum";

	std::vector<int> xLength{32,64,128,256,512};
	std::vector<int> zLength(xLength.size());
	std::vector<int> timeStepLength(xLength.size());
	for (int i = 0; i < xLength.size(); i++) {
		zLength.at(i) = xLength.at(i) * 3 / 2;
		timeStepLength.at(i) = (int)1000 * xLength.at(i)*xLength.at(i) / l0 / l0;
	}
		
	std::vector<std::vector<std::ostringstream>> filePaths;
	filePaths.resize(xLength.size());
	
	for (int j = 0; j < xLength.size(); j++) {
		filePaths.at(j).resize(timeSteps);
		for (int i = begin; i <= end; i++)
			filePaths.at(j).at(i - begin) << "C:/Users/Timon/Desktop/Auswertung_TGV_10hm3_X_CumD3Q27F3/Auswertung_TGV_10hm3_X_CumD3Q27F3/" << kernelName << "_" << xLength.at(j) << "_3_" << zLength.at(j) << "_AD_X_" << i*timeStepLength.at(j) << ".dat";

	}
	std::vector<std::vector<std::vector<double>>> dataForOneSimulationGroup;
	
	for (int j = 0; j < filePaths.size(); j++) {
		std::vector<std::vector<double>> dataForOneSimulation;
		dataForOneSimulation.resize(timeSteps);
		for (int i = 0; i < filePaths.at(j).size(); i++) {
			std::ifstream file;
			file.open(filePaths.at(j).at(i).str());

			if (file.is_open()) {
				double data = 0.0;
				while (file >> data)
					dataForOneSimulation.at(i).push_back(data);

				file.close();
			}
			else
				int stop = 1;
		}
		dataForOneSimulationGroup.push_back(dataForOneSimulation);
	}

	std::shared_ptr<FFTCalculator> calulator = FFTCalculator::getInstance();

	std::vector<double> phiDifForOneSimGroup;
	std::vector<double> nyDifForOneSimGroup;
	for (int i = 0; i < dataForOneSimulationGroup.size(); i++) {
		int timeStepLength = 1000 * xLength.at(i)*xLength.at(i) / l0 / l0;

		double phiDiff = calulator->calcPhiDiff(dataForOneSimulationGroup.at(i), false, xLength.at(i), zLength.at(i), timeStepLength);
		double ny = calulator->calcNy(dataForOneSimulationGroup.at(i), false, xLength.at(i), zLength.at(i), timeStepLength);
		double nyDiff = abs(ny - viscosity) / viscosity;
		phiDifForOneSimGroup.push_back(phiDiff);
		nyDifForOneSimGroup.push_back(nyDiff);
	}



	std::fstream dataOutPhi;
	std::string dataOutFilePathPhi = "C:/Users/Timon/Desktop/Auswertung_TGV_10hm3_X_CumD3Q27F3/NumericalTestAuswertung/" + kernelName + "_PhiDiff.dat";
	dataOutPhi.open(dataOutFilePathPhi, std::ios::out);

	std::fstream dataOutNy;
	std::string dataOutFilePathNy = "C:/Users/Timon/Desktop/Auswertung_TGV_10hm3_X_CumD3Q27F3/NumericalTestAuswertung/" + kernelName + "_NyDiff.dat";
	dataOutNy.open(dataOutFilePathNy, std::ios::out);

	for (int i = 0; i < phiDifForOneSimGroup.size(); i++) {
		dataOutPhi << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		dataOutPhi << phiDifForOneSimGroup.at(i);
		dataOutNy << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		dataOutNy << nyDifForOneSimGroup.at(i);

		if (i < phiDifForOneSimGroup.size() - 1) {
			dataOutPhi << std::endl;
			dataOutNy << std::endl;
		}
	}

	dataOutPhi.close();
}


static bool startNumericalTests(const std::string &configFile)
{
	std::shared_ptr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance(configFile);
	configReader->readConfigFile();

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

	return testQueue->getNumberOfFailedTests() > 0;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	//validateTestSuite();

	bool tests_passed{false};

	if (argc > 1)
        tests_passed = startNumericalTests(argv[1]);
	else
		std::cout << "Configuration file must be set!: lbmgm <config file>" << std::endl << std::flush;

    MPI_Finalize();

	return tests_passed;
}
