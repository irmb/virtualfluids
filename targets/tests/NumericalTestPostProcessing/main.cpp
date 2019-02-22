#include "Simulation/BasicSimulation.h"


#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/LogFileReader/LogFileReader.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantImp.h"
#include "Utilities/MathematicaFile/MathematicaFile.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactoryImp.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistantFactory/MathematicaAssistantFactoryImp.h"
#include "Utilities/MathematicaAssistant/MathematicaAssistant.h"

#include <memory>
#include <cfloat>
#include <cstdio>

#include <iostream>
#include <fstream>




int main(int argc, char **argv)
{
	std::vector<BasicSimulation> simulation;
	simulation.push_back(ShearWave);
	//simulation.push_back(TaylorGreenVortexUx);
	//simulation.push_back(TaylorGreenVortexUz);

	std::vector<Assistant> assistants;
	//assistants.push_back(Phi);
	//assistants.push_back(Ny);
	//assistants.push_back(L2Norm);
	//assistants.push_back(L2NormBetweenKernels);
	assistants.push_back(Time);

	std::vector<DataCombination> combination;
	combination.push_back(EqualSimulationsForDifferentKernels);
	//combination.push_back(EqualKernelSimulationsForDifferentViscosities);


	std::shared_ptr<LogFileReader> logFileReader = LogFileReader::getInstance();
	std::vector<std::shared_ptr<LogFileData> > logFileDataVector = logFileReader->readLogFilesInDirectoryToLogFileData("C:/Users/Timon/Documents/studienarbeitIRMB/logFiles");

	std::shared_ptr<MathematicaFile> aMathmaticaFile = MathematicaFile::getNewInstance("C:/Users/Timon/Desktop");

	std::shared_ptr<LogFileDataAssistant> assistentLogFile = LogFileDataAssistantImp::getNewInstance();

	std::shared_ptr<MathematicaFunctionFactory> functionFactory = MathematicaFunctionFactoryImp::getNewInstance();
	std::shared_ptr<MathematicaAssistantFactory> assistantFactory = MathematicaAssistantFactoryImp::getNewInstance();
	std::vector<std::shared_ptr<MathematicaAssistant> > mathematicaAssistants = assistantFactory->makeMathematicaAssistants(assistants, functionFactory);

	for (int sim = 0; sim < simulation.size(); sim++) {
		for (int comb = 0; comb < combination.size(); comb++) {
			std::vector<std::shared_ptr<LogFileDataGroup> > logFileDataSorted = assistentLogFile->findDataCombination(logFileDataVector, simulation.at(sim), combination.at(comb));
			for (int i = 0; i < logFileDataSorted.size(); i++) {
				for (int j = 0; j < mathematicaAssistants.size(); j++)
					mathematicaAssistants.at(j)->makeMathematicaOutput(logFileDataSorted.at(i), aMathmaticaFile);
			}
		}
	}
	
		
	aMathmaticaFile->finishFile();
	return 0;
}
