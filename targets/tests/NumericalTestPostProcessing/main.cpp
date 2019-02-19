#include "Simulation/BasicSimulation.h"

#include "Tests/PhiTest/MathematicaAssistant/PhiMathematicaAssistant.h"
#include "Tests/NyTest/MathematicaAssistant/NyMathematicaAssistant.h"

#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/LogFileReader/LogFileReader.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantImp.h"
#include "Utilities/MathematicaFile/MathematicaFile.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactoryImp.h"

#include <memory>
#include <cfloat>
#include <cstdio>

#include <iostream>
#include <fstream>




int main(int argc, char **argv)
{
	BasicSimulation simulation = ShearWave;
	//BasicSimulation simulation = TaylorGreenVortexUx;
	//BasicSimulation simulation = TaylorGreenVortexUz;

	std::shared_ptr<LogFileReader> logFileReader = LogFileReader::getInstance();
	//std::shared_ptr<LogFileData> logFileData = logFileReader->readLogFileToLogFileData("C:/Users/Timon/Documents/studienarbeitIRMB/logFiles/NumericalTestLogFiles/TaylorGreenVortexUx/viscosity_0.001/ux_ 0.016_Amplitude_ 0.005/CumulantAA2016CompSP27/logfile_20190211_153133_CumulantAA2016CompSP27_vis_0.001.txt");
	std::vector<std::shared_ptr<LogFileData> > logFileDataVector = logFileReader->readLogFilesInDirectoryToLogFileData("C:/Users/Timon/Desktop/logFiles");

	std::shared_ptr<LogFileDataAssistant> assistent = LogFileDataAssistantImp::getNewInstance();
	std::vector<std::shared_ptr<LogFileDataGroup> > logFileDataSorted = assistent->findEqualSimulationsForDifferentKernels(logFileDataVector, simulation);
	//std::vector<std::vector<std::shared_ptr<LogFileData> > > logFileDataSorted = assistent->findEqualKernelSimulationsForDifferentViscosities(logFileDataVector, simulation);

	std::shared_ptr<MathematicaFile> aMathmaticaFile = MathematicaFile::getNewInstance("C:/Users/Timon/Desktop");
	std::shared_ptr<MathematicaFunctionFactory> functionFactory = MathematicaFunctionFactoryImp::getNewInstance();

	std::shared_ptr<PhiMathematicaAssistant> mathematicaAssistantPhi = PhiMathematicaAssistant::getNewInstance(functionFactory);
	std::shared_ptr<NyMathematicaAssistant> mathematicaAssistantNy = NyMathematicaAssistant::getNewInstance(functionFactory);
	
	for (int i = 0; i < logFileDataSorted.size(); i++) {
		mathematicaAssistantPhi->makeMathematicaOutput(logFileDataSorted.at(i), aMathmaticaFile);
		mathematicaAssistantNy->makeMathematicaOutput(logFileDataSorted.at(i), aMathmaticaFile);
	}
		
	aMathmaticaFile->finishFile();
	return 0;
}
