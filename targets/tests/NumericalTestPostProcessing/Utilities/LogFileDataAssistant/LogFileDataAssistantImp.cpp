#include "LogFileDataAssistantImp.h"

#include "Utilities/AlmostEquals.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroupImp.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategy.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategyFactory/LogFileDataAssistantStrategyFactoryImp.h"

std::vector<std::vector<std::shared_ptr<LogFileData>>> LogFileDataAssistantImp::sortLogFileDataAfterKernels(std::vector<std::shared_ptr<LogFileData>> logFileData)
{
	std::vector<std::string> kernelNames;
	for (int i = 0; i < logFileData.size(); i++) {
		if (i == 0)
			kernelNames.push_back(logFileData.at(i)->getKernel());
		else
		{
			bool newKernel = true;
			for (int j = 0; j < kernelNames.size(); j++) {
				if (kernelNames.at(i) == logFileData.at(i)->getKernel())
					newKernel = false;
			}
			if (newKernel)
				kernelNames.push_back(logFileData.at(i)->getKernel());
		}
	}

	std::vector<std::vector<std::shared_ptr<LogFileData> > > logFileDataAfterKernels;
	logFileDataAfterKernels.resize(kernelNames.size());
	for (int i = 0; i < kernelNames.size(); i++) {
		for (int j = 0; j < logFileData.size(); j++) {
			if (kernelNames.at(i) == logFileData.at(j)->getKernel())
				logFileDataAfterKernels.at(i).push_back(logFileData.at(j));
		}
	}
	return logFileDataAfterKernels;
}

bool LogFileDataAssistantImp::checkEqualSimulationData(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
	if (logFileData1->getNumberOfTimeSteps() != logFileData2->getNumberOfTimeSteps())
		return false;
	if (logFileData1->getBasisTimeStepLength() != logFileData2->getBasisTimeStepLength())
		return false;

	return true;
}

bool LogFileDataAssistantImp::checkEqualViscosity(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
	if (!equalDouble(logFileData1->getViscosity(), logFileData2->getViscosity()))
		return false;

	return true;
}

bool LogFileDataAssistantImp::checkEqualKernel(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
	if (logFileData1->getKernel() != logFileData2->getKernel())
		return false;

	return true;
}

std::vector<std::shared_ptr<LogFileDataGroup>> LogFileDataAssistantImp::castToLogFileDataGroup(std::vector<std::shared_ptr<LogFileDataGroupImp>> data)
{
	std::vector<std::shared_ptr<LogFileDataGroup>> casted;

	for (int i = 0; i < data.size(); i++)
		casted.push_back(data.at(i));
	return casted;
}

bool LogFileDataAssistantImp::equalDouble(double num1, double num2)
{
	const FloatingPoint<double> lhs(num1), rhs(num2);

	if (lhs.AlmostEquals(rhs))
		return true;
	return false;
}

std::vector<std::shared_ptr<LogFileData>> LogFileDataAssistantImp::getSimulationGroupLogFileData(std::string simName, std::vector<std::shared_ptr<LogFileData>> allLogFileData)
{
	std::vector<std::shared_ptr<LogFileData>> simGroupLogFileData;

	for (int i = 0; i < allLogFileData.size(); i++) {
		if (allLogFileData.at(i)->getSimName() == simName)
			simGroupLogFileData.push_back(allLogFileData.at(i));
	}

	return simGroupLogFileData;
}

std::shared_ptr<LogFileDataAssistant> LogFileDataAssistantImp::getNewInstance()
{
	return std::shared_ptr<LogFileDataAssistant>(new LogFileDataAssistantImp());
}

std::vector<std::shared_ptr<LogFileDataGroup>> LogFileDataAssistantImp::findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData>> allLogFileData, BasicSimulation simulation)
{
	std::shared_ptr<LogFileDataAssistantStrategy> strategy = assistentStrategyFactory->makeLogFileDataAssistantStrategy(simulation);

	std::vector<std::shared_ptr<LogFileData>> myLogFileData = getSimulationGroupLogFileData(strategy->getSimulationName(), allLogFileData);

	std::vector<std::shared_ptr<LogFileDataGroupImp>  > kernelGroups;
	kernelGroups.push_back(LogFileDataGroupImp::getNewInstance());
	kernelGroups.at(0)->addLogFileData(myLogFileData.at(0));

	for (int i = 1; i < myLogFileData.size(); i++) {
		bool added = false;
		for (int j = 0; j < kernelGroups.size(); j++) {
			if (checkEqualSimulationData(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
				if (checkEqualViscosity(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
					if (strategy->checkSimulationParamerer(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
						kernelGroups.at(j)->addLogFileData(myLogFileData.at(i));
						added = true;
					}
				}
			}
		}
		if (!added) {
			std::shared_ptr<LogFileDataGroupImp> newGroup = LogFileDataGroupImp::getNewInstance();
			newGroup->addLogFileData(myLogFileData.at(i));
			kernelGroups.push_back(newGroup);
		}
	}

	return castToLogFileDataGroup(kernelGroups);
}

std::vector<std::shared_ptr<LogFileDataGroup> > LogFileDataAssistantImp::findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData>> allLogFileData, BasicSimulation simulation)
{
	
	std::shared_ptr<LogFileDataAssistantStrategy> strategy = assistentStrategyFactory->makeLogFileDataAssistantStrategy(simulation);

	std::vector<std::shared_ptr<LogFileData>> myLogFileData = getSimulationGroupLogFileData(strategy->getSimulationName(), allLogFileData);

	std::vector<std::shared_ptr<LogFileDataGroupImp>  > viscosityGroups;
	viscosityGroups.push_back(LogFileDataGroupImp::getNewInstance());
	viscosityGroups.at(0)->addLogFileData(myLogFileData.at(0));

	for (int i = 1; i < myLogFileData.size(); i++) {
		bool added = false;
		for (int j = 0; j < viscosityGroups.size(); j++) {
			if (checkEqualSimulationData(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
				if (checkEqualKernel(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
					if (strategy->checkSimulationParamerer(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
						viscosityGroups.at(j)->addLogFileData(myLogFileData.at(i));
						added = true;
					}
				}
			}
		}
		if (!added) {
			std::shared_ptr<LogFileDataGroupImp> newGroup = LogFileDataGroupImp::getNewInstance();
			newGroup->addLogFileData(myLogFileData.at(i));
			viscosityGroups.push_back(newGroup);
		}
	}

	return castToLogFileDataGroup(viscosityGroups);
}

LogFileDataAssistantImp::LogFileDataAssistantImp()
{
	assistentStrategyFactory = LogFileDataAssistantStrategyFactoryImp::getNewInstance();
}