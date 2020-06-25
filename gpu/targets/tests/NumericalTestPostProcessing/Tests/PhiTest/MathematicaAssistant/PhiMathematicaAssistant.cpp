#include "PhiMathematicaAssistant.h"

#include "Tests/PhiTest/LogFileData/PhiLogFileData.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"


std::shared_ptr<PhiMathematicaAssistant> PhiMathematicaAssistant::getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
	return std::shared_ptr<PhiMathematicaAssistant>(new PhiMathematicaAssistant(functionFactory));
}

void PhiMathematicaAssistant::makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
	std::shared_ptr<SortedDataPhi> mySortedData = sortLogFileData(logFileData);

	makePhiDiffMathematicaOutput(aMathmaticaFile, mySortedData);
	makeOrderOfAccuracyMathematicaOutput(aMathmaticaFile, mySortedData);

}

PhiMathematicaAssistant::PhiMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory) : MathematicaAssistantImp(functionFactory)
{
}

bool PhiMathematicaAssistant::checkTestParameter(std::shared_ptr<PhiLogFileData> logFileData1, std::shared_ptr<PhiLogFileData> logFileData2)
{
	if (logFileData1->getStartTimeStepCalculation() != logFileData2->getStartTimeStepCalculation())
		return false;
	if (logFileData1->getEndTimeStepCalculation() != logFileData2->getEndTimeStepCalculation())
		return false;
	if (logFileData1->getDataToCalc() != logFileData2->getDataToCalc())
		return false;

	return true;
}

std::shared_ptr<SortedDataPhi> PhiMathematicaAssistant::sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData)
{
	std::vector<std::vector<std::shared_ptr<PhiLogFileData> > > testLogFileData;
	std::vector<std::vector<std::string> > basicListNames;
	for (int i = 0; i < logFileData->getLogFileData(0)->getPhiLogFileData().size(); i++) {
		std::vector<std::shared_ptr<PhiLogFileData> > aTestLogFileDataGroup;
		aTestLogFileDataGroup.push_back(logFileData->getLogFileData(0)->getPhiLogFileData().at(i));
		std::vector<std::string> aListNameGroup;
		aListNameGroup.push_back(logFileData->getLogFileData(0)->getSimulationSigniture());
		basicListNames.push_back(aListNameGroup);
	}
	for (int i = 0; i < logFileData->getGroupSize(); i++) {
		for (int j = 0; j < logFileData->getLogFileData(i)->getPhiLogFileData().size(); j++) {
			std::string dataToCalc = logFileData->getLogFileData(i)->getPhiLogFileData().at(j)->getDataToCalc();
			bool added = false;
			for (int k = 0; k < testLogFileData.size(); k++) {
				if (checkTestParameter(logFileData->getLogFileData(i)->getPhiLogFileData().at(j), testLogFileData.at(k).at(0))) {
					testLogFileData.at(k).push_back(logFileData->getLogFileData(i)->getPhiLogFileData().at(j));
					basicListNames.at(k).push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
					added = true;
				}
			}
			if (!added) {
				std::vector<std::shared_ptr<PhiLogFileData> > aTestLogFileDataGroup;
				aTestLogFileDataGroup.push_back(logFileData->getLogFileData(i)->getPhiLogFileData().at(j));
				testLogFileData.push_back(aTestLogFileDataGroup);
				std::vector<std::string> aListNameGroup;
				aListNameGroup.push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
				basicListNames.push_back(aListNameGroup);
			}
		}
	}
	std::shared_ptr<SortedDataPhi> mySortedData = std::shared_ptr<SortedDataPhi>(new SortedDataPhi);
	mySortedData->basicListNames = basicListNames;
	mySortedData->testLogFileData = testLogFileData;

	return mySortedData;
}

void PhiMathematicaAssistant::makePhiDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData)
{
	for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
		std::vector<std::vector<double> > gridLengths;
		std::vector<std::vector<double> > phiDiff;
		std::vector<std::string> aBasicListNamesList;
		for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
			gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
			phiDiff.push_back(sortedData->testLogFileData.at(i).at(j)->getPhiDiff());
			aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
		}
		std::vector<std::string> finalListNames = finalizeListNames(aBasicListNamesList, "PhiDiff", sortedData->testLogFileData.at(i).at(0)->getDataToCalc());
		addSecondOrderOfAccuracyRef(gridLengths, phiDiff, finalListNames);
		addFourthOrderOfAccuracyRef(gridLengths, phiDiff, finalListNames);
		addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNames, gridLengths, phiDiff, "L[dx]", "Err Phi[-]");
	}
}

void PhiMathematicaAssistant::makeOrderOfAccuracyMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData)
{
	for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
		for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
			std::vector<std::vector<double>> ooA = sortedData->testLogFileData.at(i).at(j)->getOrderOfAccuracy();
			std::string basicListName = sortedData->basicListNames.at(i).at(j);
			std::string dataToCalc = sortedData->testLogFileData.at(i).at(j)->getDataToCalc();
			std::string finalListName = finalizeListName(basicListName, "PhiDiffOrderOfAccuracy", dataToCalc);

			addListOfListsToMathematicaFile(aMathmaticaFile, finalListName, ooA);
		}
	}
}

PhiMathematicaAssistant::PhiMathematicaAssistant()
{
}
