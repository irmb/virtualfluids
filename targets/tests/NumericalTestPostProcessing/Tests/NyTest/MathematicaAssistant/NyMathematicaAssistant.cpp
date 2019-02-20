#include "NyMathematicaAssistant.h"

#include "Tests/NyTest/LogFileData/NyLogFileData.h"

#include "Utilities/DataPoint/DataPoint.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactory.h"

std::shared_ptr<NyMathematicaAssistant> NyMathematicaAssistant::getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
	return std::shared_ptr<NyMathematicaAssistant>(new NyMathematicaAssistant(functionFactory));
}

void NyMathematicaAssistant::makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
	std::shared_ptr<SortedDataNy> mySortedData = sortLogFileData(logFileData);

	makeNyDiffMathematicaOutput(aMathmaticaFile, mySortedData);
}

NyMathematicaAssistant::NyMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory) : functionFactory(functionFactory)
{
}

bool NyMathematicaAssistant::checkTestParameter(std::shared_ptr<NyLogFileData> logFileData1, std::shared_ptr<NyLogFileData> logFileData2)
{
	if (logFileData1->getStartTimeStepCalculation() != logFileData2->getStartTimeStepCalculation())
		return false;
	if (logFileData1->getEndTimeStepCalculation() != logFileData2->getEndTimeStepCalculation())
		return false;
	if (logFileData1->getDataToCalc() != logFileData2->getDataToCalc())
		return false;

	return true;
}

std::shared_ptr<SortedDataNy> NyMathematicaAssistant::sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData)
{
	std::vector<std::vector<std::shared_ptr<NyLogFileData> > > testLogFileData;
	std::vector<std::vector<std::string> > basicListNames;
	for (int i = 0; i < logFileData->getLogFileData(0)->getNyLogFileData().size(); i++) {
		std::vector<std::shared_ptr<NyLogFileData> > aTestLogFileDataGroup;
		aTestLogFileDataGroup.push_back(logFileData->getLogFileData(0)->getNyLogFileData().at(i));
		std::vector<std::string> aListNameGroup;
		aListNameGroup.push_back(logFileData->getLogFileData(0)->getSimulationSigniture());
		basicListNames.push_back(aListNameGroup);
	}
	for (int i = 0; i < logFileData->getGroupSize(); i++) {
		for (int j = 0; j < logFileData->getLogFileData(i)->getNyLogFileData().size(); j++) {
			std::string dataToCalc = logFileData->getLogFileData(i)->getNyLogFileData().at(j)->getDataToCalc();
			bool added = false;
			for (int k = 0; k < testLogFileData.size(); k++) {
				if (checkTestParameter(logFileData->getLogFileData(i)->getNyLogFileData().at(j), testLogFileData.at(k).at(0))) {
					testLogFileData.at(k).push_back(logFileData->getLogFileData(i)->getNyLogFileData().at(j));
					basicListNames.at(k).push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
					added = true;
				}
			}
			if (!added) {
				std::vector<std::shared_ptr<NyLogFileData> > aTestLogFileDataGroup;
				aTestLogFileDataGroup.push_back(logFileData->getLogFileData(i)->getNyLogFileData().at(j));
				testLogFileData.push_back(aTestLogFileDataGroup);
				std::vector<std::string> aListNameGroup;
				aListNameGroup.push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
				basicListNames.push_back(aListNameGroup);
			}
		}
	}
	std::shared_ptr<SortedDataNy> mySortedData = std::shared_ptr<SortedDataNy>(new SortedDataNy);
	mySortedData->basicListNames = basicListNames;
	mySortedData->testLogFileData = testLogFileData;

	return mySortedData;
}

void NyMathematicaAssistant::addListLogLogPlotToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::string> listNames, std::vector<std::vector<double> > xAxesData, std::vector<std::vector<double> > yAxesData, std::string labelXAxes, std::string labelYAxes)
{
	std::vector<std::vector<std::shared_ptr<DataPoint> > > dataPointGroup;

	for (int i = 0; i < xAxesData.size(); i++) {
		std::vector<std::shared_ptr<DataPoint> > dataPoints;
		for (int j = 0; j < xAxesData.at(i).size(); j++)
			dataPoints.push_back(DataPoint::getNewInstance(xAxesData.at(i).at(j), yAxesData.at(i).at(j)));
		dataPointGroup.push_back(dataPoints);
	}
	std::vector<std::shared_ptr<MathematicaPointList> > pointList;
	for (int i = 0; i < dataPointGroup.size(); i++) {
		std::shared_ptr<MathematicaPointList> aPointList = functionFactory->makeMathematicaPointList(aMathmaticaFile, listNames.at(i), dataPointGroup.at(i));
		pointList.push_back(aPointList);
	}
	std::shared_ptr<MathematicaListPlot> listLogLogPlot = functionFactory->makeMathematicaListPlot(aMathmaticaFile, pointList, "ListLogLogPlot", labelXAxes, labelYAxes);
}

std::vector<std::string> NyMathematicaAssistant::finalizeListNames(std::vector<std::string> basicListNames, std::string dataToCalc, std::string testName)
{
	std::vector<std::string> finalListNames;

	for (int i = 0; i < basicListNames.size(); i++)
		finalListNames.push_back(testName + basicListNames.at(i) + dataToCalc);

	return finalListNames;
}

void NyMathematicaAssistant::makeNyDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataNy> sortedData)
{
	for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
		std::vector<std::vector<double> > gridLengths;
		std::vector<std::vector<double> > nyDiff;
		std::vector<std::string> aBasicListNamesList;
		for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
			gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
			nyDiff.push_back(sortedData->testLogFileData.at(i).at(j)->getNyDiff());
			aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
		}
		std::vector<std::string> finalListNames = finalizeListNames(aBasicListNamesList, sortedData->testLogFileData.at(i).at(0)->getDataToCalc(), "NyDiff");
		addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNames, gridLengths, nyDiff, "L[dx]", "Err Ny");
	}
}

NyMathematicaAssistant::NyMathematicaAssistant()
{
}
