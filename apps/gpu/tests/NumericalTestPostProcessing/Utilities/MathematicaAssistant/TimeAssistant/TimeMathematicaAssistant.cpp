#include "TimeMathematicaAssistant.h"

#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"
#include "Utilities/LogFileData/LogFileData.h"

#include "Utilities/DataPoint/DataPoint.h"
#include "Utilities/MathematicaFunctionFactory/MathematicaFunctionFactory.h"

std::shared_ptr<TimeMathematicaAssistant> TimeMathematicaAssistant::getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
    return std::shared_ptr<TimeMathematicaAssistant>(new TimeMathematicaAssistant(functionFactory));
}

TimeMathematicaAssistant::TimeMathematicaAssistant()
{

}

TimeMathematicaAssistant::TimeMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory) : MathematicaAssistantImp(functionFactory)
{

}

void TimeMathematicaAssistant::makeSimulationTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    std::vector<std::vector<double> > simTimes(logFileData->getGroupSize());
    for (int i = 0; i < logFileData->getGroupSize(); i++)
        for(int j = 0; j < logFileData->getLogFileData(i)->getSimTime().size(); j++)
            simTimes.at(i).push_back((double)logFileData->getLogFileData(i)->getSimTime().at(j));

    makeTimeMathematicaOutput(logFileData, aMathmaticaFile, simTimes, "SimulationTime");
}

void TimeMathematicaAssistant::makeTestTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    std::vector<std::vector<double> > testTimes;
    for (int i = 0; i < logFileData->getGroupSize(); i++)
        testTimes.push_back(logFileData->getLogFileData(i)->getTestTime());

    makeTimeMathematicaOutput(logFileData, aMathmaticaFile, testTimes, "TestTime");
}

void TimeMathematicaAssistant::makeAnalyticalWriteTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    std::vector<std::vector<double> > analyticalWriteTime(logFileData->getGroupSize());
    for (int i = 0; i < logFileData->getGroupSize(); i++)
        for (int j = 0; j < logFileData->getLogFileData(i)->getAnalyticalVTKWritingTime().size(); j++)
            analyticalWriteTime.at(i).push_back((double)logFileData->getLogFileData(i)->getAnalyticalVTKWritingTime().at(j));

    makeTimeMathematicaOutput(logFileData, aMathmaticaFile, analyticalWriteTime, "AnalyticalVTKWritingTime");
}

void TimeMathematicaAssistant::makeTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::vector<double>> times, std::string timeName)
{
    std::vector<std::vector<double> > grids;
    std::vector<std::string> listNames;
    for (int i = 0; i < logFileData->getGroupSize(); i++) {
        grids.push_back(logFileData->getLogFileData(i)->getBasicGridLengths());
        listNames.push_back(logFileData->getLogFileData(i)->getSimulationSigniture() + timeName);
    }

    std::vector<std::vector<std::shared_ptr<DataPoint> > > dataPointGroup;

    for (int j = 0; j < grids.size(); j++) {
        std::vector<std::shared_ptr<DataPoint> > dataPoints;
        for (int k = 0; k < times.at(j).size(); k++)
            dataPoints.push_back(DataPoint::getNewInstance(grids.at(j).at(k), times.at(j).at(k)));
        dataPointGroup.push_back(dataPoints);
    }
    std::vector<std::shared_ptr<MathematicaPointList> > pointList;
    for (int j = 0; j < dataPointGroup.size(); j++) {
        std::shared_ptr<MathematicaPointList> aPointList = functionFactory->makeMathematicaPointList(aMathmaticaFile, listNames.at(j), dataPointGroup.at(j));
        pointList.push_back(aPointList);
    }
    std::shared_ptr<MathematicaListPlot> listLogLogPlot = functionFactory->makeMathematicaListPlot(aMathmaticaFile, pointList, "ListLogLogPlot", "L[dx]", "Time [sec]");
}

void TimeMathematicaAssistant::makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    makeSimulationTimeMathematicaOutput(logFileData, aMathmaticaFile);
    makeTestTimeMathematicaOutput(logFileData, aMathmaticaFile);
    makeAnalyticalWriteTimeMathematicaOutput(logFileData, aMathmaticaFile);
}