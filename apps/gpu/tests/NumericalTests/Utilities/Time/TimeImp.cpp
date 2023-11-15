#include "TimeImp.h"

#include <sstream>

std::shared_ptr<TimeImp> TimeImp::getNewInstance()
{
    return std::shared_ptr<TimeImp>(new TimeImp());
}

TimeImp::TimeImp()
{

}

void TimeImp::setSimulationStartTime()
{
    simulationStartTime = time(NULL);
}

void TimeImp::setSimulationEndTime()
{
    simulationEndTime = time(NULL);
}

void TimeImp::setTestStartTime()
{
    testStartTime = clock();
}

void TimeImp::setTestEndTime()
{
    testEndTime = clock();
}

void TimeImp::setAnalyticalResultWriteStartTime()
{
    analyticalResultWriteStartTime = time(NULL);
}

void TimeImp::setAnalyticalResultWriteEndTime()
{
    analyticalResultWriteEndTime = time(NULL);
}

void TimeImp::setResultCheckStartTime()
{
    resultCheckStartTime = clock();
}

void TimeImp::setResultCheckEndTime()
{
    resultCheckEndTime = clock();
}

std::string TimeImp::getSimulationTime()
{
    std::ostringstream oss;
    oss << calcSimulationTime() << "sec";
    return oss.str();
}

std::string TimeImp::getResultCheckTime()
{
    std::ostringstream oss;
    oss << calcResultCheckTime() << "sec";
    return oss.str();
}

std::string TimeImp::getTestTime()
{
    std::ostringstream oss;
    oss << calcTestTime() << "sec";
    return oss.str();
}

std::string TimeImp::getAnalyticalResultWriteTime()
{
    std::ostringstream oss;
    oss << calcAnalyticalResultWriteTime() << "sec";
    return oss.str();
}

double TimeImp::calcSimulationTime()
{
    return difftime(simulationEndTime, simulationStartTime);
}

float TimeImp::calcResultCheckTime()
{
    float timeInMiliSec = ((float)(resultCheckEndTime - resultCheckStartTime) / CLOCKS_PER_SEC);
    return timeInMiliSec;
}

float TimeImp::calcTestTime()
{
    float timeInMiliSec = ((float)(testEndTime - testStartTime) / CLOCKS_PER_SEC);
    return timeInMiliSec;
}

double TimeImp::calcAnalyticalResultWriteTime()
{
    return difftime(analyticalResultWriteEndTime, analyticalResultWriteStartTime);
}
