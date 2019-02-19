#include "LogFileDataAssistantStrategyFactoryImp.h"

#include "Simulation/ShearWave/LogFileDataAssistantStrategy/ShearWaveLogFileDataAssistantStrategy.h"
#include "Simulation/TaylorGreenVortexUx/LogFileDataAssistantStrategy/TaylorGreenVortexUxLogFileDataAssistantStrategy.h"
#include "Simulation/TaylorGreenVortexUz/LogFileDataAssistantStrategy/TaylorGreenVortexUzLogFileDataAssistantStrategy.h"

std::shared_ptr<LogFileDataAssistantStrategyFactory> LogFileDataAssistantStrategyFactoryImp::getNewInstance()
{
	return std::shared_ptr<LogFileDataAssistantStrategyFactory>(new LogFileDataAssistantStrategyFactoryImp());
}

std::shared_ptr<LogFileDataAssistantStrategy> LogFileDataAssistantStrategyFactoryImp::makeLogFileDataAssistantStrategy(BasicSimulation sim)
{
	std::shared_ptr<LogFileDataAssistantStrategy> assistentStrategy;
	switch (sim)
	{
	case ShearWave:
		assistentStrategy = ShearWaveLogFileDataAssistantStrategy::getNewInstance();
		break;
	case TaylorGreenVortexUx:
		assistentStrategy = TaylorGreenVortexUxLogFileDataAssistantStrategy::getNewInstance();
		break;
	case TaylorGreenVortexUz:
		assistentStrategy = TaylorGreenVortexUzLogFileDataAssistantStrategy::getNewInstance();
		break;
	default:
		break;
	}
	return assistentStrategy;
}

LogFileDataAssistantStrategyFactoryImp::LogFileDataAssistantStrategyFactoryImp()
{
}
