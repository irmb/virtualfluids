#include "TestImp.h"

#include "Utilities\PostProcessingStrategy\PostProcessingStrategy.h"
#include "Utilities\NumericalTestSimulation\NumericalTestSimulation.h"

void TestImp::update()
{
	for (int i = 0; i < simulations.size(); i++)
	{
		if(simulationRun.at(i) == false)
		{
			if (simulations.at(i)->getSimulationRun())
			{
				simulationRun.at(i) = true;
				postProStrategies.at(i)->evaluate();
			}
		}
	}

	if (CheckAllSimulationRun())
		evaluate();				
}

void TestImp::addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PostProcessingStrategy> postProStrategy)
{
	simulations.push_back(sim);
	simInfos.push_back(simInfo);
	postProStrategies.push_back(postProStrategy);
	simulationRun.push_back(false);
}

std::string TestImp::getSimulationName()
{
	return simulationName;
}

TestImp::TestImp(std::shared_ptr<ColorConsoleOutput> colorOutput) : colorOutput(colorOutput)
{
	simulationRun.resize(0);
	simulations.resize(0);
	simInfos.resize(0);
}

bool TestImp::CheckAllSimulationRun()
{
	for(int i=0; i< simulationRun.size(); i++)
		if(simulationRun.at(i)==false)
			return false;
	
	return true;
}