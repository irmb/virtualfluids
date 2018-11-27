#include "TestImp.h"

#include "Utilities\PostProcessingResults\PostProcessingResults.h"
#include "Utilities\TestSimulation\TestSimulation.h"

void TestImp::update()
{
	for (int i = 0; i < simulations.size(); i++)
	{
		if(simulationRun.at(i) == false)
		{
			if (simulations.at(i)->getSimulationRun())
			{
				simulationRun.at(i) = true;
				postProResults.at(i)->evaluate();
			}
		}
	}

	if (CheckAllSimulationRun())
		evaluate();				
}

void TestImp::addSimulation(std::shared_ptr<TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo, std::shared_ptr< PostProcessingResults> postProResult)
{
	simulations.push_back(sim);
	simInfos.push_back(simInfo);
	postProResults.push_back(postProResult);
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
