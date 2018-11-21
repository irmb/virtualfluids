#ifndef TEST_IMP_H
#define TEST_IMP_H

#include "Utilities\Test\Test.h"

#include <vector>

class TestSimulation;
class SimulationResults;
class SimulationInfo;

class TestImp : public Test
{
public:
	void update();
	void addSimulation(std::shared_ptr< TestSimulation> sim, std::shared_ptr< SimulationInfo> simInfo);

	virtual void evaluate() = 0;
	
protected:
	TestImp();
	bool CheckAllSimulationRun();


	std::vector< bool> simulationRun;
	std::vector< std::shared_ptr< TestSimulation>> simulations;
	std::vector< std::shared_ptr< SimulationResults>> simResults;
	std::vector< std::shared_ptr< SimulationInfo>> simInfos;
};
#endif 