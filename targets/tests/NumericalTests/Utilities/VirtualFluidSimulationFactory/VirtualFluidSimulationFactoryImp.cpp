#include "VirtualFluidSimulationFactoryImp.h"

#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\VirtualFluidSimulation\VirtualFluidSimulationImp.h"

std::shared_ptr<VirtualFluidSimulationFactory> VirtualFluidSimulationFactoryImp::getNewInstance()
{
	return std::shared_ptr< VirtualFluidSimulationFactory>(new VirtualFluidSimulationFactoryImp());
}

VirtualFluidSimulationFactoryImp::VirtualFluidSimulationFactoryImp()
{

}

std::vector<std::shared_ptr<VirtualFluidSimulation>> VirtualFluidSimulationFactoryImp::makeVirtualFluidSimulations(std::vector< std::shared_ptr< SimulationParameter> > simPara)
{
	std::vector< std::shared_ptr<VirtualFluidSimulation> > vfSimulations;

	for (int i = 0; i < simPara.size(); i++) {
		std::shared_ptr< VirtualFluidSimulationImp> vfSim = VirtualFluidSimulationImp::getNewInstance();

		vfSim->initParameter(simPara.at(i)->getViscosity(), simPara.at(i)->getGridPath(), simPara.at(i)->getFilePath(), simPara.at(i)->getNumberOfGridLevels(), simPara.at(i)->getEndTime(), simPara.at(i)->getTimeStepLength(), simPara.at(i)->getDevices(), simPara.at(i)->getMaxVelocity());
		vfSim->initInitialConditions(simPara.at(i)->getInitialCondition());
		vfSim->initGridProvider();
		vfSim->initCalculator(simPara.at(i)->getCalculator());
		vfSim->initSimulationResults(simPara.at(i)->getLx(), simPara.at(i)->getLz(), simPara.at(i)->getTimeStepLength());
		vfSim->setTestResults(simPara.at(i)->getTestResults());
		vfSim->initDataWriter(simPara.at(i)->getYSliceForCalculation(), simPara.at(i)->getStartTimeCalculation(), simPara.at(i)->getEndTime(), simPara.at(i)->getTimeStepLength(), simPara.at(i)->getWriteFiles(), simPara.at(i)->getStartTimeDataWriter());
		vfSimulations.push_back(vfSim);
	}

	return vfSimulations;
}
