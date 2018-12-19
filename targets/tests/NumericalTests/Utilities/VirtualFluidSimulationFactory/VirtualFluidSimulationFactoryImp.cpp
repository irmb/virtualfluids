#include "VirtualFluidSimulationFactoryImp.h"

#include "Utilities\TestSimulation\TestSimulation.h"
#include "Utilities\SimulationParameter\SimulationParameter.h"
#include "Utilities\VirtualFluidSimulation\VirtualFluidSimulationImp.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"

std::shared_ptr<VirtualFluidSimulationFactory> VirtualFluidSimulationFactoryImp::getNewInstance()
{
	return std::shared_ptr< VirtualFluidSimulationFactory>(new VirtualFluidSimulationFactoryImp());
}

VirtualFluidSimulationFactoryImp::VirtualFluidSimulationFactoryImp()
{

}

std::vector<std::shared_ptr<VirtualFluidSimulation>> VirtualFluidSimulationFactoryImp::makeVirtualFluidSimulations(std::vector< std::shared_ptr< TestSimulation> > testSim)
{
	std::vector< std::shared_ptr< VirtualFluidSimulation>> vfSimulations;

	for (int i = 0; i < testSim.size(); i++) {
		std::shared_ptr< SimulationParameter> simPara = testSim.at(i)->getSimulationParameter();
		std::shared_ptr< VirtualFluidSimulationImp> vfSim = VirtualFluidSimulationImp::getNewInstance();
		std::shared_ptr< Parameter> para = Parameter::make();

		testSim.at(i)->setParameter(para);
		vfSim->initParameter(para, simPara->getKernelConfiguration() ,simPara->getViscosity(), simPara->getGridPath(), simPara->getFilePath(), simPara->getNumberOfGridLevels(), simPara->getEndTime(), simPara->getTimeStepLength(), simPara->getDevices(), simPara->getMaxVelocity());
		vfSim->initInitialConditions(simPara->getInitialCondition());
		vfSim->initGridProvider();
		vfSim->setDataWriter(testSim.at(i)->getDataWriter());
		vfSim->setTestSimulation(testSim.at(i));
		vfSimulations.push_back(vfSim);		
	}

	return vfSimulations;
}
