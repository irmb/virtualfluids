#include "NumericalTestGridReader.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"

#include "Utilities/InitialCondition/InitialCondition.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::shared_ptr<NumericalTestGridReader> NumericalTestGridReader::getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager)
{
	return std::shared_ptr<NumericalTestGridReader>(new NumericalTestGridReader(para, initialCondition, cudaManager));
}

void NumericalTestGridReader::setInitalNodeValues(const int numberOfNodes, const int level) const
{
	initialCondition->init(level);
	for (int j = 0; j <= numberOfNodes; j++){
		para->getParH(level)->vx_SP[j] = initialCondition->getInitVX(j, level);
		para->getParH(level)->vy_SP[j] = initialCondition->getInitVY(j, level);
		para->getParH(level)->vz_SP[j] = initialCondition->getInitVZ(j, level);
		para->getParH(level)->rho_SP[j] = initialCondition->getInitROH(j, level);
		para->getParH(level)->press_SP[j] = initialCondition->getInitPRESS(j, level);
	}
}

NumericalTestGridReader::NumericalTestGridReader(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager) : GridReader(true, para, cudaManager), initialCondition(initialCondition)
{

}