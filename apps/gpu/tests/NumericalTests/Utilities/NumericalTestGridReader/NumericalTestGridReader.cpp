#include "NumericalTestGridReader.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"

#include "Utilities/InitialCondition/InitialCondition.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::shared_ptr<NumericalTestGridReader> NumericalTestGridReader::getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager)
{
	return std::shared_ptr<NumericalTestGridReader>(new NumericalTestGridReader(para, initialCondition, cudaManager));
}

void NumericalTestGridReader::setInitalNodeValues(uint numberOfNodes, int level) const
{
	initialCondition->init(level);
	for (uint j = 0; j <= numberOfNodes; j++){
		para->getParH(level)->velocityX[j] = initialCondition->getInitVX(j, level);
		para->getParH(level)->velocityY[j] = initialCondition->getInitVY(j, level);
		para->getParH(level)->velocityZ[j] = initialCondition->getInitVZ(j, level);
		para->getParH(level)->rho[j] = initialCondition->getInitROH(j, level);
		para->getParH(level)->pressure[j] = initialCondition->getInitPRESS(j, level);
	}
}

NumericalTestGridReader::NumericalTestGridReader(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager) : GridReader(FILEFORMAT::BINARY, para, cudaManager), initialCondition(initialCondition)
{

}