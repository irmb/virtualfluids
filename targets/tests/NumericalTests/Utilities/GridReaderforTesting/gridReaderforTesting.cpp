#include "gridReaderforTesting.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"

#include "Utilities/InitialCondition/InitialCondition.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::shared_ptr<GridReaderforTesting> GridReaderforTesting::getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition)
{
	return std::shared_ptr<GridReaderforTesting>(new GridReaderforTesting(para, initialCondition));
}

void GridReaderforTesting::setInitalNodeValues(const int numberOfNodes, const int level) const
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

GridReaderforTesting::GridReaderforTesting(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition) :GridReader(true, para), initialCondition(initialCondition)
{

}