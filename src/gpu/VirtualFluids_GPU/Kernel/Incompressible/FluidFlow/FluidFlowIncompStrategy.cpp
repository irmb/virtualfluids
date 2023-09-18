#include "FluidFlowIncompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<FluidFlowIncompStrategy> FluidFlowIncompStrategy::getInstance()
{
    static std::shared_ptr<FluidFlowIncompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<FluidFlowIncompStrategy>(new FluidFlowIncompStrategy());
	return uniqueInstance;
}

bool FluidFlowIncompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (para->getCompOn())
		return false;
	else
		return true;
}

FluidFlowIncompStrategy::FluidFlowIncompStrategy() {
}
