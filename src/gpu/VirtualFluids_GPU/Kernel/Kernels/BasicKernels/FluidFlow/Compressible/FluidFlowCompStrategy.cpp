#include "FluidFlowCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<FluidFlowCompStrategy> FluidFlowCompStrategy::getInstance()
{
    static std::shared_ptr<FluidFlowCompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<FluidFlowCompStrategy>(new FluidFlowCompStrategy());
	return uniqueInstance;
}

bool FluidFlowCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

FluidFlowCompStrategy::FluidFlowCompStrategy() {
}
