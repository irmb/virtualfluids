#include "WaleFluidFlowCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<WaleFluidFlowCompStrategy> WaleFluidFlowCompStrategy::getInstance()
{
    static std::shared_ptr<WaleFluidFlowCompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<WaleFluidFlowCompStrategy>(new WaleFluidFlowCompStrategy());
	return uniqueInstance;
}

bool WaleFluidFlowCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (!para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

WaleFluidFlowCompStrategy::WaleFluidFlowCompStrategy() {
}
