#include "AMDFluidFlowCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<AMDFluidFlowCompStrategy> AMDFluidFlowCompStrategy::getInstance()
{
    static std::shared_ptr<AMDFluidFlowCompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<AMDFluidFlowCompStrategy>(new AMDFluidFlowCompStrategy());
	return uniqueInstance;
}

bool AMDFluidFlowCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (!para->getUseAMD())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

AMDFluidFlowCompStrategy::AMDFluidFlowCompStrategy() {
}
