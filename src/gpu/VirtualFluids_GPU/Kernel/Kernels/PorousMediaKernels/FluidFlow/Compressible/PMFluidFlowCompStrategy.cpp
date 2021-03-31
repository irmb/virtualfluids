#include "PMFluidFlowCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<PMFluidFlowCompStrategy> PMFluidFlowCompStrategy::getInstance()
{
    static std::shared_ptr<PMFluidFlowCompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<PMFluidFlowCompStrategy>(new PMFluidFlowCompStrategy());
	return uniqueInstance;
}

bool PMFluidFlowCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else if (!para->getSimulatePorousMedia())
		return false;
	else
		return true;
}

PMFluidFlowCompStrategy::PMFluidFlowCompStrategy(){
}
