#include "TurbulentViscosityFluidFlowCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<TurbulentViscosityFluidFlowCompStrategy> TurbulentViscosityFluidFlowCompStrategy::getInstance()
{
    static std::shared_ptr<TurbulentViscosityFluidFlowCompStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<TurbulentViscosityFluidFlowCompStrategy>(new TurbulentViscosityFluidFlowCompStrategy());
	return uniqueInstance;
}

bool TurbulentViscosityFluidFlowCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (!para->getUseTurbulentViscosity())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

TurbulentViscosityFluidFlowCompStrategy::TurbulentViscosityFluidFlowCompStrategy() {}
