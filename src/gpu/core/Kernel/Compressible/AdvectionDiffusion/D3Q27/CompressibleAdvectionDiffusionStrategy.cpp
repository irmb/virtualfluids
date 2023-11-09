#include "CompressibleAdvectionDiffusionStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<CompressibleAdvectionDiffusionStrategy> CompressibleAdvectionDiffusionStrategy::getInstance()
{
	std::shared_ptr<CompressibleAdvectionDiffusionStrategy> uniqueInstance;
	if(!uniqueInstance)
		uniqueInstance = std::shared_ptr<CompressibleAdvectionDiffusionStrategy>(new CompressibleAdvectionDiffusionStrategy());
	return uniqueInstance;
}

bool CompressibleAdvectionDiffusionStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else if (!para->getDiffOn())
		return false;
	else if (para->getDiffMod() == 7)
		return false;
	else
		return true;
}

CompressibleAdvectionDiffusionStrategy::CompressibleAdvectionDiffusionStrategy()
{
}
