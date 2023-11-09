#include "IncompressibleAdvectionDiffusionStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<IncompressibleAdvectionDiffusionStrategy> IncompressibleAdvectionDiffusionStrategy::getInstance()
{
	static std::shared_ptr<IncompressibleAdvectionDiffusionStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<IncompressibleAdvectionDiffusionStrategy>(new IncompressibleAdvectionDiffusionStrategy());
	return uniqueInstance;
}

bool IncompressibleAdvectionDiffusionStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (para->getCompOn())
		return false;
	else if (!para->getDiffOn())
		return false;
	else if (para->getDiffMod() == 7)
		return false;
	else
		return true;
}

IncompressibleAdvectionDiffusionStrategy::IncompressibleAdvectionDiffusionStrategy()
{
}
