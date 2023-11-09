#include "IncompressibleAdvectionDiffusionD3Q7Strategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<IncompressibleAdvectionDiffusionD3Q7Strategy> IncompressibleAdvectionDiffusionD3Q7Strategy::getInstance()
{
	static std::shared_ptr<IncompressibleAdvectionDiffusionD3Q7Strategy> uniqueInstance;
	if(!uniqueInstance)
		uniqueInstance = std::shared_ptr<IncompressibleAdvectionDiffusionD3Q7Strategy>(new IncompressibleAdvectionDiffusionD3Q7Strategy());
	return uniqueInstance;
}

bool IncompressibleAdvectionDiffusionD3Q7Strategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (para->getCompOn())
		return false;
	else if (!para->getDiffOn())
		return false;
	else if (para->getDiffMod() == 27)
		return false;
	else
		return true;
}

IncompressibleAdvectionDiffusionD3Q7Strategy::IncompressibleAdvectionDiffusionD3Q7Strategy()
{
}