#include "ADMod7IncompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<ADMod7IncompStrategy> ADMod7IncompStrategy::getInstance()
{
	static std::shared_ptr<ADMod7IncompStrategy> uniqueInstance;
	if(!uniqueInstance)
		uniqueInstance = std::shared_ptr<ADMod7IncompStrategy>(new ADMod7IncompStrategy());
	return uniqueInstance;
}

bool ADMod7IncompStrategy::checkParameter(std::shared_ptr<Parameter> para)
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

ADMod7IncompStrategy::ADMod7IncompStrategy()
{
}