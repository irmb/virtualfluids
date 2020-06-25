#include "ADMod7CompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<ADMod7CompStrategy> ADMod7CompStrategy::getInstance()
{
	static std::shared_ptr<ADMod7CompStrategy> uniqueInstance;
	if(!uniqueInstance)
		uniqueInstance = std::shared_ptr<ADMod7CompStrategy>(new ADMod7CompStrategy());
	return uniqueInstance;
}

bool ADMod7CompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else if (!para->getDiffOn())
		return false;
	else if (para->getDiffMod() == 27)
		return false;
	else
		return true;
}

ADMod7CompStrategy::ADMod7CompStrategy()
{
}
