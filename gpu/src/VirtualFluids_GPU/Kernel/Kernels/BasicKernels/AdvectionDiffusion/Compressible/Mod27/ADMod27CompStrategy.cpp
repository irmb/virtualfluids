#include "ADMod27CompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<ADMod27CompStrategy> ADMod27CompStrategy::getInstance()
{
	std::shared_ptr<ADMod27CompStrategy> uniqueInstance;
	if(!uniqueInstance)
		uniqueInstance = std::shared_ptr<ADMod27CompStrategy>(new ADMod27CompStrategy());
	return uniqueInstance;
}

bool ADMod27CompStrategy::checkParameter(std::shared_ptr<Parameter> para)
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

ADMod27CompStrategy::ADMod27CompStrategy()
{
}
