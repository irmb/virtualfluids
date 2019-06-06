#include "ADMod27IncompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<ADMod27IncompStrategy> ADMod27IncompStrategy::getInstance()
{
	static std::shared_ptr<ADMod27IncompStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<ADMod27IncompStrategy>(new ADMod27IncompStrategy());
	return uniqueInstance;
}

bool ADMod27IncompStrategy::checkParameter(std::shared_ptr<Parameter> para)
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

ADMod27IncompStrategy::ADMod27IncompStrategy()
{
}
