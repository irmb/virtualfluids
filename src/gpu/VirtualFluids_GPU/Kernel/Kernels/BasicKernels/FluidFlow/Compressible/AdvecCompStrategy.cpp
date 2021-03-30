#include "AdvecCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<AdvecCompStrategy> AdvecCompStrategy::getInstance()
{
	static std::shared_ptr<AdvecCompStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<AdvecCompStrategy>(new AdvecCompStrategy());
	return uniqueInstance;
}

bool AdvecCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

AdvecCompStrategy::AdvecCompStrategy()
{
}
