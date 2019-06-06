#include "WaleAdvecCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<WaleAdvecCompStrategy> WaleAdvecCompStrategy::getInstance()
{
	static std::shared_ptr<WaleAdvecCompStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<WaleAdvecCompStrategy>(new WaleAdvecCompStrategy());
	return uniqueInstance;
}

bool WaleAdvecCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (!para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

WaleAdvecCompStrategy::WaleAdvecCompStrategy()
{
}
