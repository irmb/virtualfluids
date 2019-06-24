#include "AdvecIncompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<AdvecIncompStrategy> AdvecIncompStrategy::getInstance()
{
	static std::shared_ptr<AdvecIncompStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<AdvecIncompStrategy>(new AdvecIncompStrategy());
	return uniqueInstance;
}

bool AdvecIncompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (para->getCompOn())
		return false;
	else
		return true;
}

AdvecIncompStrategy::AdvecIncompStrategy()
{
}
