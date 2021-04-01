#include "PMAdvecCompStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<PMAdvecCompStrategy> PMAdvecCompStrategy::getInstance()
{
	static std::shared_ptr<PMAdvecCompStrategy> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<PMAdvecCompStrategy>(new PMAdvecCompStrategy());
	return uniqueInstance;
}

bool PMAdvecCompStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else if (!para->getSimulatePorousMedia())
		return false;
	else
		return true;
}

PMAdvecCompStrategy::PMAdvecCompStrategy()
{
}
