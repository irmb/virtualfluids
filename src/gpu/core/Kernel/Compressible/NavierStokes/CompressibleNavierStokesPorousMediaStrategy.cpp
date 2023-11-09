#include "CompressibleNavierStokesPorousMediaStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<CompressibleNavierStokesPorousMediaStrategy> CompressibleNavierStokesPorousMediaStrategy::getInstance()
{
    static std::shared_ptr<CompressibleNavierStokesPorousMediaStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<CompressibleNavierStokesPorousMediaStrategy>(new CompressibleNavierStokesPorousMediaStrategy());
	return uniqueInstance;
}

bool CompressibleNavierStokesPorousMediaStrategy::checkParameter(std::shared_ptr<Parameter> para)
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

CompressibleNavierStokesPorousMediaStrategy::CompressibleNavierStokesPorousMediaStrategy(){
}
