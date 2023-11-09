#include "CompressibleNavierStokesStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<CompressibleNavierStokesStrategy> CompressibleNavierStokesStrategy::getInstance()
{
    static std::shared_ptr<CompressibleNavierStokesStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<CompressibleNavierStokesStrategy>(new CompressibleNavierStokesStrategy());
	return uniqueInstance;
}

bool CompressibleNavierStokesStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

CompressibleNavierStokesStrategy::CompressibleNavierStokesStrategy() {
}
