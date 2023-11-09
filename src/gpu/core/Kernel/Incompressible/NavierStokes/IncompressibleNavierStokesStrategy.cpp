#include "IncompressibleNavierStokesStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<IncompressibleNavierStokesStrategy> IncompressibleNavierStokesStrategy::getInstance()
{
    static std::shared_ptr<IncompressibleNavierStokesStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<IncompressibleNavierStokesStrategy>(new IncompressibleNavierStokesStrategy());
	return uniqueInstance;
}

bool IncompressibleNavierStokesStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (para->getUseWale())
		return false;
	else if (para->getCompOn())
		return false;
	else
		return true;
}

IncompressibleNavierStokesStrategy::IncompressibleNavierStokesStrategy() {
}
