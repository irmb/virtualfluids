#include "CompressibleNavierStokesWaleStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<CompressibleNavierStokesWaleStrategy> CompressibleNavierStokesWaleStrategy::getInstance()
{
    static std::shared_ptr<CompressibleNavierStokesWaleStrategy> uniqueInstance;
	if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<CompressibleNavierStokesWaleStrategy>(new CompressibleNavierStokesWaleStrategy());
	return uniqueInstance;
}

bool CompressibleNavierStokesWaleStrategy::checkParameter(std::shared_ptr<Parameter> para)
{
	if (!para->getUseWale())
		return false;
	else if (!para->getCompOn())
		return false;
	else
		return true;
}

CompressibleNavierStokesWaleStrategy::CompressibleNavierStokesWaleStrategy() {
}
