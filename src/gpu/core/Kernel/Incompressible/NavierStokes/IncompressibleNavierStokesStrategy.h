#ifndef FLUID_FLOW_INCOMP_STRATEGY_H
#define FLUID_FLOW_INCOMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class IncompressibleNavierStokesStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<IncompressibleNavierStokesStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    IncompressibleNavierStokesStrategy();

};
#endif 