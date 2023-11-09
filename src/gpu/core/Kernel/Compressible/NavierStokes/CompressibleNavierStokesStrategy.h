#ifndef FLUID_FLOW_COMP_STRATEGY_H
#define FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class CompressibleNavierStokesStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<CompressibleNavierStokesStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    CompressibleNavierStokesStrategy();

};
#endif 
