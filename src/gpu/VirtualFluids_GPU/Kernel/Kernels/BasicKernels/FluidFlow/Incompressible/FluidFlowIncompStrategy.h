#ifndef FLUID_FLOW_INCOMP_STRATEGY_H
#define FLUID_FLOW_INCOMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class FluidFlowIncompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<FluidFlowIncompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    FluidFlowIncompStrategy();

};
#endif 