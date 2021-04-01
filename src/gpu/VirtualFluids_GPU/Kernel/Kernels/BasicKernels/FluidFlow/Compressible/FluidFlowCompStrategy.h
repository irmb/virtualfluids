#ifndef FLUID_FLOW_COMP_STRATEGY_H
#define FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class FluidFlowCompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<FluidFlowCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    FluidFlowCompStrategy();

};
#endif 