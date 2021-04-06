#ifndef WALE_FLUID_FLOW_COMP_STRATEGY_H
#define WALE_FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class WaleFluidFlowCompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<WaleFluidFlowCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    WaleFluidFlowCompStrategy();

};
#endif 