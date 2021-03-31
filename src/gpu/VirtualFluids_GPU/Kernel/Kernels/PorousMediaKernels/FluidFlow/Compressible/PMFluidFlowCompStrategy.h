#ifndef PM_FLUID_FLOW_COMP_STRATEGY_H
#define PM_FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class PMFluidFlowCompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<PMFluidFlowCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    PMFluidFlowCompStrategy();

};
#endif 