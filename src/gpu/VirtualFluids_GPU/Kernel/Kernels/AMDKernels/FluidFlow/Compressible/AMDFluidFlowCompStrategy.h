#ifndef AMD_FLUID_FLOW_COMP_STRATEGY_H
#define AMD_FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class AMDFluidFlowCompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<AMDFluidFlowCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    AMDFluidFlowCompStrategy();

};
#endif 