#ifndef AMD_FLUID_FLOW_COMP_STRATEGY_H
#define AMD_FLUID_FLOW_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class TurbulentViscosityFluidFlowCompStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<TurbulentViscosityFluidFlowCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    TurbulentViscosityFluidFlowCompStrategy();

};
#endif 