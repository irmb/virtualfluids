#ifndef ADVEC_COMP_STRATEGY_H
#define ADVEC_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class AdvecCompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<AdvecCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	AdvecCompStrategy();

};
#endif 