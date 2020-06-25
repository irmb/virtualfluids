#ifndef WALE_ADVEC_COMP_STRATEGY_H
#define WALE_ADVEC_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class WaleAdvecCompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<WaleAdvecCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	WaleAdvecCompStrategy();

};
#endif 