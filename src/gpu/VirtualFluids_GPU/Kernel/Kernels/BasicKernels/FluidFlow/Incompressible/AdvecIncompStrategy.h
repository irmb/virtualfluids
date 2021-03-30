#ifndef ADVEC_INCOMP_STRATEGY_H
#define ADVEC_INCOMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class AdvecIncompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<AdvecIncompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	AdvecIncompStrategy();

};
#endif 