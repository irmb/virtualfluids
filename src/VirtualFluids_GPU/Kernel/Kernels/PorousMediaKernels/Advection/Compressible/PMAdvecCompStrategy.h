#ifndef PM_ADVEC_COMP_STRATEGY_H
#define PM_ADVEC_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class PMAdvecCompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<PMAdvecCompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	PMAdvecCompStrategy();

};
#endif 