#ifndef AD_MOD7_COMP_STRATEGY_H
#define AD_MOD7_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class ADMod7CompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<ADMod7CompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	ADMod7CompStrategy();

};
#endif 