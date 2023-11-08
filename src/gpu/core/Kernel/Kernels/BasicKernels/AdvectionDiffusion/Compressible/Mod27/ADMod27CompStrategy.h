#ifndef AD_MOD27_COMP_STRATEGY_H
#define AD_MOD27_COMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class ADMod27CompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<ADMod27CompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	ADMod27CompStrategy();

};
#endif 