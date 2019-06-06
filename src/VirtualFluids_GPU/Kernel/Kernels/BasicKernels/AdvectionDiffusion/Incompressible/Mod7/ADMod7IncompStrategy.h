#ifndef AD_MOD7_INCOMP_STRATEGY_H
#define AD_MOD7_INCOMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class ADMod7IncompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<ADMod7IncompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	ADMod7IncompStrategy();

};
#endif 