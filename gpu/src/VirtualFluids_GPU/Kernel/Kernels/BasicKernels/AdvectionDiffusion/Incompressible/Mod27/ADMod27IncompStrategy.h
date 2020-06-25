#ifndef AD_MOD27_INCOMP_STRATEGY_H
#define AD_MOD27_INCOMP_STRATEGY_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class ADMod27IncompStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<ADMod27IncompStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	ADMod27IncompStrategy();

};
#endif 