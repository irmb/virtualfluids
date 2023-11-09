#ifndef CompressibleAdvectionDiffusionStrategy_H
#define CompressibleAdvectionDiffusionStrategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class CompressibleAdvectionDiffusionStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<CompressibleAdvectionDiffusionStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	CompressibleAdvectionDiffusionStrategy();

};
#endif 