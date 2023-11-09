#ifndef IncompressibleAdvectionDiffusionStrategy_H
#define IncompressibleAdvectionDiffusionStrategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class IncompressibleAdvectionDiffusionStrategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<IncompressibleAdvectionDiffusionStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	IncompressibleAdvectionDiffusionStrategy();

};
#endif 