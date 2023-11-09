#ifndef CompressibleAdvectionDiffusionD3Q7Strategy_H
#define CompressibleAdvectionDiffusionD3Q7Strategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class CompressibleAdvectionDiffusionD3Q7Strategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<CompressibleAdvectionDiffusionD3Q7Strategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	CompressibleAdvectionDiffusionD3Q7Strategy();

};
#endif 