#ifndef IncompressibleAdvectionDiffusionD3Q7Strategy_H
#define IncompressibleAdvectionDiffusionD3Q7Strategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"

class IncompressibleAdvectionDiffusionD3Q7Strategy : public CheckParameterStrategy
{
public:
	static std::shared_ptr<IncompressibleAdvectionDiffusionD3Q7Strategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
	IncompressibleAdvectionDiffusionD3Q7Strategy();

};
#endif 