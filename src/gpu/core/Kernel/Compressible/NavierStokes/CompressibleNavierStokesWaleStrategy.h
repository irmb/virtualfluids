#ifndef CompressibleNavierStokesWaleStrategy_H
#define CompressibleNavierStokesWaleStrategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class CompressibleNavierStokesWaleStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<CompressibleNavierStokesWaleStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    CompressibleNavierStokesWaleStrategy();

};
#endif 