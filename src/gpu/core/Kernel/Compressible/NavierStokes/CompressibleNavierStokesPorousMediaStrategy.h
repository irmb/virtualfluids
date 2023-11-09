#ifndef CompressibleNavierStokesPorousMediaStrategy_H
#define CompressibleNavierStokesPorousMediaStrategy_H

#include "Kernel/Utilities/CheckParameterStrategy/CheckParameterStrategy.h"


class CompressibleNavierStokesPorousMediaStrategy : public CheckParameterStrategy
{
public:
    static std::shared_ptr<CompressibleNavierStokesPorousMediaStrategy> getInstance();

	bool checkParameter(std::shared_ptr<Parameter> para);

private:
    CompressibleNavierStokesPorousMediaStrategy();

};
#endif 