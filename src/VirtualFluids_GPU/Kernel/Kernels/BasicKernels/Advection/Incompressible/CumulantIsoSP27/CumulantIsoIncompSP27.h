#ifndef CUMULANT_ISO_INCOMP_SP27_H
#define CUMULANT_ISO_INCOMP_SP27_H

#include "Kernel/KernelImp.h"

class CumulantIsoIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantIsoIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantIsoIncompSP27();
	CumulantIsoIncompSP27(std::shared_ptr< Parameter> para, int level);

};

#endif 