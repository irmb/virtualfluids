#ifndef CUMULANT_K15_INCOMP_H
#define CUMULANT_K15_INCOMP_H

#include "Kernel\KernelImp.h"

class CumulantK15Incomp : public KernelImp
{
public:
	static std::shared_ptr<CumulantK15Incomp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK15Incomp();
	CumulantK15Incomp(std::shared_ptr< Parameter> para, int level);
};
#endif 