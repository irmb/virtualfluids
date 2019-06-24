#ifndef CUMULANT_1H_INCOMP_SP27_H
#define CUMULANT_1H_INCOMP_SP27_H

#include "Kernel\KernelImp.h"

class Cumulant1hIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<Cumulant1hIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	Cumulant1hIncompSP27();
	Cumulant1hIncompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 