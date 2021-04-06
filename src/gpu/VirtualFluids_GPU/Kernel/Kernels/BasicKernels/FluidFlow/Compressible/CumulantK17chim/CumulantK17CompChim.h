#ifndef CUMULANT_K17_COMP_CHIM_H
#define CUMULANT_K17_COMP_CHIM_H

#include "Kernel/KernelImp.h"

class CumulantK17CompChim : public KernelImp
{
public:
	static std::shared_ptr<CumulantK17CompChim> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
    CumulantK17CompChim();
    CumulantK17CompChim(std::shared_ptr<Parameter> para, int level);
};

#endif 
