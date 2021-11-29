#ifndef AMDCUMULANT_K17_COMP_CHIM_H
#define AMDCUMULANT_K17_COMP_CHIM_H

#include "Kernel/KernelImp.h"

class AMDCumulantK17CompChim : public KernelImp
{
public:
	static std::shared_ptr<AMDCumulantK17CompChim> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
    AMDCumulantK17CompChim();
    AMDCumulantK17CompChim(std::shared_ptr<Parameter> para, int level);
};

#endif 
