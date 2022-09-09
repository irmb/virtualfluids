#ifndef TurbulentViscosityCUMULANT_K17_COMP_CHIM_H
#define TurbulentViscosityCUMULANT_K17_COMP_CHIM_H

#include "Kernel/KernelImp.h"
#include "Parameter/Parameter.h"

template<TurbulenceModel turbulenceModel> 
class TurbulentViscosityCumulantK17CompChim : public KernelImp
{
public:
	static std::shared_ptr< TurbulentViscosityCumulantK17CompChim<turbulenceModel> > getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
    TurbulentViscosityCumulantK17CompChim();
    TurbulentViscosityCumulantK17CompChim(std::shared_ptr<Parameter> para, int level);
};

#endif 
