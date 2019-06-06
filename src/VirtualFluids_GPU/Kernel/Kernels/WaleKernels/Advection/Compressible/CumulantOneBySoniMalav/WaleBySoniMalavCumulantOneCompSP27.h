#ifndef WALE_CUMULANT_ONE_COMP_BY_SONI_MALAV_SP27_H
#define WALE_CUMULANT_ONE_COMP_BY_SONI_MALAV_SP27_H

#include "Kernel\KernelImp.h"

class WaleBySoniMalavCumulantOneCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<WaleBySoniMalavCumulantOneCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleBySoniMalavCumulantOneCompSP27();
	WaleBySoniMalavCumulantOneCompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 