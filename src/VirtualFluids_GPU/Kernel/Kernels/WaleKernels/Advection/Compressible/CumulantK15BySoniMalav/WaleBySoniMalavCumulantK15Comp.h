#ifndef WALE_CUMULANT_K15_COMP_BY_SONI_MALAV_H
#define WALE_CUMULANT_K15_COMP_BY_SONI_MALAV_H

#include "Kernel\KernelImp.h"

class WaleBySoniMalavCumulantK15Comp : public KernelImp
{
public:
	static std::shared_ptr<WaleBySoniMalavCumulantK15Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleBySoniMalavCumulantK15Comp();
	WaleBySoniMalavCumulantK15Comp(std::shared_ptr< Parameter> para, int level);

};
#endif 