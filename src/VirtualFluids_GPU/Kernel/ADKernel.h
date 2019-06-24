#ifndef ADV_DIFF_KERNEL_H
#define ADV_DIFF_KERNEL_H

#include "Kernel\KernelImp.h"

class ADKernel : public KernelImp
{
	virtual void run() = 0;
};
#endif 