#ifndef CUMULANT_K15_UNIFIED_COMP_H
#define CUMULANT_K15_UNIFIED_COMP_H

#include "Kernel/KernelImp.h"

class CumulantK15Unified : public KernelImp
{
public:
    static std::shared_ptr<CumulantK15Unified> getNewInstance(std::shared_ptr<Parameter> para, int level);
    void run();

    CumulantK15Unified(std::shared_ptr<Parameter> para, int level);
};
#endif
