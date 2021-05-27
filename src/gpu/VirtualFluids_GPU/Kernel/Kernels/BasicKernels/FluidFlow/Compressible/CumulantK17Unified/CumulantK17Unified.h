#ifndef CUMULANT_K17_UNIFIED_H
#define CUMULANT_K17_UNIFIED_H

#include "Kernel/KernelImp.h"

class CumulantK17Unified : public KernelImp
{
public:
    static std::shared_ptr<CumulantK17Unified> getNewInstance(std::shared_ptr<Parameter> para, int level);
    void run();

    CumulantK17Unified(std::shared_ptr<Parameter> para, int level);
};

#endif
