#ifndef GPU_BKGUnified_H
#define GPU_BKGUnified_H

#include "Kernel/KernelImp.h"

class BGKUnified : public KernelImp
{
public:
    static std::shared_ptr<BGKUnified> getNewInstance(std::shared_ptr<Parameter> para, int level);
    void run();

    BGKUnified(std::shared_ptr<Parameter> para, int level);
};

#endif
