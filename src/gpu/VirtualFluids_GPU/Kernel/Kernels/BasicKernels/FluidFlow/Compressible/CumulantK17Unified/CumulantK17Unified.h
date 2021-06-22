#ifndef CUMULANT_K17_UNIFIED_H
#define CUMULANT_K17_UNIFIED_H

#include "Kernel/KernelImp.h"

namespace vf
{
namespace gpu
{


class CumulantK17Unified : public KernelImp
{
public:
    CumulantK17Unified(std::shared_ptr<Parameter> para, int level);

    void run();
};

}
}

#endif
