#ifndef CUMULANT_K15_UNIFIED_COMP_H
#define CUMULANT_K15_UNIFIED_COMP_H

#include "Kernel/KernelImp.h"

namespace vf
{
namespace gpu
{
class CumulantK15Unified : public KernelImp
{
public:
    CumulantK15Unified(std::shared_ptr<Parameter> para, int level);
    
    void run();
};

}
}

#endif
