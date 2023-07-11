#ifndef B15CompressibleNavierStokesBGKplusUnified_H
#define B15CompressibleNavierStokesBGKplusUnified_H

#include "Kernel/KernelImp.h"

namespace vf
{
namespace gpu
{

class B15CompressibleNavierStokesBGKplusUnified : public KernelImp
{
public:
    B15CompressibleNavierStokesBGKplusUnified(std::shared_ptr<Parameter> para, int level);

    void run();
};

}
}

#endif
