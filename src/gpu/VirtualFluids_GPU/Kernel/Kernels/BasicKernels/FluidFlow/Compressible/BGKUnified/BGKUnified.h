#ifndef GPU_BKGUnified_H
#define GPU_BKGUnified_H

#include "Kernel/KernelImp.h"

namespace vf
{
namespace gpu
{

class BGKUnified : public KernelImp
{
public:
    BGKUnified(std::shared_ptr<Parameter> para, int level);

    void run();
};

}
}

#endif
