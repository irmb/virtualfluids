#ifndef K15CompressibleNavierStokesUnified_H
#define K15CompressibleNavierStokesUnified_H

#include "Kernel/KernelImp.h"

namespace vf
{
namespace gpu
{
class K15CompressibleNavierStokesUnified : public KernelImp
{
public:
    K15CompressibleNavierStokesUnified(std::shared_ptr<Parameter> para, int level);
    
    void run();
};

}
}

#endif
