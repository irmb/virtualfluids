#ifndef K15CompressibleNavierStokes_H
#define K15CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class K15CompressibleNavierStokes : public KernelImp
{
public:
    static std::shared_ptr< K15CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
    void run();

private:
    K15CompressibleNavierStokes();
    K15CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 