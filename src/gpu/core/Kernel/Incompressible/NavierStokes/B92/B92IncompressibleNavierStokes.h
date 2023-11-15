#ifndef B92IncompressibleNavierStokes_H
#define B92IncompressibleNavierStokes_H

#include "Kernel/KernelImp.h"


class B92IncompressibleNavierStokes : public KernelImp
{
public:
    static std::shared_ptr<B92IncompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
    void run();

private:
    B92IncompressibleNavierStokes();
    B92IncompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 