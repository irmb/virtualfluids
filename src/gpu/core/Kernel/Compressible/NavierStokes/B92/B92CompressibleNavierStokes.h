#ifndef B92CompressibleNavierStokes_H
#define B92CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class B92CompressibleNavierStokes : public KernelImp
{
public:
    static std::shared_ptr<B92CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
    void run();
    
private:
    B92CompressibleNavierStokes();
    B92CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 