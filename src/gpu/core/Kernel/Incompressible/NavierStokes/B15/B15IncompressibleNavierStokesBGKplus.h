#ifndef B15IncompressibleNavierStokesBGKplus_H
#define B15IncompressibleNavierStokesBGKplus_H

#include "Kernel/KernelImp.h"

class B15IncompressibleNavierStokesBGKplus : public KernelImp
{
public:
    static std::shared_ptr<B15IncompressibleNavierStokesBGKplus> getNewInstance(std::shared_ptr< Parameter> para, int level);
    void run();

private:
    B15IncompressibleNavierStokesBGKplus();
    B15IncompressibleNavierStokesBGKplus(std::shared_ptr< Parameter> para, int level);

};
#endif 