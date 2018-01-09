#ifndef GridCpuBuilder_H
#define GridCpuBuilder_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

#include "../GridBuilderImp.h"


class GridCpuBuilder : public GridBuilderImp
{
public:
    VF_PUBLIC GridCpuBuilder();
    virtual VF_PUBLIC ~GridCpuBuilder();

protected:
    void createGridKernels(std::string distribution);

};

#endif

