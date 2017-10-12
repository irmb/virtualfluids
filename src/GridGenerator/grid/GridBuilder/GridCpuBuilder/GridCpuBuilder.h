#ifndef GridCpuBuilder_H
#define GridCpuBuilder_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <vector>
#include <string>
#include <memory>

#include "../GridBuilderImp.h"


class GridCpuBuilder : public GridBuilderImp
{
public:
    GridGenerator_EXPORT GridCpuBuilder();
    virtual GridGenerator_EXPORT ~GridCpuBuilder();

protected:
    void createGridKernels(std::string distribution);

};

#endif

