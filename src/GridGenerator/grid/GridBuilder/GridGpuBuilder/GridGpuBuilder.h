#ifndef GridGpuBuilder_H
#define GridGpuBuilder_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <vector>
#include <string>
#include <memory>

#include "../GridBuilderImp.h"


class GridGpuBuilder : public GridBuilderImp
{
public:
	GridGenerator_EXPORT GridGpuBuilder();
	virtual GridGenerator_EXPORT ~GridGpuBuilder();

protected:
    void createGridKernels(std::string distribution);

};

#endif

