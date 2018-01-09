#ifndef GridGpuBuilder_H
#define GridGpuBuilder_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

#include "../GridBuilderImp.h"


class GridGpuBuilder : public GridBuilderImp
{
public:
	VF_PUBLIC GridGpuBuilder();
	virtual VF_PUBLIC ~GridGpuBuilder();

protected:
    void createGridKernels(std::string distribution);

};

#endif

