#ifndef PreCollisionInteractor_H
#define PreCollisionInteractor_H

#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

#include <cassert>

class Parameter;
class GridProvider;
class CudaMemoryManager;

class VIRTUALFLUIDS_GPU_EXPORT PreCollisionInteractor
{
private:
    SPtr<Parameter> para;

protected:
    PreCollisionInteractor()
    {
        this->updateInterval = 1;
    }

public:
    virtual ~PreCollisionInteractor()
    {
    }

    virtual void init(Parameter *para, GridProvider *gridProvider, CudaMemoryManager *cudaMemoryManager) = 0;
    virtual void interact(Parameter *para, CudaMemoryManager *cudaMemoryManager, int level, uint t) = 0;
    virtual void free(Parameter *para, CudaMemoryManager *cudaMemoryManager) = 0;

protected:
    uint updateInterval;
};

#endif