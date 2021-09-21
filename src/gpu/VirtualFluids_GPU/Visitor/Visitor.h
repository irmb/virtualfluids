#ifndef Visitor_H
#define Visitor_H

#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"
#include "VirtualFluids_GPU_export.h"

#include <cassert>

class Parameter;
class GridProvider;
class CudaMemoryManager;

class VIRTUALFLUIDS_GPU_EXPORT Visitor
{
private:
    SPtr<Parameter> para;

protected:
    Visitor()
    {
        this->updateInterval = 1;
    }

public:
    virtual ~Visitor()
    {
    }

    virtual void init(Parameter *para, GridProvider *gridProvider, CudaMemoryManager *cudaManager) = 0;
    virtual void visit(Parameter *para, CudaMemoryManager *cudaManager, int level, uint t) = 0;
    virtual void free(Parameter *para, CudaMemoryManager *cudaManager) = 0;

protected:
    uint updateInterval;
};

#endif