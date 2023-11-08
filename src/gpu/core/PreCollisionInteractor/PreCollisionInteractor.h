#ifndef PreCollisionInteractor_H
#define PreCollisionInteractor_H

#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>


#include <cassert>

class Parameter;
class GridProvider;
class CudaMemoryManager;

class PreCollisionInteractor
{

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
    virtual void getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider) = 0;

protected:
    uint updateInterval;
};

#endif