#ifndef PREPROCESSOR_STRATEGY_H
#define PREPROCESSOR_STRATEGY_H

#include <DataTypes.h>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>


class PreProcessorStrategy
{
public:
    virtual void init(int level) = 0;
    virtual bool checkParameter() = 0;
};
#endif