#ifndef CHECK_PARAMETER_STRATEGY_H
#define CHECK_PARAMETER_STRATEGY_H

#include <memory>

class Parameter;

class CheckParameterStrategy
{
public:
    virtual ~CheckParameterStrategy() = default;
    virtual bool checkParameter(std::shared_ptr<Parameter> para) = 0;
};
#endif
