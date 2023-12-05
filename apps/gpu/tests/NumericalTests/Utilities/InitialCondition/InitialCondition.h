#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include "gpu/core/Calculation/Calculation.h"

#include <vector>
#include <memory>

class Parameter;

class InitialCondition
{
public:
    virtual ~InitialCondition() = default;
    virtual void setParameter(std::shared_ptr<Parameter> para) = 0;
    virtual void init(const int level) = 0;
    virtual real getInitVX(int i, int level) = 0;
    virtual real getInitVY(int i, int level) = 0;
    virtual real getInitVZ(int i, int level) = 0;
    virtual real getInitROH(int i, int level) = 0;
    virtual real getInitPRESS(int i, int level) = 0;

};
#endif