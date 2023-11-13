#ifndef Init_COMP_AD_27_H
#define Init_COMP_AD_27_H

#include "PreProcessor/PreProcessorStrategy/PreProcessorStrategy.h"

#include <memory>

class Parameter;

class InitCompAD27 : public PreProcessorStrategy
{
public:
    static std::shared_ptr<PreProcessorStrategy> getNewInstance(std::shared_ptr< Parameter> para);
    void init(int level);
    bool checkParameter();

private:
    InitCompAD27();
    InitCompAD27(std::shared_ptr< Parameter> para);
    std::shared_ptr< Parameter> para;
};

#endif 