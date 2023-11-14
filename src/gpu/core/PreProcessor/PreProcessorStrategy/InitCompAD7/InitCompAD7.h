#ifndef INIT_COMP_AD_7_H
#define INIT_COMP_AD_7_H

#include "PreProcessor/PreProcessorStrategy/PreProcessorStrategy.h"

#include <memory>

class Parameter;

class InitCompAD7 : public PreProcessorStrategy
{
public:
    static std::shared_ptr<InitCompAD7> getNewInstance(std::shared_ptr< Parameter> para);
    void init(int level);
    bool checkParameter();

private:
    InitCompAD7();
    InitCompAD7(std::shared_ptr< Parameter> para);
    std::shared_ptr< Parameter> para;
};
#endif 