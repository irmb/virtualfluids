#ifndef MATHEMATICA_ASSISTANT_FACTORY_IMP_H
#define MATHEMATICA_ASSISTANT_FACTORY_IMP_H

#include "MathematicaAssistantFactory.h"

class MathematicaAssistantFactoryImp : public MathematicaAssistantFactory
{
public:
    static std::shared_ptr<MathematicaAssistantFactory> getNewInstance();

    std::vector<std::shared_ptr<MathematicaAssistant> > makeMathematicaAssistants(std::vector<Assistant> types, std::shared_ptr<MathematicaFunctionFactory> functionFactory);

private: 
    MathematicaAssistantFactoryImp();
};
#endif