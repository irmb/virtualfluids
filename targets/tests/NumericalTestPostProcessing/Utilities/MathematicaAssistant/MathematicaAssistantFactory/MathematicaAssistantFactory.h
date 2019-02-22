#ifndef MATHEMATICA_ASSISTANT_FACTORY_H
#define MATHEMATICA_ASSISTANT_FACTORY_H

#include "../MathematicaAssistant.h"

#include <memory>
#include <vector>

class MathematicaFunctionFactory;

class MathematicaAssistantFactory
{
public:
	virtual std::vector<std::shared_ptr<MathematicaAssistant> > makeMathematicaAssistants(std::vector<Assistant> types, std::shared_ptr<MathematicaFunctionFactory> functionFactory) = 0;

};
#endif