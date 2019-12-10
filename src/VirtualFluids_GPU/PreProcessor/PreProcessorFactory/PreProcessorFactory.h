#ifndef PREPROCESSOR_FACTORY_H
#define PREPROCESSOR_FACTORY_H

#include "PreProcessor/PreProcessorType.h"

#include <memory>
#include <vector>

class PreProcessor;
class Parameter;

class VF_PUBLIC PreProcessorFactory
{
public:
	virtual std::shared_ptr<PreProcessor> makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para) = 0;

};
#endif