#ifndef PREPROCESSOR_FACTORY_H
#define PREPROCESSOR_FACTORY_H

#include "PreProcessor/PreProcessorType.h"

#include <memory>
#include <vector>



class PreProcessor;
class Parameter;

class PreProcessorFactory
{
public:
    virtual ~PreProcessorFactory() = default;
	virtual std::shared_ptr<PreProcessor> makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para) = 0;

};
#endif