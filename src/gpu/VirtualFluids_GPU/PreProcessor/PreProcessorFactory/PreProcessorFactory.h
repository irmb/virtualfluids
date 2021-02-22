#ifndef PREPROCESSOR_FACTORY_H
#define PREPROCESSOR_FACTORY_H

#include "PreProcessor/PreProcessorType.h"

#include <memory>
#include <vector>

#include "VirtualFluids_GPU_export.h"

class PreProcessor;
class Parameter;

class VIRTUALFLUIDS_GPU_EXPORT PreProcessorFactory
{
public:
    virtual ~PreProcessorFactory() = default;
	virtual std::shared_ptr<PreProcessor> makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para) = 0;

};
#endif