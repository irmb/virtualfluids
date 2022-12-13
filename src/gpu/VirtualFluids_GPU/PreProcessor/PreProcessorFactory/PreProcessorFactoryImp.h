#ifndef PREPROCESSOR_FACTORY_IMP_H
#define PREPROCESSOR_FACTORY_IMP_H

#include "PreProcessorFactory.h"

class PreProcessorStrategy;

class VIRTUALFLUIDS_GPU_EXPORT PreProcessorFactoryImp : public PreProcessorFactory
{
public:
	std::shared_ptr<PreProcessor> makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para);

	std::shared_ptr<PreProcessorStrategy> makePreProcessorStrategy(PreProcessorType preProcessorTyp, std::shared_ptr<Parameter> para);

};
#endif