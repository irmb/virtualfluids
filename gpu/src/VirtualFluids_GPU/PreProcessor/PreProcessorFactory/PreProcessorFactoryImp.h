#ifndef PREPROCESSOR_FACTORY_IMP_H
#define PREPROCESSOR_FACTORY_IMP_H

#include "PreProcessorFactory.h"

class PreProcessorStrategy;

class VF_PUBLIC PreProcessorFactoryImp : public PreProcessorFactory
{
public:
	static std::shared_ptr< PreProcessorFactoryImp> getInstance();
	std::shared_ptr<PreProcessor> makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para);

private:
	PreProcessorFactoryImp();

	std::shared_ptr<PreProcessorStrategy> makePreProcessorStrategy(PreProcessorType preProcessorTyp, std::shared_ptr<Parameter> para);

};
#endif