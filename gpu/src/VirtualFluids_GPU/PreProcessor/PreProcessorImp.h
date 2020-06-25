#ifndef PREPROCESSOR_IMP_H
#define PREPROCESSOR_IMP_H

#include "PreProcessor.h"

#include <memory>
#include <vector>

class PreProcessorStrategy;

class PreProcessorImp : public PreProcessor
{
public:
	static std::shared_ptr<PreProcessorImp> getNewInstance();
	void addStrategy(std::shared_ptr<PreProcessorStrategy> strategy);

	void init(std::shared_ptr<Parameter> para, int level);
	bool checkParameter();

private:
	PreProcessorImp();

	std::vector< std::shared_ptr<PreProcessorStrategy>> strategies;
	
};
#endif