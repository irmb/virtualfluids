#ifndef INIT_COMP_SP27_H
#define INIT_COMP_SP27_H

#include "PreProcessor\PreProcessorStrategy\PreProcessorStrategy.h"

#include <memory>

class Parameter;

class InitCompSP27 : public PreProcessorStrategy
{
public:
	static std::shared_ptr<PreProcessorStrategy> getNewInstance(std::shared_ptr< Parameter> para);
	void init(int level);
	bool checkParameter();

private:
	InitCompSP27();
	InitCompSP27(std::shared_ptr< Parameter> para);
	std::shared_ptr< Parameter> para;
};

#endif 