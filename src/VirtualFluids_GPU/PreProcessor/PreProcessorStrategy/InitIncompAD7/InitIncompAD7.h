#ifndef INIT_INCOMP_AD7_H
#define INIT_INCOMP_AD7_H

#include "PreProcessor/PreProcessorStrategy/PreProcessorStrategy.h"

#include <memory>

class Parameter;

class InitIncompAD7 : public PreProcessorStrategy
{
public:
	static std::shared_ptr<PreProcessorStrategy> getNewInstance(std::shared_ptr< Parameter> para);
	void init(int level);
	bool checkParameter();

private:
	InitIncompAD7();
	InitIncompAD7(std::shared_ptr< Parameter> para);
	std::shared_ptr< Parameter> para;
};

#endif 