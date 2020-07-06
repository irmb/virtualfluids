#ifndef INITF3_H
#define INITF3_H

#include "PreProcessor/PreProcessorStrategy/PreProcessorStrategy.h"

#include <memory>

class Parameter;

class InitF3 : public PreProcessorStrategy
{
public:
	static std::shared_ptr<PreProcessorStrategy> getNewInstance(std::shared_ptr< Parameter> para);
	void init(int level);
	bool checkParameter();

private:
	InitF3();
	InitF3(std::shared_ptr< Parameter> para);
	std::shared_ptr< Parameter> para;
};

#endif 