#include "PreProcessorImp.h"

#include "PreProcessorStrategy/PreProcessorStrategy.h"

#include "Parameter/Parameter.h"

std::shared_ptr<PreProcessorImp> PreProcessorImp::getNewInstance()
{
	return std::shared_ptr<PreProcessorImp>(new PreProcessorImp());
}

void PreProcessorImp::addStrategy(std::shared_ptr<PreProcessorStrategy> strategy)
{
	strategies.push_back(strategy);
}

void PreProcessorImp::init(std::shared_ptr<Parameter> para, int level)
{
	para->getParD(level)->evenOrOdd = false;
	for (std::size_t i = 0; i < strategies.size(); i++)
		strategies.at(i)->init(level);

	para->getParD(level)->evenOrOdd = true;
	for (std::size_t i = 0; i < strategies.size(); i++)
		strategies.at(i)->init(level);
}

bool PreProcessorImp::checkParameter()
{
	for (std::size_t i = 0; i < strategies.size(); i++) {
		if (!strategies.at(i)->checkParameter())
			return false;
	}
	return true;
}

PreProcessorImp::PreProcessorImp()
{
	strategies.resize(0);
}
