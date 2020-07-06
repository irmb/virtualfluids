#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <memory>

class Parameter;

class PreProcessor
{
public:
	virtual void init(std::shared_ptr<Parameter> para, int level) = 0;
	virtual bool checkParameter() = 0;
	
};
#endif