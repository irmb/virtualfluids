#include "ADKernelMapper.h"

std::shared_ptr<ADKernelMapper> ADKernelMapper::getInstance()
{
	static std::shared_ptr<ADKernelMapper> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<ADKernelMapper>(new ADKernelMapper());
	return uniqueInstance;
}

std::string ADKernelMapper::getString(ADKernelType enumeration)
{
	return myEnumMapper.getString(enumeration);
}

ADKernelType ADKernelMapper::getEnum(std::string name)
{
	return myEnumMapper.getEnum(name);
}

ADKernelMapper::ADKernelMapper()
{
	myEnumMapper.addEnum(LB_ADComp27, "ADComp27");
	myEnumMapper.addEnum(LB_ADComp7, "ADComp7");
	myEnumMapper.addEnum(LB_ADIncomp27, "ADIncomp27");
	myEnumMapper.addEnum(LB_ADIncomp7, "ADIncomp7");
}