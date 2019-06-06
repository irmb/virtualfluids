#ifndef AD_KERNEL_MAPPPER_H
#define AD_KERNEL_MAPPER_H

#include "Utilities/EnumMapper/EnumMapperImp.h"
#include "Kernel//Utilities/KernelType.h"

#include <memory>

class VF_PUBLIC ADKernelMapper
{
public:
	static std::shared_ptr<ADKernelMapper> getInstance();

	std::string getString(ADKernelType enumeration);
	ADKernelType getEnum(std::string name);

private:
	ADKernelMapper();

	EnumMapperImp<ADKernelType> myEnumMapper;
};
#endif