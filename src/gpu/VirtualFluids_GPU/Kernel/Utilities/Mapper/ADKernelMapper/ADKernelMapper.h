#ifndef AD_KERNEL_MAPPPER_H
#define AD_KERNEL_MAPPER_H

#include "Utilities/EnumMapper/EnumMapperImp.h"
#include "Kernel//Utilities/KernelType.h"

#include "VirtualFluids_GPU_export.h"

#include <memory>
#include <string>

class VIRTUALFLUIDS_GPU_EXPORT ADKernelMapper
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