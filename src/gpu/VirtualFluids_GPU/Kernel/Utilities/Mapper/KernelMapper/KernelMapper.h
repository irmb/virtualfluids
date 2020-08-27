#ifndef KERNEL_MAPPPER_H
#define KERNEL_MAPPER_H

#include "Utilities/EnumMapper/EnumMapperImp.h"
#include "Kernel//Utilities/KernelType.h"

#include "VirtualFluids_GPU_export.h"

#include <memory>
#include <string>

class VIRTUALFLUIDS_GPU_EXPORT KernelMapper
{
public:
	static std::shared_ptr<KernelMapper> getInstance(); 

	std::string getString(KernelType enumeration);
	KernelType getEnum(std::string name);

private:
	KernelMapper();

	EnumMapperImp<KernelType> myEnumMapper;
};
#endif