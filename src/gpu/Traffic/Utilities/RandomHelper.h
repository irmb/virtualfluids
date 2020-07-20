#pragma once

#include <random>
#include <VirtualFluidsDefinitions.h>

class VIRTUALFLUIDS_GPU_EXPORT RandomHelper
{
public:
	static std::mt19937 make_engine();
};

