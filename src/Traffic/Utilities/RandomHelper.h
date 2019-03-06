#pragma once
#include <random>
#include <VirtualFluidsDefinitions.h>

class VF_PUBLIC RandomHelper
{
public:
	static std::mt19937 make_engine();
};

