#pragma once

#include <random>

#include "Traffic_export.h"

class TRAFFIC_EXPORT RandomHelper
{
public:
	static std::mt19937 make_engine();
};

