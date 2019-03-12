#pragma once

#include <stdexcept>
#include "Core/DataTypes.h"


static uint castSizeT_Uint(size_t number) {
	if (number > UINT_MAX)
	{
		throw std::overflow_error("number is larger than UINT_MAX");
	}
	return static_cast<uint>(number);
}

static int castSizeT_Int(size_t number) {
	if (number > INT_MAX)
	{
		throw std::overflow_error("number is larger than INT_MAX");
	}
	return static_cast<uint>(number);
}