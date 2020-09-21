#pragma once


// using standard exceptions
#include <stdexcept>

#include "Traffic_export.h"

class TRAFFIC_EXPORT invalidInput_error :
	public std::runtime_error
{
public:
	invalidInput_error(char const* const message) throw();
	virtual char const* what() const throw();
};

