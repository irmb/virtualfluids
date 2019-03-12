#pragma once
#include <VirtualFluidsDefinitions.h>

// using standard exceptions
#include <iostream>
#include <stdexcept>

class VF_PUBLIC invalidInput_error :
	public std::runtime_error
{
public:
	invalidInput_error(char const* const message) throw();
	virtual char const* what() const throw();
};

