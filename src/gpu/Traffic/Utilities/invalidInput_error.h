#pragma once
#include <VirtualFluidsDefinitions.h>

// using standard exceptions
#include <stdexcept>

class VIRTUALFLUIDS_GPU_EXPORT invalidInput_error :
	public std::runtime_error
{
public:
	invalidInput_error(char const* const message) throw();
	virtual char const* what() const throw();
};

