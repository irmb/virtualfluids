#include "invalidInput_error.h"

#include <iostream>


invalidInput_error::invalidInput_error(char const * const message) throw() : runtime_error(message)
{
}

char const * invalidInput_error::what() const throw()
{
	return runtime_error::what();
}

