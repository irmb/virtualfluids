#include "RandomHelper.h"



mt19937 RandomHelper::make_engine()
{
	std::random_device r;
	std::seed_seq seed{r()};
	return std::mt19937(seed);
}