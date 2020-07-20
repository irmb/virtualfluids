#pragma once

#include <vector>

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

class VIRTUALFLUIDS_GPU_EXPORT VectorHelper
{
public:
	static void fillVector(std::vector<int>& vector, int insertNumber);
	static void fillVector(std::vector<std::vector<int> > &vector, int insertNumber);

	static void dispVector(const std::vector<int> &vector);
	static void dispVector(const std::vector<std::vector<int> >& vector);
	static void dispVectorColour(const std::vector<int> &vector);
	static void dispVectorColour(const std::vector<std::vector<int> >& vector);
	static void makeVectorOutputColourful(const int outputNumber);
};

