#pragma once

#include <vector>


#include "Core/DataTypes.h"

#include "Traffic_export.h"

class TRAFFIC_EXPORT VectorHelper
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

