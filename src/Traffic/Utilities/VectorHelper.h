#pragma once
#include <vector>
#include <iostream>

#include <windows.h> //for colourful console output
#include <iomanip>	//formatting output streams

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

class VF_PUBLIC VectorHelper
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

